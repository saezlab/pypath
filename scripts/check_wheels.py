#!/usr/bin/env python3
"""Build and test each package with only its declared, locked dependencies.

No resource build or database connection is made. Wheels and disposable virtual
environments live in a temporary directory; use --keep to retain them for review.
"""

import argparse
from pathlib import Path
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("core", "resolver", "build", "subsets", "postgres", "api", "client")
DEPENDENCIES = {
    "core": (),
    "resolver": ("core",),
    "build": ("core", "resolver"),
    "subsets": ("core",),
    "postgres": ("core", "subsets"),
    "api": ("core",),
    # The client reads its SQL and contract from its own package; nothing of the workspace.
    "client": (),
}
# Package directories under packages/, where not the package's name.
DIRECTORIES = {"client": "query/python"}


def run(*command, cwd=ROOT):
    subprocess.run(command, cwd=cwd, check=True)


def check(directory, packages, *, offline=False):
    wheels = directory / "wheels"
    wheels.mkdir()
    offline_args = ("--offline",) if offline else ()
    required = set(packages)
    for package in packages:
        required.update(DEPENDENCIES[package])

    def build_wheel(source_args, name):
        # Build the wheel from a fresh sdist so stale build/lib files in an
        # editable checkout cannot enter a released wheel.
        run(
            "uv",
            "build",
            *offline_args,
            *source_args,
            "--sdist",
            "--build-constraints",
            "deploy/build-constraints.txt",
            "--out-dir",
            str(wheels),
        )
        archive = next(wheels.glob(f"{name}-*.tar.gz"))
        run(
            "uv",
            "build",
            *offline_args,
            str(archive),
            "--wheel",
            "--build-constraints",
            "deploy/build-constraints.txt",
            "--out-dir",
            str(wheels),
        )

    for package in PACKAGES:
        if package in required:
            build_wheel(("--package", f"omnipath-{package}"), f"omnipath_{package}")
    if required & {"build", "postgres"}:
        build_wheel(("pypath",), "pypath_omnipath")
    for package in required:
        wheel = next(wheels.glob(f"omnipath_{package}-*.whl"))
        with zipfile.ZipFile(wheel) as archive:
            names = set(archive.namelist())
            assert not any(
                name.startswith(("tests/", "scripts/", "src/", "rust/")) for name in names
            )
            source = ROOT / "packages" / DIRECTORIES.get(package, package) / f"omnipath_{package}"
            expected = {
                f"omnipath_{package}/{path.relative_to(source).as_posix()}": path
                for path in source.rglob("*")
                if path.is_file() and path.suffix in {".py", ".sql", ".yaml", ".zstd"}
                or (package == "client" and path.suffix == ".json")
            }
            actual = {
                name
                for name in names
                if name.startswith(f"omnipath_{package}/")
                and (
                    Path(name).suffix in {".py", ".sql", ".yaml", ".zstd"}
                    or (package == "client" and Path(name).suffix == ".json")
                )
            }
            assert actual == set(expected), {
                "extra": actual - set(expected),
                "missing": set(expected) - actual,
            }
            for name, path in expected.items():
                assert archive.read(name) == path.read_bytes(), name
    for package in packages:
        environment = directory / package
        run("uv", "venv", str(environment))
        python = environment / "bin/python"
        requirements = directory / f"{package}-requirements.txt"
        run(
            "uv",
            "export",
            "--frozen",
            "--package",
            f"omnipath-{package}",
            "--no-dev",
            "--no-emit-workspace",
            "--no-emit-package",
            "pypath-omnipath",
            "--no-hashes",
            "--output-file",
            str(requirements),
        )
        install = [
            str(next(wheels.glob(f"omnipath_{name}-*.whl")))
            for name in (*DEPENDENCIES[package], package)
        ]
        if package in {"build", "postgres"}:
            install.append(str(next(wheels.glob("pypath_omnipath-*.whl"))))
        run(
            "uv",
            "pip",
            "install",
            *offline_args,
            "--python",
            str(python),
            "--no-deps",
            "--build-constraints",
            str(ROOT / "deploy/build-constraints.txt"),
            "-r",
            str(requirements),
            *install,
        )
        # -I removes both cwd and PYTHONPATH; subprocess cwd is outside the repo.
        code = f"import omnipath_{package}; import importlib.util\n"
        forbidden = {
            "core": ("omnipath_build", "omnipath_postgres", "omnipath_api"),
            "resolver": ("omnipath_build", "pypath", "duckdb"),
            "subsets": (
                "omnipath_postgres",
                "omnipath_build",
                "pypath",
                "psycopg",
                "omnipath_subsets.compatibility",
            ),
            "postgres": ("omnipath_build", "omnipath_resolver", "psycopg"),
            "api": ("omnipath_build", "omnipath_resolver", "pypath"),
            "client": (
                "omnipath_core",
                "omnipath_build",
                "omnipath_resolver",
                "pypath",
                "omnipath_api",
                "fastapi",
                "pyarrow",
                "biolink_model",
            ),
        }.get(package, ())
        for name in forbidden:
            code += f"assert importlib.util.find_spec({name!r}) is None, {name!r}\n"
        if package in {"build", "subsets", "postgres", "api", "client"}:
            code += f"from omnipath_{package}.cli import main\n"
            run(
                str(python),
                "-I",
                "-c",
                f"from omnipath_{package}.cli import main; main(['--help'])",
                cwd=directory,
            )
        if package == "build":
            code += (
                "from omnipath_build.provenance import runtime_provenance\n"
                "software = runtime_provenance()['software']['packages']\n"
                "assert set(software) == {'omnipath-core', 'omnipath-build', 'omnipath-resolver', 'pypath-omnipath', 'biolink-model'}\n"
                "assert software['omnipath-resolver']['native_files']\n"
                "assert all(software[name]['vocabulary']['content_sha256'] for name in ('omnipath-core', 'pypath-omnipath', 'biolink-model'))\n"
            )
        if package == "resolver":
            code += (
                "from omnipath_resolver import EntityResolver, LibraryMatcher, IdentityRuntime\n"
            )
        if package == "subsets":
            code += "from omnipath_subsets.scientific import run_product\n"
        run(str(python), "-I", "-c", code, cwd=directory)
        print(f"PASS isolated {package} wheel", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packages", nargs="+", choices=PACKAGES, default=PACKAGES)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--keep", type=Path)
    args = parser.parse_args()
    if args.keep:
        args.keep.mkdir(parents=True, exist_ok=False)
        check(args.keep.resolve(), args.packages, offline=args.offline)
    else:
        with tempfile.TemporaryDirectory(prefix="omnipath-wheel-check-") as directory:
            check(Path(directory), args.packages, offline=args.offline)


if __name__ == "__main__":
    main()
