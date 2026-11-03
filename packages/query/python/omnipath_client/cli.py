"""Command line: install resources into the local data directory, or list them."""

from __future__ import annotations

import argparse

from .installation import DIMENSIONS, available, data_directory, init, read_installation


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="omnipath-client", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    def common(command):
        command.add_argument(
            "--release", default="latest", help="OmniPath release (default: latest)"
        )
        command.add_argument(
            "--data-dir", help="data directory (default: OMNIPATH_DATA_DIR or per user)"
        )
        command.add_argument("--api-url", help="API to read the catalog from")

    install = commands.add_parser(
        "init",
        help="download resources and prepare them for querying",
        description="Download resources and prepare them for querying. Without a selection, "
        "every resource of the release. Running it again updates what is installed.",
    )
    install.add_argument("resources", nargs="*", help="resource names (default: all)")
    for dimension in DIMENSIONS:
        install.add_argument(
            f"--{dimension.replace('_', '-')}",
            nargs="+",
            metavar="TAG",
            help=f"only resources with one of these {dimension.replace('_', ' ')} tags",
        )
    install.add_argument(
        "--evidence", action="store_true", help="also download raw evidence payloads"
    )
    install.add_argument(
        "--materialize",
        action=argparse.BooleanOptionalAction,
        help="native DuckDB tables (the default for a new directory; about twice the files' "
        "size in addition) or, with --no-materialize, views over the files (small, but "
        "the prebuilt queries slow down with many resources); kept on later runs",
    )
    common(install)
    listing = commands.add_parser(
        "list", help="list the resources of a release, their tags, sizes and what is installed"
    )
    common(listing)
    args = parser.parse_args(argv)

    if args.command == "init":
        init(
            args.resources or None,
            **{dimension: getattr(args, dimension) for dimension in DIMENSIONS},
            release=args.release,
            data_dir=args.data_dir,
            include_evidence=args.evidence,
            materialize=args.materialize,
            api_url=args.api_url,
            progress=lambda message: print(message, flush=True),
        )
        return
    installation = read_installation(data_directory(args.data_dir))
    installed = {
        r["resource_id"]: r["version"] for r in (installation or {"resources": []})["resources"]
    }
    for record in available(release=args.release, api_url=args.api_url):
        name = record["resource_id"]
        size = sum(f["size_bytes"] for f in record["files"]) / 1e6
        tags = " ".join(sorted(tag["id"] for tag in record.get("tags") or [] if "id" in tag))
        state = (
            "installed"
            if installed.get(name) == record["version"]
            else (f"installed {installed[name]}" if name in installed else "")
        )
        print(f"{name:24s} {record['version']:12s} {size:8.1f} MB  {state:22s} {tags}")


if __name__ == "__main__":
    main()
