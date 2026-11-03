"""CLI for omnipath-utils."""

from __future__ import annotations

import sys


def main():
    """Main CLI entry point."""
    if len(sys.argv) < 2:
        print('Usage: omnipath-utils <command> [options]')
        print('Commands: build, serve, export-resolver')
        sys.exit(1)

    cmd = sys.argv[1]

    if cmd == 'build':
        from omnipath_utils.cli._build import build_cmd

        build_cmd(sys.argv[2:])
    elif cmd == 'serve':
        from omnipath_utils.cli._serve import serve_cmd

        serve_cmd(sys.argv[2:])
    elif cmd == 'export-resolver':
        from omnipath_utils.cli._export_resolver import export_resolver_cmd

        export_resolver_cmd(sys.argv[2:])
    else:
        print(f'Unknown command: {cmd}')
        sys.exit(1)


if __name__ == '__main__':
    main()
