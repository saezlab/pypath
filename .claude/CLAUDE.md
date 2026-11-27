## Python Environment

This project uses `uv` for Python dependency management. When running Python commands, always use:

```bash
uv run python <script or -c "code">
```

Instead of calling `python` or `python3` directly.

## Python coding style

Python code follows the lab guide, `saezverse/guidelines/python-coding-style.md`
(in the `saezverse` repository). Read it before you write or change Python. The
main rules:

- Spaces around `=` in keyword arguments and defaults: `f(a = 1)`.
- Single quotes.
- Blank lines around blocks and after a `def` line.
- One argument per line, with a trailing comma, when a call breaks.
- Napoleon docstrings, with the quotes on their own lines.

`ruff format` does not follow the guide: it removes those spaces and switches the
quotes. Do not run it over code written in the lab style.
