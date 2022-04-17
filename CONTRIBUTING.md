# Contribution guidelines

Here are some guidelines to effectively contribute to `pypath` development.
Please make sure you follow them as much as possible.

## Adding a new resource

- Create a new module in the [`pypath.inputs`][10] module, or a new function
  if a module for the resource already exists. Follow other input modules as
  examples and use the module infrastructure, especially [`share.curl`][7]
  for downloads and [`share.settings`][8] for configuration. Configuration
  is stored in [`settings.yaml`][9].
- Add all the required information in [`resources.json`][3] and
  [`descriptions.py`][4]. These information will be shown also on the
  OmniPath web page after the next update of the server.
- Most important to add the license, for license short names see this
  [directory][5]
- Check the loaded data directly in the output of the new input method
- Create input definitions for the relevant core databases (network,
  enzyme-substrate, complexes, annotations, intercell); this might be
  not straightforward, don't hesitate to ask for help
- Check the data loaded into the core databases
- Our [daily builds][6] might reveal issues and also here we can see if
  the new resource is loaded correctly into the web service data frames.

## Adding a new class, method or function

- Write the docstrings - description of the object, inputs, outputs, examples
...
- Rebuild the documentation.
- Include in `__all__` module variable if applicable.
- Write the unit test(s) necessary.
- Add dependencies to `setup.py` if applicable.

## Coding style

The `pypath` aims to follow the style standards most widely shared in the
Python community such as [PEP8][1] and the [Google style guide][2]. `pypath`
has undergone large refactoring efforts since 2019 and it has still around 30%
of its working codebase left with a sometimes messy and ugly style. Please
forgiving when you see these parts and feel free to improvethem if you want.
To see a more consistent style characteristic to the module, we recommend to
look at some new core modules, e.g. `pypath.core.network` or
`pypath.core.interaction`. As we mentioned we follow the standards, we have
exceptions in a few points. In the document below we only discuss these
alterations and additional rules. For general guidelines you can refer to the
[PEP8][1] and the [Google style guide][2].

### Blank lines

Vertical empty space does not cost anything and helps a lot with readability.
We use it extensively. Between, before and after methods and classes we leave
2 empty lines. Before and after blocks we leave an empty line. Inside blocks
we add blank lines between short segments of logically strongly connected
elements (usually one or a few lines).

```python

def function1(a = 1, b = 2):

    a = a + 1
    b = b - 1

    return a * b


def function2(a, b):

    for i in a:

        if i // b:

            yield i

```

### Indentation of argument lists

If the list of arguments or anything inside parentheses needs to be broken to
multiple lines we start newline after the opening parenthese and each element
should have its own line one level deeper indented than the opening line.
We add commas after each element including the last one, except if it is
`**kwargs` because comma after `**kwargs` is syntax error in Python 2.
If the elements are statements separated by operators (e.g. `+` or `and`)
we leave the operator at the end of the line and start the next statement
in a new line. The closing parenthesis returns to the upper indentation level
and left hanging in its own line.

```python

a = function(
    foo = 'foo',
    bar = 'bar',
    baz = 'baz',
)

a = function(
    foo = 'foo',
    bar = 'bar',
    **kwargs
)

if (
    condition1 and
    condition2
):

    do_something()

```

Anywhere else if the vertical alignment helps readability is good to use line
breaks. For example when similar statements are repeated.

```python

a = (
    record[i],
    record[i + 7],
)

```

### Long comprehensions

If the contents of a list or other comprehension does not fit into one line
or even if it fits but looks better if broken into multiple pieces, we start
the value, each `for` and `if` statements in a new line. The closing
parenthesis has the original indentation level in a new line.

```python

foo = [
    bar
    for bar in baz
    if bar[0]
]

```

### Docstrings

After the triple double quotes the docstring starts in a new line and also
the closing quotes have their own line. For docstrings we follow the
Napoleon (Google) standard:
https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

Thanks to type hints we don't have to add the type of arguments and return
values in parentheses (enough to write `arg:` instead of `arg (int):`).

```python

def function(a: int) -> List[int]:
    """
    Does something.

    Args:
        a: A number. I just make this sentence longer, hopefully
            enough long to show how to make a line break: indent the rest
            of the lines, so the argument name itself sticks out.

    Returns:
        A list of integers. I just make this sentence longer,
        hopefully enough long to show how to make a line break.
    """

    return [a, a + 1, a + 2]

```

### Python 3 features

By the end of 2021 we stopped supporting Python 2. Since then we use Python
3 only features such as [type hinting][11].

You can see many uses of `future.utils.iteritems()`, `past.builtins.xrange()`
instead of simply `.items()` and `range()`, or `map`, `filter`, `reduce`
(`past.builtins.reduce`) which all aimed to maintain Python 2 compatibility.
These are already not necessary and should not be used in new code.

### Spaces around operators

Operators never stick to the operands, always surrounded by spaces.
Also there is a space after the comma when listing elements of tuples and
lists within a line.

```python

def function(a = 1):

    a = a * 4 + 3
    a = (1, 2, 3)

```

### Quotes

We use single quotes unless there is any strong reason to use double.

### Imports

We group imports the following way:

- Compatibility imports (e.g. `future.utils.iteritems`; not used any more)
- Standard library imports
- Typing imports
- Imports from other third party dependencies (e.g. `numpy`)
- Imports from `pypath` itself

Between these sections we leave an empty line

## Logging

`pypath` has its own session and log manager built in. A class either should
inherit from `pypath.share.session.Logger` and then have its `_log` method
or a logger instance should be created at the module level and used by
methods and classes in that module. Each logger instance have a default label
so in the log one can see which message comes from which part of the module.


[1]: https://www.python.org/dev/peps/pep-0008/
[2]: https://github.com/google/styleguide/blob/gh-pages/pyguide.md
[3]: pypath/resources/data/resources.json
[4]: pypath/resources/descriptions.py
[5]: pypath/data/licenses
[6]: https://status.omnipathdb.org/inputs/
[7]: pypath/share/curl.py
[8]: pypath/share/settings.py
[9]: pypath/data/settings.yaml
[10]: pypath/inputs
[11]: https://realpython.com/lessons/type-hinting/
