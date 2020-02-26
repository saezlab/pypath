# Contribution guidelines

Here are some guidelines to effectively contribute to PyPath development.
Please make sure you follow them as much as possible.

## Adding a new resource

- Add all the required information in `descriptions.py`
- Make sure which type of licensing the resource has.
- Double check the formatting of the data upon import - both in PyPath itself
as in the export table (e.g. web service).
- Rebuild the webpage HTML.

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
has undergone large refactoring efforts over the past year (2019) and it has
still around 30% of its working codebase left with a sometimes messy and ugly
style. Please be forgiving when you see these parts and feel free to improve
them if you want. To see a more consistent style characteristic to this module
we recommend to look at some new core modules, e.g. `pypath.core.network` or
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
in a new line. The closing parenthese returns to the upper indentation level
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
parenthese has the original indentation level in a new line.

```python

foo = [
    bar
    for bar in baz
    if bar[0]
]

```

### Docstrings

After the triple double quotes the docstring starts in a new line and also
the closing quotes have their own line. Arguments are described as
`:param type,type name:` and the description follows in a new line, indented.

```python

def function(a):
    """
    Does something.
    
    :param int,float a:
        A number.
    
    :returns:
        A list of integers.
    """
    
    return [a, a + 1, a + 2]

```

### Spaces around operators

Operators never stick to the operands, always surrounded by spaces.
Also there is a space after the comma when listing elements of tuples and
lists within a line.

```python

def function(a = 1):
    
    a = a * 4 + 3

```

### Quotes

We use single quotes unless there is any strong reason to use double.


## Logging

`pypath` has its own session and log manager built in. A class either should
inherit from `pypath.share.session.Logger` and then have its `_log` method
or a logger instance should be created at the module level and used by
methods and classes in that module. Each logger instance have a default label
so in the log one can see which message comes from which part of the module.


[1] https://www.python.org/dev/peps/pep-0008/
[2] https://github.com/google/styleguide/blob/gh-pages/pyguide.md
