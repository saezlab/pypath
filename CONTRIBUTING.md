# Contribution guidelines

Here are some guidelines to affectively contribute to PyPath development.
Please make sure you follow them as much as possible.

## Adding a new resource

- Add all the required information in `descriptions.py`
- Make sure which type of licensing the resource has.
- Double check the formating of the data upon import - both in PyPath itself
as in the export table (e.g. web service).
- Rebuild the webpage HTML.

## Adding a new class, method or function

- Write the docstrings - description of the object, inputs, outputs, examples
...
- Rebuild the documentation.
- Include in `__all__` module variable if applicable.
- Write the unit test(s) necessary.
- Add dependencies to `setup.py` if applicable.
