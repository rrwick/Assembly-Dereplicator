# Assembly dereplicator tests

Assembly dereplicator comes with a few automated tests to help with development and spotting bugs. You'll need [pytest](https://docs.pytest.org/en/latest/) installed to run them.

To run the tests, execute this command from Assembly dereplicator's root directory:
```
python3 -m pytest
```

Or if you have [Coverage.py](https://coverage.readthedocs.io/en/coverage-4.5.1a/) installed, you can run the tests through it to get some code coverage stats:
```
coverage run -m pytest && coverage report -m
```
