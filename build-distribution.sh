# After updates to the repository, run these commands to deploy Persimmon to Pypi

pip install -U twine # MUST UPDATE - OLD VERSION DOESN'T WORK
python setup.py sdist
pip install wheel
python setup.py bdist_wheel

# twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# twine upload dist/*