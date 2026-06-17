install:
	pip install -e .[dev,notebook]

test:
	pytest
