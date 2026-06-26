install:
	python -m pip install -e '.[dev]'

test:
	python -m pytest tests
