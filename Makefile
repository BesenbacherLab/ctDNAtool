.PHONY: update-deps init update test

test:
	poetry run pytest

lint:
	poetry run flake8 src/ctDNAtool test/

black:
	poetry run black src/ctDNAtool test/

clean:
	find . -prune -name ".egg-info" -type d -exec rm -rf {} ';'
	find . -prune -name ".eggs" -type d -exec rm -rf {} ';'
	find . -prune -name "__pycache__" -type d -exec rm -rf {} ';'
	rm -rf build/ dist/ .pytest_cache .egg .nox
