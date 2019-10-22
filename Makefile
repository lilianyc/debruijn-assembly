default: help


env: ## Create conda env
	conda env create -f environment.yml
.PHONY: env

test: ## Run tests
	pytest tests
.PHONY: test

test-coverage: ## Run tests with coverage
	pytest --cov
.PHONY: test-coverage

lint: ## Lint code
	pylint debruijn
.PHONY: lint

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
