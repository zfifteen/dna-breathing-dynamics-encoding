# Makefile for DNA Breathing Dynamics Framework
# Supports Apple Silicon optimization and comprehensive testing

.PHONY: help deps venv install build test clean format lint typecheck smoke bench all

# Default target
.DEFAULT_GOAL := help

# Detect platform
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Python configuration
PYTHON := python3
VENV := .venv
VENV_BIN := $(VENV)/bin
PIP := $(VENV_BIN)/pip

# Homebrew prefix detection (Apple Silicon vs Intel)
ifeq ($(UNAME_M),arm64)
    HOMEBREW_PREFIX := /opt/homebrew
    PLATFORM_FLAGS := -DAPPLE_SILICON
else
    HOMEBREW_PREFIX := /usr/local
    PLATFORM_FLAGS :=
endif

# Build directories
BUILD_DIR := build
DIST_DIR := dist
DOCS_BUILD := docs/_build

help:  ## Show this help message
	@echo "DNA Breathing Dynamics Framework - Build System"
	@echo ""
	@echo "Platform: $(UNAME_S) $(UNAME_M)"
	@echo "Homebrew: $(HOMEBREW_PREFIX)"
	@echo ""
	@echo "Available targets:"
	@awk 'BEGIN {FS = ":.*##"; printf "\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Dependencies

deps:  ## Install system dependencies (macOS Homebrew)
ifeq ($(UNAME_S),Darwin)
	@echo "Installing Homebrew dependencies..."
	brew install mpfr gmp libomp
	@echo "System dependencies installed."
else
	@echo "Non-macOS system detected. Please install mpfr, gmp manually."
endif

venv:  ## Create Python virtual environment
	@echo "Creating virtual environment..."
	$(PYTHON) -m venv $(VENV)
	@echo "Virtual environment created at $(VENV)"

install: venv  ## Install Python dependencies
	@echo "Installing Python dependencies..."
	$(PIP) install --upgrade pip setuptools wheel
	$(PIP) install -r requirements.txt
	$(PIP) install -e .
	@echo "Dependencies installed."

##@ Build

build: install  ## Build C/C++ extensions
	@echo "Building C/C++ extensions..."
	$(VENV_BIN)/python setup.py build_ext --inplace
	@echo "Build complete."

##@ Testing

test: build  ## Run all tests
	@echo "Running test suite..."
	$(VENV_BIN)/pytest tests/ -v

smoke:  ## Run smoke tests (<5s quick validation)
	@echo "Running smoke tests..."
	$(VENV_BIN)/pytest tests/ -v -m smoke --tb=short

unit:  ## Run unit tests only
	@echo "Running unit tests..."
	$(VENV_BIN)/pytest tests/unit/ -v

integration:  ## Run integration tests
	@echo "Running integration tests..."
	$(VENV_BIN)/pytest tests/integration/ -v

performance:  ## Run performance tests
	@echo "Running performance tests..."
	$(VENV_BIN)/pytest tests/performance/ -v -m performance

validation:  ## Run scientific validation tests
	@echo "Running validation tests..."
	$(VENV_BIN)/pytest tests/ -v -m validation

bench:  ## Run performance benchmarks
	@echo "Running benchmarks..."
	$(VENV_BIN)/pytest tests/performance/ --benchmark-only --benchmark-autosave

##@ Code Quality

format:  ## Format code with Black and isort
	@echo "Formatting code..."
	$(VENV_BIN)/black src/ tests/ examples/
	$(VENV_BIN)/isort src/ tests/ examples/
	@echo "Code formatted."

lint:  ## Lint code with Flake8
	@echo "Linting code..."
	$(VENV_BIN)/flake8 src/ tests/ examples/ --max-line-length=88 --extend-ignore=E203,W503

typecheck:  ## Type check with MyPy
	@echo "Type checking..."
	$(VENV_BIN)/mypy src/

check: format lint typecheck  ## Run all code quality checks

##@ Cleanup

clean:  ## Remove build artifacts and caches
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR) $(DIST_DIR) $(DOCS_BUILD)
	rm -rf *.egg-info
	rm -rf .pytest_cache .coverage htmlcov .mypy_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.so" -delete
	find . -type f -name "*.c" ! -path "*/src/extensions/*" -delete  # Remove Cython-generated .c files
	@echo "Clean complete."

clean-venv:  ## Remove virtual environment
	@echo "Removing virtual environment..."
	rm -rf $(VENV)

clean-all: clean clean-venv  ## Remove all generated files including venv

##@ Development

dev: install  ## Set up development environment
	@echo "Development environment ready."
	@echo "Activate with: source $(VENV_BIN)/activate"

shell:  ## Start IPython shell with project loaded
	$(VENV_BIN)/ipython

jupyter:  ## Start Jupyter notebook
	$(VENV_BIN)/jupyter notebook

##@ Comprehensive Targets

all: deps install build test check  ## Full build and test cycle

ci: build test lint typecheck  ## CI/CD pipeline (without deps)

.PHONY: docs
docs:  ## Build documentation
	@echo "Building documentation..."
	@echo "Documentation build not yet configured."
