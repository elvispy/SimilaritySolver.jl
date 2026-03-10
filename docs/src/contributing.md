# Contributing to SimilaritySolver.jl

We welcome contributions from the scientific computing community! Whether you're reporting a bug, proposing a new feature, or improving the documentation, your help is appreciated.

## Reporting Bugs

If you encounter an issue, please open an issue on GitHub. To help us resolve it quickly, please include:
- A **Minimal Working Example (MWE)** that reproduces the error.
- Your Julia version and package versions (`Pkg.status()`).
- A brief description of the expected vs. actual behavior.

## Proposing Features

We are particularly interested in features that expand the scope of automated symmetry analysis. If you have a proposal:
1. Open an issue to discuss the feature's scope and design.
2. Provide a motivating example (e.g., a PDE that the current dilation method cannot solve).

## Development Setup

To contribute code:

1. **Fork the repository** on GitHub.
2. **Clone your fork** and set up the development environment:
   ```bash
   git clone https://github.com/YOUR_USERNAME/SimilaritySolver.jl.git
   cd SimilaritySolver.jl
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```
3. **Run the tests** to ensure everything is working:
   ```bash
   julia --project test/runtests.jl
   ```

## Contribution Workflow

1. Create a new branch for your feature or fix: `git checkout -b feature/your-feature-name`.
2. Implement your changes and **add tests** for any new functionality.
3. Ensure the full test suite passes.
4. Open a Pull Request (PR) with a clear description of your changes.

## Documentation Improvements

Documentation is built with `Documenter.jl`. To build and preview changes locally:
```bash
julia --project=docs docs/make.jl
```
Open `docs/build/index.html` in your browser.

## Community Standards

We follow the [Julia Community Standards](https://julialang.org/community/standards/). Please be respectful and professional in all interactions.
