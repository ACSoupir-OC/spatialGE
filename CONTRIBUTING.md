# Contributing to spatialGE

Thank you for your interest in contributing to spatialGE! This document provides guidelines for contributing to the package.

## Types of Contributions

### 1. Report Bugs

Report bugs at https://github.com/FridleyLab/spatialGE/issues

When filing an issue, include:
- Your operating system and R version
- Steps to reproduce the bug
- Expected vs actual behavior
- A minimal reproducible example if possible

### 2. Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

### 3. Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

### 4. Write Documentation

You can never have enough documentation! Contributing to documentation includes:
- Improving function examples in Roxygen comments
- Writing vignettes for new workflows
- Fixing typos in documentation
- Adding FAQ entries

### 5. Submit Feedback

The best way to send feedback is to file an issue at https://github.com/FridleyLab/spatialGE/issues

If you are proposing a feature:
- Explain in detail how it would work
- Keep the scope as narrow as possible
- Remember that this is a volunteer-driven project

## Get Started!

Ready to contribute? Here's how to set up spatialGE for local development:

1. Fork the spatialGE repo on GitHub
2. Clone your fork locally:
   ```bash
   git clone git@github.com:your_name_here/spatialGE.git
   cd spatialGE
   ```

3. Install dependencies:
   ```r
   install.packages("devtools")
   devtools::install_deps()
   ```

4. Create a branch for local development:
   ```bash
   git checkout -b name-of-your-bugfix-or-feature
   ```

5. Make your changes following the guidelines below

6. Test your changes:
   ```r
   devtools::test()
   devtools::check()
   ```

7. Commit your changes and push your branch to GitHub:
   ```bash
   git add .
   git commit -m "Your detailed description of your changes."
   git push origin name-of-your-bugfix-or-feature
   ```

8. Submit a pull request through the GitHub website

## Development Guidelines

### Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use spaces, not tabs
- Keep lines under 80 characters
- Use meaningful variable names

### Documentation

- All exported functions must have Roxygen documentation
- Include @examples for every function
- Examples should be <20 lines and use test data
- Use `\dontrun{}` for examples that require internet or long runtime

### Testing

- Write unit tests for new functions using testthat
- Aim for >80% code coverage
- All tests must pass before submitting a PR
- Run `devtools::test()` before committing

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):
- `feat:` - New features
- `fix:` - Bug fixes
- `docs:` - Documentation changes
- `refactor:` - Code refactoring (no behavior change)
- `test:` - Test additions or modifications
- `chore:` - Maintenance tasks

Example:
```
feat: add spatial gradient detection to STgradient

- Added Spearman correlation for gradient detection
- Added log-distance transformation option
- Updated documentation with examples
```

## Pull Request Process

1. Update the NEWS.md file with your changes
2. Update documentation with `devtools::document()`
3. Ensure all tests pass
4. Request review from maintainers
5. Address review comments
6. Merge after approval

## Code of Conduct

Please note that this project is released with a [Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Questions?

Contact the maintainers:
- Oscar Ospina: oscar.ospina@jhmi.edu
- Alex Soupir: Alex.Soupir@moffitt.org

Thank you for contributing to spatialGE! 🎉
