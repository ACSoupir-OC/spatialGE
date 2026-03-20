# spatialGE Development Guidelines

## Version Control with Git

### Commit Message Convention

Use descriptive commit messages following the format:

```
[area] brief description of change

More detailed explanation if needed:
- What changed
- Why it changed
- Any implications
```

**Areas:** `core`, `plotting`, `tests`, `docs`, `utils`, `build`

**Examples:**
- `core: fix SThet Moran calculation for sparse matrices`
- `plotting: update ggplot2 compatibility for visualization functions`
- `tests: add test coverage for STgradient edge cases`

### Before Making Changes

1. **Check current branch status**: `git status`
2. **Pull latest changes**: `git pull`
3. **Create feature branch**: `git checkout -b fix-ggplot2-compatibility`
4. **Make incremental commits** with meaningful messages
5. **Test before committing** to ensure nothing breaks

### When Committing

- **Always test** before committing
- **Describe intent** clearly in commit message
- **Include test results** if relevant (pass/fail counts)
- **Reference related issues** if applicable

## Development Workflow

1. Clone repo or pull latest changes
2. Create feature branch for your work
3. Make changes incrementally
4. Run tests after each meaningful change
5. Commit with descriptive messages
6. Push and create PR for review

## Testing Before Commit

Run the test suite:
```bash
Rscript run_tests.R
```

Document results in commit message if significant changes.

## Key Files

- `R/` - Core R code (main refactoring target)
- `src/` - C++ code (Rcpp)
- `tests/testthat/` - Test files
- `RefactoringPlan.md` - Refactoring strategy
- `DEVELOPMENT.md` - This file

---

*Last updated: 2026-02-27*