# Support and Help

## Documentation

- **Package vignettes**: https://fridleylab.github.io/spatialGE/articles/
- **Reference manual**: https://fridleylab.github.io/spatialGE/reference/
- **Quick start guide**: `vignette("getting_started", package="spatialGE")`

## Getting Help

### 1. Check Documentation

Before filing an issue, check:
- Function documentation: `?function_name` or `help(function_name)`
- Vignettes: `vignette(package="spatialGE")`
- FAQ in README.md

### 2. GitHub Issues

File issues at: https://github.com/FridleyLab/spatialGE/issues

**When filing an issue, include:**
- R version (`R.version.string`)
- spatialGE version (`packageVersion("spatialGE")`)
- Operating system
- Minimal reproducible example
- Expected vs actual behavior

### 3. Email

For questions not suitable for public issues, contact:
- Oscar Ospina: oscar.ospina@jhmi.edu
- Alex Soupir: Alex.Soupir@moffitt.org

### 4. Bioconductor Support

For general Bioconductor questions:
- Bioconductor Support Site: https://support.bioconductor.org/
- Tag posts with `spatialGE`

## Reporting Bugs

Bugs should be reported at: https://github.com/FridleyLab/spatialGE/issues

**Bug report checklist:**
- [ ] Confirmed it's not user error
- [ ] Checked existing issues
- [ ] Created minimal reproducible example
- [ ] Included session info (`sessionInfo()`)

## Feature Requests

Feature requests are welcome! File an issue with:
- Clear description of the feature
- Use case / motivation
- Example of how it would work
- References to similar tools (if applicable)

## Training Materials

### Vignettes

1. **Getting Started**: Basic workflow and STlist creation
2. **STlist Object Structure**: Understanding the data structure
3. **Spatial Heterogeneity**: Using SThet for Moran's I analysis
4. **Spatial Clustering**: Using STclust for tissue domain detection
5. **Differential Expression**: Using STdiff for spatial DE analysis
6. **Gene Set Enrichment**: Using STenrich for pathway analysis
7. **Spatial Gradients**: Using STgradient for gradient detection

### External Resources

- Spatial transcriptomics review: [DOI:10.1093/bioinformatics/btac145](https://doi.org/10.1093/bioinformatics/btac145)
- 10X Visium tutorials: https://support.10xgenomics.com/
- Seurat spatial vignettes: https://satijalab.org/seurat/articles/

## Version Support

| Version | Status | Support |
|---------|--------|---------|
| 2.0.x   | Current | Full support |
| 1.2.x   | Legacy | Bug fixes only |
| <1.2    | Deprecated | No support |

## Response Time

We aim to respond to:
- Bug reports: Within 1 week
- Feature requests: Within 2 weeks
- General questions: Within 1 week

## Acknowledgments

When using spatialGE in your research, please cite:

> Ospina O, Soupir A, Fridley B (2024). spatialGE: An R package for visualization and analysis of spatial heterogeneity in spatially-resolved gene expression. Bioinformatics, 40(3), btae123. doi:10.1093/bioinformatics/btae123

---

**Thank you for using spatialGE!**
