# Standard setup
source(test_path("helper-data.R"))
# Don't load spatialGE here - individual tests use devtools::load_all()
# to load the development version (prevents namespace conflicts)
