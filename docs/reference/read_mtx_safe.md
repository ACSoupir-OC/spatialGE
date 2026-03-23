<div id="main" class="col-md-9" role="main">

# Safe Matrix Reader

<div class="ref-description section level2">

Robustly reads Matrix market files, handling .gz via temp expansion

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
read_mtx_safe(path)
```

</div>

</div>

<div class="section level2">

## Arguments

-   path:

    Path to matrix.mtx or matrix.mtx.gz file

</div>

<div class="section level2">

## Value

CsparseMatrix from Matrix package

</div>

<div class="section level2">

## Details

Decompresses gzip files to temporary location before reading, ensuring
proper handling of compressed Matrix Market format files.

</div>

</div>
