<div id="main" class="col-md-9" role="main">

# load_images: Place tissue images within STlist

<div class="ref-description section level2">

Loads the images from tissues to the appropriate STlist slot.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
load_images(x = NULL, images = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist

-   images:

    a string indicating a folder to load images from

</div>

<div class="section level2">

## Value

an STlist with images

</div>

<div class="section level2">

## Details

This function looks for `.PNG` or `.JPG` files within a folder matching
the sample names in an existing STlist. Then, loads the images to the
STlist which can be used for plotting along with other spatialGE plots.

</div>

</div>
