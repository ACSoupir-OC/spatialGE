##
# @title quilt_p
# @description Creates a quilt plot from ST data.
# @details
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values (expression or cell scores). The data
# frame has column names 'x_pos', 'y_pos', and 'values'. It also takes a color
# palette name from the 'khroma' or RColorBrewer packages. Finally, it takes a
# name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the values to be plotted.
# @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return, a ggplot object.
#
#' @import ggplot2
#
#
quilt_p <- function(data_f=NULL, color_pal="BuRd", leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T, ptsize=0.5){

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=xpos, y=ypos, color=values)) +
    geom_point(size=ptsize) +
    scale_color_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(color=leg_name, title=title_name) +
    theme_void() +
    theme(legend.position="right")
    #ggplot2::theme_classic()

  if(visium){
    p <- p + scale_y_reverse()
    # scale_y_reverse() + coord_fixed(ratio=1.7)
  }

  p <- p + coord_equal()

  return(p)
}


##
# @title krige_p
# @description Creates a kriging plot from ST data.
# @details
# Function to produce a "kriging plot" from a data frame with three columns:
# x coordinates, y coordinates, and predicted kriging values. The data frame has
# column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
# SpatialPolygons object to mask the predicted grid to the area of the tissue.
# It also takes a color palette name from the 'khroma' package. Finally, it
# takes a name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the kriging prediction values to be plotted.
# @param mask, an object of class SpatialPolygons containing a large polygon
# encasing all the predicted grid, and a smaller polygon drawing the concave hull
# of the tissue shape.
# @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @return, a ggplot object.
#
#' @import ggplot2
#' @import ggpolypath
#
#
krige_p <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', minvalue=minvalue, maxvalue=maxvalue,
                    visium=T){

  #requireNamespace('sf')

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Convert the SpatialPolygon mask into a data frame.
  #mask_df <- ggplot2::fortify(mask)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos)) +
    geom_raster(aes(fill=krige), interpolate=F) +
    scale_fill_gradientn(colors=p_palette, limits=c(minvalue, maxvalue), oob=scales::squish) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name, title=title_name) +
    geom_sf(data=mask, color='white', fill="white", linewidth=2, inherit.aes=F) +
    theme_void() +
    theme(legend.position="right")

  return(p)
}

