#' .fov_image_theme
#' @description
#' internal function to setup the theme for the fov background on the whole
#' image
#'
#' @param back.color not used
#' @param back.border color for the borders of the background (default=NA)
#' @param title.col character indicating the color of the title
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect element_text
#' @keywords internal
.fov_image_theme <- function(back.color="black", back.border=NA,
                             title.col="white")
{
    theme(panel.border=element_blank(),
          legend.key=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          panel.grid=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill = "transparent", colour = NA),
          plot.title=element_text(color=title.col, hjust=0.5, face="bold"),
          plot.background=element_rect(fill="transparent", colour=back.border))
}

#' .centroid_image_theme
#' @description
#' internal function to setup the theme for the centroid plot background
#'
#' @param back.color not used
#' @param back.border color for the borders of the background (default=NA)
#' @param title.col character indicating the color of the title
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect
#' @keywords internal
.centroid_image_theme <- function(back.border=NA){
    theme(panel.border=element_rect(color = "black"),
          legend.key=element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          panel.grid=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill = "transparent", colour = NA),
          plot.title=element_blank(),
          plot.background=element_rect(fill="transparent", colour=back.border))
}


#' .negative_image_theme
#' @description
#' internal function to setup the theme for the negative background for
#' negative plots
#'
#' @param fill_color color to fill the element_rect (default is "black")
#' @param fore_color color for all the other elements (default is "white")
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_line element_rect element_text element_blank
#' @keywords internal
.negative_image_theme <- function(fill_color="black", fore_color="white"){
    theme(panel.border = element_rect(color = fore_color),
          panel.background=element_rect(fill=fill_color, color=NA),
          plot.background=element_rect(fill=fill_color, color=NA),
          title=element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          legend.background=element_rect(fill=fill_color, color=NA),
          legend.text=element_text(color=fore_color),
          legend.title=element_text(color=fore_color),
          panel.grid=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank())
}

#' .light_theme
#' @description
#' internal function to setup the white background theme for the First Filter plot
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect
#' @keywords internal
.light_theme <- function(){
    theme(panel.border=element_rect(color = "black", linewidth = 0.1),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          panel.grid=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank())
}

#' .dark_theme
#' @description
#' internal function to setup the black background theme for the First Filter plot
#' #' @param fill_color color to fill the element_rect (default is "black")
#' @param fore_color color for all the other elements (default is "white")
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect element_text
#' @keywords internal
.dark_theme <- function(fill_color="black", fore_color="white"){
    theme(panel.border=element_blank(),
          panel.background=element_rect(fill=fill_color, color=NA),
          plot.background=element_rect(fill=fill_color, color=NA),
          title=element_text(color=fore_color),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          legend.background=element_rect(fill=fill_color, color=NA),
          legend.text=element_text(color=fore_color),
          panel.grid=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank())
}

