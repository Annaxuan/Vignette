library(imager) 
library(tidyverse) 
library(tidymodels) 
library(sp) 
library(scales)
library(cowplot) 
#devtools::install_github("sharlagelfand/dmc") 
library(dmc)

change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size.
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
}

process_image <- function(image_file_name, k_list){
  ## process_image(image_file_name, k_list) take an image and process this image to
  ## compute clusterings with different k values. Where k is in k_list. Gives 
  ## sufficient information for other functions. Referenced some codes in week2 tut.
  ##
  ## Input: 
  ## - image_file_name : a PNG or JPEG image.
  ## - k_list : list of the number of centres in the clustering
  ##
  ## Output: 
  ## - cluster_info: A tibble of information derived from the k_means method. Contains 
  ##   the original output of the kclust calls for each k, the tidied clusters,
  ##   their associated RGB values and their nearest DMC thread color information.
  ##
  ## Example:
  ##   library(imager) 
  ##   library(dmc)
  ##   library(tidyverse) 
  ##   library(tidymodels)
  ##   library(scales)
  ##   image_file_name = "image10.png"
  ##   k_list = 2:8
  ##   cluster_info <- process_image(image_file_name, k_list)
  
  im <- imager::load.image(image_file_name)
  tidy_dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  dat <- select(tidy_dat,c(-x,-y))
  kclusts <-
    tibble(k = c(k_list)) %>%
    mutate(
      kclust = map(k, ~kmeans(x = dat , centers = .x, nstart = 4)),
      glanced = map(kclust, glance)
    )
  
  # clusterings contains the tidied clusters.
  clusterings <- 
    kclusts %>%
    unnest(cols = c(glanced)) %>%
    mutate(centers = map(kclust, tidy), 
           tidy_dat = map(kclust, ~augment(.x, tidy_dat) %>%
                            rename(cluster = .cluster)))
  
  for (k in k_list) {
    center <- clusterings[clusterings$k==k,]$centers[[1]]
    center <- center %>%
      mutate(col = rgb(R, G, B)) %>%
      mutate(dmc = map(col, ~dmc(.x)))
    clusterings[clusterings$k==k,]$centers[[1]] <- center
  }
  return(clusterings)
}



scree_plot <- function(cluster_info){
  ## scree_plot(cluster_info) produces a ratio version scree plot. Referenced some 
  ## codes in extra_help-2.Rmd, which is provided by professor.
  ##
  ## Input: 
  ## - cluster_info: A tibble of information derived from the k_means method. Contains 
  ##   the original output of the kclust calls for each k, the tidied clusters,
  ##   their associated RGB values and their nearest DMC thread color information.
  ##
  ## Output: 
  ## - scree_graph: A ratio version scree plot that helps finding the right number of 
  ##   clusters (k value).
  ##
  ## Example:
  ##   library(imager) 
  ##   library(dmc)
  ##   library(tidyverse) 
  ##   library(tidymodels)
  ##   library(scales)
  ##   library(cowplot) 
  ##   library(sp) 
  ##   image_file_name = "image10.png"
  ##   k_list = 2:8
  ##   cluster_info <- process_image(image_file_name, k_list)
  ##   scree_plot(cluster_info)
  
  # ratio version scree plot.
  nclust = length(cluster_info$k)
  ratio = rep(NA, nclust-1)
  
  for (kk in 2:nclust){
    ratio[kk-1] = cluster_info$tot.withinss[kk]/cluster_info$tot.withinss[kk-1]
  }
  
  plot_data <- data.frame(k = cluster_info$k[2:nclust],ratio)
  ggplot(plot_data,aes(x=k, y = ratio)) + geom_line() + 
    ggtitle("Ratio version scree plot") + 
    labs(y="Ratio of total within sum squares", x = "number of clusters(k)") + 
    theme(plot.title = element_text(hjust = 0.5))
}


color_strips <- function(cluster_info){
  ## color_strips(cluster_info) produces color strips with the DMC color closest 
  ## to the cluster center color. Referenced some codes in extra_help-2.Rmd, 
  ## which is provided by professor.
  ##
  ## Input: 
  ## - cluster_info: A tibble of information derived from the k_means method. Contains 
  ##   the original output of the kclust calls for each k, the tidied clusters,
  ##   their associated RGB values and their nearest DMC thread color information.
  ##
  ## Output: 
  ## - A color strips has k rows. Each row shows plot for different k values
  ##   that contains clusterings and its center's closet DMC color value.
  ##
  ## Example:
  ##   library(imager) 
  ##   library(dmc)
  ##   library(tidyverse) 
  ##   library(tidymodels)
  ##   library(scales)
  ##   library(cowplot) 
  ##   library(sp) 
  ##   image_file_name = "image10.png"
  ##   k_list = 2:8
  ##   cluster_info <- process_image(image_file_name, k_list)
  ##   color_strips(cluster_info)
  
  color_lst <- list()
  for (index in 1:nrow(cluster_info)){
    colours2 = sapply(cluster_info$centers[[index]]$dmc, "[[", 3)
    n_col = length(colours2)
  
    rect_dat <- tibble(x1 = c(0:(n_col-1)), 
                     x2 = c(1:n_col), 
                     y1 = rep(0,n_col),
                     y2 =rep(1,n_col), 
                     colour = colours2)
  
    color <- rect_dat %>% ggplot() + coord_fixed() + 
      geom_rect(aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=colour), color="black") +
      geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=colour), size=3) + 
      scale_fill_manual(values = rect_dat$colour, breaks = rect_dat$colour) + 
      theme_void() + theme(legend.position = "none") 
  
    color_lst[[index]] = color
  }
  
  plot_grid(
    plotlist = color_lst,
    ncol = 1, 
    labels = c(cluster_info$k)
    )
}


make_pattern <- function(cluster_info, k, x_size, black_white, background_colour){
  ## make_pattern(cluster_info, k, x_size, black_white, background_colour) plots 
  ## the vignette pattern. Referenced some codes in extra_help-2.Rmd, which is
  ## provided by professor.
  ##
  ## Input: 
  ## - cluster_info: A tibble of information derived from the k_means method. Contains 
  ##   the original output of the kclust calls for each k, the tidied clusters,
  ##   their associated RGB values and their nearest DMC thread color information.
  ## - k: the number of clusters.
  ## - x_size: the number pixels of height/width. This is a square pattern, where 
  ##   height equals to width.
  ## - black_white: A boolean value, where TRUE means showing the pattern in black
  ##   or white, FALSE means showing the pattern in different color.
  ## - background_colour: A NULL value or a DMC code like "#XXXXXX" in string. 
  ##   Where NULL means that the background color is shown, where a DMC code "#XXXXXX"
  ##   means that the background "#XXXXXX" can not be shown in the pattern.
  ##
  ## Output: 
  ## - A vignette pattern of the image, with k number clusters, width and height of x_size, 
  ##   shown in black and white or colored.
  ##
  ## Precondition:
  ## - background_colour must be the correct DMC code shows in the color strip. 
  ##   The format of this string is "#XXXXXX".
  ## 
  ## Example:
  ##   library(imager) 
  ##   library(dmc)
  ##   library(tidyverse) 
  ##   library(tidymodels)
  ##   library(scales)
  ##   library(cowplot) 
  ##   library(sp) 
  ##   image_file_name = "image10.png"
  ##   k_list = 2:8
  ##   cluster_info <- process_image(image_file_name, k_list)
  ##   make_pattern(cluster_info, 4, 80, FALSE, NULL)
  
  k_cluster_info <- cluster_info[cluster_info$k==k,]
  
  image_df <- k_cluster_info$tidy_dat[[1]]
  
  true_df <- change_resolution(image_df, x_size)
  
  clus_col <- k_cluster_info$centers[[1]]
  
  clus_dmc <- clus_col$dmc
  
  clus_col$dmc_num <- sapply(clus_dmc, "[[", 1)
  
  true_df$dmc <- NA
  for (index in 1:nrow(true_df)){
    true_df[index,]$dmc <- clus_col[clus_col$cluster == true_df[index,]$cluster,]$dmc_num
  }
  
  dmc_t <- tibble(dmc = sapply(clus_dmc, "[[", 1), 
                  name = sapply(clus_dmc, "[[", 2), 
                  hex = sapply(clus_dmc, "[[", 3))
  
  # shape value list
  shape_l <- c()
  i = 0
  dmc_background = ""
  # extract background colour
  if (!is.null(background_colour)){
    background <- dmc_t %>% filter(dmc_t$hex==background_colour)
    dmc_background <- background$dmc
    dmc_t <- dmc_t %>% filter(dmc_t$hex!=background_colour)
    while (i < k-1){
      shape_l <- c(shape_l, i)
      i = i+1
    }
  }else{
    while (i < k){
      shape_l <- c(shape_l, i)
      i = i+1
    }
  }
  
  # plot the pattern
  if (black_white){
    graph <- true_df %>% filter(true_df$dmc!=dmc_background) %>% 
      ggplot(aes(x, y)) + geom_point(aes(shape = factor(dmc))) +
      scale_colour_manual(values = dmc_t %>% select(dmc, hex) %>% deframe,
                          label =  dmc_t %>% select(dmc, name) %>% deframe) + 
      scale_shape_manual(values = shape_l, labels =  dmc_t %>% select(dmc, name) %>% deframe) +
      scale_y_reverse() + theme_void() +
      background_grid(major = "xy",minor = "xy", size.major = 0.8, size.minor = 0.3, 
                      color.major = "black") + 
      ggtitle("Vignette pattern in black and white") + theme(plot.title = element_text(hjust = 0.5))
  }else{
    graph <- true_df %>% filter(true_df$dmc!=dmc_background) %>% ggplot(aes(x, y)) + 
      geom_point(aes(shape = factor(dmc), col = factor(dmc))) +
      scale_colour_manual(values = dmc_t %>% select(dmc, hex) %>% deframe,
                          labels =  dmc_t %>% select(dmc, name) %>% deframe) +
      scale_shape_manual(values = shape_l, labels =  dmc_t %>% select(dmc, name) %>% deframe) +
      scale_y_reverse() + theme_void() +
      background_grid(major = "xy",minor = "xy", size.major = 0.8, size.minor = 0.3, 
                      color.major = "black", color.minor = "grey") + 
      ggtitle("Vignette pattern with four colors") + theme(plot.title = element_text(hjust = 0.5))
  }
  graph
}
