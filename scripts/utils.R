# All code in this R script has written by Evan Collins. Any questions can be sent to evantomcollins@gmail.com.

# Load packages
library(dplyr)
library(plotly)
library(ggplot2)
library(Rmisc)
library(png)
library(grid)
library(ggpubr)

# Load files
# Yale Brain Atlas parcel dictionary
parcel_dict <- read.csv("../data/misc/parcel_dict.csv")
parcel_names <- parcel_dict$Name
vertex_number_by_parcel <- parcel_dict$Num_vertex
# Yale Brain Atlas whole brain positions and indices
atlas_whole_positions <- read.csv("../data/atlas/atlas_whole_positions.csv")
atlas_whole_indices <- read.csv("../data/atlas/atlas_whole_indices.csv")
# Yale Brain Atlas left hemisphere positions and indices
atlas_LH_positions <- read.csv("../data/atlas/atlas_LH_positions.csv")
atlas_LH_indices <- read.csv("../data/atlas/atlas_LH_indices.csv")
# Yale Brain Atlas right hemisphere positions and indices
atlas_RH_positions <- read.csv("../data/atlas/atlas_RH_positions.csv")
atlas_RH_indices <- read.csv("../data/atlas/atlas_RH_indices.csv")
# Yale Brain Atlas deep left hemisphere (amygdala, hippocampus, insula) positions and indices
atlas_deepLH_positions <- read.csv("../data/atlas/atlas_LH_positions.csv")
atlas_deepLH_indices <- read.csv("../data/atlas/atlas_LH_indices.csv")
# Yale Brain Atlas deep right hemisphere (amygdala, hippocampus, insula) positions and indices
atlas_deepRH_positions <- read.csv("../data/atlas/atlas_deepRH_positions.csv")
atlas_deepRH_indices <- read.csv("../data/atlas/atlas_deepRH_indices.csv")



# Function to find the closest point among known points for input coordinates
find_closest_points <- function(input_x, input_y, input_z, known_points) {
  closest_indices <- numeric(length(input_x))
  for (i in seq_along(input_x)) {
    min_distance <- Inf
    closest_point_index <- -1
    for (j in seq_along(known_points$x)) {
      distance <- sqrt((known_points$x[j] - input_x[i])^2 + (known_points$y[j] - input_y[i])^2 + (known_points$z[j] - input_z[i])^2)
      if (distance < min_distance) {
        min_distance <- distance
        closest_point_index <- j
      }
    }
    closest_indices[i] <- closest_point_index
  }
  return(closest_indices)
}


# Function to assign data in MNI152 space to Yale Brain Atlas parcels
get_parcel_labels_for_coords <- function(x_coords, y_coords, z_coords){

  stopifnot("x_coords must be a vector containing numeric or integer values." = is.vector(x_coords) & (all(is.numeric(x_coords) | is.integer(x_coords))))
  stopifnot("y_coords must be a vector containing numeric or integer values." = is.vector(y_coords) & (all(is.numeric(y_coords) | is.integer(y_coords))))
  stopifnot("z_coords must be a vector containing numeric or integer values." = is.vector(z_coords) & (all(is.numeric(z_coords) | is.integer(z_coords))))

  # Load atlas coordinates with parcel labels
  atlas_whole_positions <- read.csv("../data/atlas/atlas_whole_positions.csv")

  # Compute closest index in atlas for each inputted point
  closest_indices <- find_closest_points(x_coords, y_coords, z_coords, atlas_whole_positions)
  # Translate index into parcel label
  parcel_labels <- atlas_whole_positions$parcel[closest_indices]
  # Output dataframe containing inputted coordinates and parcel labels
  labeled_coord_df <- data.frame("x" = x_coords,
                                 "y" = y_coords,
                                 "z" = z_coords,
                                 "parcel" = parcel_labels)
  return(labeled_coord_df)
}


# Function to compute average value by parcel
calculate_average_value_by_parcel <- function(values, parcel_labels){

  stopifnot("values must be a vector containing numeric or integer values." = is.vector(values) & (all(is.numeric(values) | is.integer(values))))
  stopifnot("parcel_labels must be a vector containing character values of parcel names." = is.vector(parcel_labels) & (all(is.character(parcel_labels) & all(parcel_labels %in% parcel_names))))
  stopifnot("values must be same length as parcel_labels." = length(values) == length(parcel_labels))

  data_df <- data.frame("value" = values,
                        "parcel_label" = parcel_labels)
  average_value_by_parcel_df <- data.frame("parcel" = parcel_names,
                                           "N" = 0,
                                           "average_value" = 0,
                                           "stdev" = 0)
  data_stats_df <- summarySE(data_df, measurevar = "value", groupvars = c("parcel_label"))
  data_stats_df$sd[is.na(data_stats_df$sd)] <- 0

  for (i in 1:nrow(data_stats_df)){
    average_value_by_parcel_df$N[average_value_by_parcel_df$parcel == data_stats_df$parcel_label[i]] <- data_stats_df$N[i]
    average_value_by_parcel_df$average_value[average_value_by_parcel_df$parcel == data_stats_df$parcel_label[i]] <- data_stats_df$value[i]
    average_value_by_parcel_df$stdev[which(average_value_by_parcel_df$parcel == data_stats_df$parcel_label[i])] <- data_stats_df$sd[i]
  }
  return(average_value_by_parcel_df)
}


# Function to get RGB values for Yale Brain Atlas parcels and ggplot legend based on inputted values along defined lower and upper bound colors
get_rgb_lu_for_parcel_values <- function(parcel_values, lower_bound_color = c(255, 255, 255), upper_bound_color = c(255, 0, 0), legend = FALSE, legend_title = NULL, legend_values = FALSE){

  stopifnot("parcel_values must be a vector containing numeric or integer values." = is.vector(parcel_values) & (all(is.numeric(parcel_values) | is.integer(parcel_values))))
  stopifnot("lower_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(lower_bound_color) & (all(is.numeric(lower_bound_color) | is.integer(lower_bound_color))) & (length(lower_bound_color) == 3) & all(lower_bound_color >= 0 & lower_bound_color <= 255))
  stopifnot("upper_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(upper_bound_color) & (all(is.numeric(upper_bound_color) | is.integer(upper_bound_color))) & (length(upper_bound_color) == 3) & all(upper_bound_color >= 0 & upper_bound_color <= 255))
  stopifnot("legend must be either TRUE (1) or FALSE (0)." = legend == TRUE | legend == FALSE)
  stopifnot("legend_title must be a single character element." = is.null(legend_title) | class(legend_title) == "call" | (is.character(legend_title) & (length(legend_title) == 1)))
  stopifnot("legend_values must be either TRUE (1) or FALSE (0)." = legend_values == TRUE | legend_values == FALSE)

  parcel_rgb_df <- data.frame("values" = parcel_values)
  parcel_rgb_legend_fig <- ggplot(parcel_rgb_df,
                                  aes(x=values,
                                      y=values,
                                      color=values)) +
    geom_point() +
    scale_color_gradientn(limits = c(min(parcel_rgb_df$values), max(parcel_rgb_df$values)),
                          colors = c(rgb(lower_bound_color[1]/255, lower_bound_color[2]/255, lower_bound_color[3]/255), rgb(upper_bound_color[1]/255, upper_bound_color[2]/255, upper_bound_color[3]/255)),
                          name = legend_title) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.title=element_text(size=18),
          legend.text=element_text(size=15),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(2, 'cm'))

  if (legend_values == FALSE){
    parcel_rgb_legend_fig <- parcel_rgb_legend_fig + theme(legend.text = element_blank())
  }

  parcel_rgb_df$R <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[1,]
  parcel_rgb_df$G <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[2,]
  parcel_rgb_df$B <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[3,]

  parcel_rgb_df_legend <- parcel_rgb_df

  if (legend == TRUE){
    parcel_rgb_legend_data <- get_legend(parcel_rgb_legend_fig)
    parcel_rgb_legend <- as_ggplot(parcel_rgb_legend_data)
    parcel_rgb_df_legend <- list(parcel_rgb_df, parcel_rgb_legend)
  }

  return(parcel_rgb_df_legend)
}


# Function to get RGB values for Yale Brain Atlas parcels and ggplot legend based on inputted values along defined lower, middle, upper bound colors
get_rgb_lmu_for_parcel_values <- function(parcel_values, lower_bound_color = c(0, 0, 255), middle_bound_color = c(255, 255, 255), upper_bound_color = c(255, 0, 0), zeroed = FALSE, legend = FALSE, legend_title = NULL, legend_values = FALSE){

  stopifnot("parcel_values must be a vector containing numeric or integer values." = is.vector(parcel_values) & (all(is.numeric(parcel_values) | is.integer(parcel_values))))
  stopifnot("lower_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(lower_bound_color) & (all(is.numeric(lower_bound_color) | is.integer(lower_bound_color))) & (length(lower_bound_color) == 3) & all(lower_bound_color >= 0 & lower_bound_color <= 255))
  stopifnot("middle_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(middle_bound_color) & (all(is.numeric(middle_bound_color) | is.integer(middle_bound_color))) & (length(middle_bound_color) == 3) & all(middle_bound_color >= 0 & middle_bound_color <= 255))
  stopifnot("upper_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(upper_bound_color) & (all(is.numeric(upper_bound_color) | is.integer(upper_bound_color))) & (length(upper_bound_color) == 3) & all(upper_bound_color >= 0 & upper_bound_color <= 255))
  stopifnot("zeroed must be either TRUE (1) or FALSE (0)." = zeroed == TRUE | zeroed == FALSE)
  stopifnot("legend must be either TRUE (1) or FALSE (0)." = legend == TRUE | legend == FALSE)
  stopifnot("legend_title must be a single character element." = is.null(legend_title) | (is.character(legend_title) & (length(legend_title) == 1)))
  stopifnot("legend_values must be either TRUE (1) or FALSE (0)." = legend_values == TRUE | legend_values == FALSE)

  parcel_rgb_df <- data.frame("values" = parcel_values)
  if (zeroed == FALSE){
    parcel_rgb_legend_fig <- ggplot(parcel_rgb_df,
                                    aes(x=values,
                                        y=values,
                                        color=values)) +
      geom_point() +
      scale_color_gradientn(limits = c(min(parcel_rgb_df$values), max(parcel_rgb_df$values)),
                            colors = c(rgb(lower_bound_color[1]/255, lower_bound_color[2]/255, lower_bound_color[3]/255), rgb(middle_bound_color[1]/255, middle_bound_color[2]/255, middle_bound_color[3]/255), rgb(upper_bound_color[1]/255, upper_bound_color[2]/255, upper_bound_color[3]/255)),
                            name = legend_title) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=15),
            legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(2, 'cm'))
  }

  if (zeroed == TRUE){
    parcel_rgb_legend_fig <- ggplot(parcel_rgb_df,
                                    aes(x=values,
                                        y=values,
                                        color=values)) +
      geom_point() +
      scale_color_gradientn(limits = c(-pmax(abs(min(parcel_rgb_df$values)),
                                             abs(max(parcel_rgb_df$values))),
                                       pmax(abs(min(parcel_rgb_df$values)),
                                            abs(max(parcel_rgb_df$values)))),
                            colors = c(rgb(lower_bound_color[1]/255, lower_bound_color[2]/255, lower_bound_color[3]/255), rgb(middle_bound_color[1]/255, middle_bound_color[2]/255, middle_bound_color[3]/255), rgb(upper_bound_color[1]/255, upper_bound_color[2]/255, upper_bound_color[3]/255)),
                            name = legend_title) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=15),
            legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(2, 'cm'))
  }

  if (legend_values == FALSE){
    parcel_rgb_legend_fig <- parcel_rgb_legend_fig + theme(legend.text = element_blank())
  }

  parcel_rgb_df$R <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[1,]
  parcel_rgb_df$G <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[2,]
  parcel_rgb_df$B <- col2rgb(ggplot_build(parcel_rgb_legend_fig)$data[[1]][[1]])[3,]

  parcel_rgb_df_legend <- parcel_rgb_df

  if (legend == TRUE){
    parcel_rgb_legend_data <- get_legend(parcel_rgb_legend_fig)
    parcel_rgb_legend <- as_ggplot(parcel_rgb_legend_data)
    parcel_rgb_df_legend <- list(parcel_rgb_df, parcel_rgb_legend)
  }

  return(parcel_rgb_df_legend)
}


# Function to generate 3D plot of Yale Brain Atlas
plot_brain_3d <- function(R_values = parcel_dict$YBA_R_color, G_values = parcel_dict$YBA_G_color, B_values = parcel_dict$YBA_B_color, brain_section = "whole", extra_text_varname = "Region", extra_text_values = parcel_dict$Region, hide_axes = TRUE) {

  stopifnot("R_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(R_values) | (is.vector(R_values) & (all(is.numeric(R_values) | is.integer(R_values))) & (length(R_values) == 690) & all(R_values >= 0 & R_values <= 255)))
  stopifnot("G_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(G_values) | (is.vector(G_values) & (all(is.numeric(G_values) | is.integer(G_values))) & (length(G_values) == 690) & all(G_values >= 0 & G_values <= 255)))
  stopifnot("B_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(B_values) | (is.vector(B_values) & (all(is.numeric(B_values) | is.integer(B_values))) & (length(B_values) == 690) & all(B_values >= 0 & B_values <= 255)))
  stopifnot("brain_section must be either 'whole', 'LH', 'RH', 'deepLH', or 'deepRH' as a single character element." = (is.character(brain_section) & (brain_section %in% c('whole', 'LH', 'RH', 'deepLH', 'deepRH'))))
  stopifnot("extra_text_varname must be a single character element." = is.null(extra_text_varname) | (is.character(extra_text_varname) & (length(extra_text_varname) == 1)))
  stopifnot("extra_text_values must be a vector of length 690." = is.null(extra_text_values) | (is.vector(extra_text_values) & (length(extra_text_values) == 690)))
  stopifnot("hide_axes must be either TRUE (1) or FALSE (0)." = hide_axes == TRUE | hide_axes == FALSE)

  if (brain_section == "whole"){
    atlas_section_positions <- atlas_whole_positions
    atlas_section_indices <- atlas_whole_indices
  } else if (brain_section == "LH"){
    atlas_section_positions <- atlas_LH_positions
    atlas_section_indices <- atlas_LH_indices
  } else if (brain_section == "RH"){
    atlas_section_positions <- atlas_RH_positions
    atlas_section_indices <- atlas_RH_indices
  } else if (brain_section == "deepLH"){
    atlas_section_positions <- atlas_deepLH_positions
    atlas_section_indices <- atlas_deepLH_indices
  } else if (brain_section == "deepRH"){
    atlas_section_positions <- atlas_deepRH_positions
    atlas_section_indices <- atlas_deepRH_indices
  }

  brain_3d <- plot_ly(
    type = "mesh3d",
    x = atlas_section_positions$x,
    y = atlas_section_positions$y,
    z = atlas_section_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict[,brain_section] == 1)]/255, G_values[which(parcel_dict[,brain_section] == 1)]/255, B_values[which(parcel_dict[,brain_section] == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)]),
    opacity = 1,
    i = atlas_section_indices$i,
    j = atlas_section_indices$j,
    k = atlas_section_indices$k,
    hovertext = paste0(
      "Parcel: ",
      rep(parcel_dict$Name[which(parcel_dict[,brain_section] == 1)], vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)]),
      "<br>", extra_text_varname,": ",
      rep(extra_text_values[which(parcel_dict[,brain_section] == 1)], vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)])
    )
  ) %>%
    layout(hoverlabel = list(
      bgcolor = "white",
      bordercolor = "darkblue",
      font = list(size = 15, family = "Verdana")
    )) %>%
    config(displaylogo = F, displayModeBar = T) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    )))

  if (hide_axes == TRUE){
    brain_3d <- brain_3d %>%
      layout(scene = list(xaxis = list(visible = FALSE),
                          yaxis = list(visible = FALSE),
                          zaxis = list(visible = FALSE)))
  }

  return(brain_3d)
}


# Function to generate a single image of six views of Yale Brain Atlas
plot_brain_map <- function(R_values = parcel_dict$YBA_R_color, G_values = parcel_dict$YBA_G_color, B_values = parcel_dict$YBA_B_color, plot_title = NULL, legend_figure = NULL){

  stopifnot("R_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(R_values) | (is.vector(R_values) & (all(is.numeric(R_values) | is.integer(R_values))) & (length(R_values) == 690) & all(R_values >= 0 & R_values <= 255)))
  stopifnot("G_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(G_values) | (is.vector(G_values) & (all(is.numeric(G_values) | is.integer(G_values))) & (length(G_values) == 690) & all(G_values >= 0 & G_values <= 255)))
  stopifnot("B_values must be a vector of length 690 containing numeric or integer values on the RGB scale from 0 to 255." = is.null(B_values) | (is.vector(B_values) & (all(is.numeric(B_values) | is.integer(B_values))) & (length(B_values) == 690) & all(B_values >= 0 & B_values <= 255)))
  stopifnot("plot_title must be a single character element." = is.null(plot_title) | (is.character(plot_title) & (length(plot_title) == 1)))
  stopifnot("legend_figure must be a ggplot object." = is.null(legend_figure) | is.ggplot(legend_figure))

  brain_list <- list()
  brain_list[[1]] <- plot_ly(
    type = "mesh3d",
    x = atlas_whole_positions$x,
    y = atlas_whole_positions$y,
    z = atlas_whole_positions$z,
    vertexcolor = rep(
      rgb(R_values/255, G_values/255, B_values/255, 1),
      vertex_number_by_parcel
    ),
    i = atlas_whole_indices$i,
    j = atlas_whole_indices$j,
    k = atlas_whole_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 0, y = 0, z = 1.4),
      up = list(x = 0, y = 1, z = 0),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[2]] <- plot_ly(
    type = "mesh3d",
    x = atlas_whole_positions$x,
    y = atlas_whole_positions$y,
    z = atlas_whole_positions$z,
    vertexcolor = rep(
      rgb(R_values/255, G_values/255, B_values/255, 1),
      vertex_number_by_parcel
    ),
    i = atlas_whole_indices$i,
    j = atlas_whole_indices$j,
    k = atlas_whole_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 0, y = 0, z = -1.5),
      up = list(x = 0, y = 1, z = 0),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[3]] <- plot_ly(
    type = "mesh3d",
    x = atlas_LH_positions$x,
    y = atlas_LH_positions$y,
    z = atlas_LH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$LH == 1)]/255, G_values[which(parcel_dict$LH == 1)]/255, B_values[which(parcel_dict$LH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$LH == 1)]),
    opacity = 1,
    i = atlas_LH_indices$i,
    j = atlas_LH_indices$j,
    k = atlas_LH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = -1.3, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[4]] <- plot_ly(
    type = "mesh3d",
    x = atlas_LH_positions$x,
    y = atlas_LH_positions$y,
    z = atlas_LH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$LH == 1)]/255, G_values[which(parcel_dict$LH == 1)]/255, B_values[which(parcel_dict$LH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$LH == 1)]),
    opacity = 1,
    i = atlas_LH_indices$i,
    j = atlas_LH_indices$j,
    k = atlas_LH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[5]] <- plot_ly(
    type = "mesh3d",
    x = atlas_RH_positions$x,
    y = atlas_RH_positions$y,
    z = atlas_RH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$RH == 1)]/255, G_values[which(parcel_dict$RH == 1)]/255, B_values[which(parcel_dict$RH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$RH == 1)]),
    i = atlas_RH_indices$i,
    j = atlas_RH_indices$j,
    k = atlas_RH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.3, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[6]] <- plot_ly(
    type = "mesh3d",
    x = atlas_RH_positions$x,
    y = atlas_RH_positions$y,
    z = atlas_RH_positions$z,
    vertexcolor = rep(
      rgb((R_values[which(parcel_dict$RH == 1)])/255, G_values[which(parcel_dict$RH == 1)]/255, B_values[which(parcel_dict$RH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$RH == 1)]),
    i = atlas_RH_indices$i,
    j = atlas_RH_indices$j,
    k = atlas_RH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = -1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  suppressWarnings(orca(brain_list[[1]], "../data/brain_map/brain1.png", width = 1500, height = 1500))
  suppressWarnings(orca(brain_list[[2]], "../data/brain_map/brain2.png", width = 1500, height = 1500))
  suppressWarnings(orca(brain_list[[3]], "../data/brain_map/brain3.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[4]], "../data/brain_map/brain4.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[5]], "../data/brain_map/brain5.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[6]], "../data/brain_map/brain6.png", width = 1800, height = 1300))


  brain_img_1 <- readPNG("../data/brain_map/brain1.png")
  brain_img_2 <- readPNG("../data/brain_map/brain2.png")
  brain_img_3 <- readPNG("../data/brain_map/brain3.png")
  brain_img_4 <- readPNG("../data/brain_map/brain4.png")
  brain_img_5 <- readPNG("../data/brain_map/brain5.png")
  brain_img_6 <- readPNG("../data/brain_map/brain6.png")

  brain_figs <- ggplot() +
    theme_classic() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.line=element_blank()) +
    annotation_custom(rasterGrob(brain_img_1, width = 1, height = 1),
                      xmin = 0, xmax = 0.294,
                      ymin = 0.5, ymax = 1) +
    annotation_custom(rasterGrob(brain_img_2, width = 1, height = 1),
                      xmin = 0, xmax = 0.294,
                      ymin = 0, ymax = 0.5) +
    annotation_custom(rasterGrob(brain_img_3, width = 1, height = 1),
                      xmin = 0.294, xmax = 0.647,
                      ymin = 0.5, ymax = 0.9333) +
    annotation_custom(rasterGrob(brain_img_4, width = 1, height = 1),
                      xmin = 0.294, xmax = 0.647,
                      ymin = 0.0667, ymax = 0.5) +
    annotation_custom(rasterGrob(brain_img_5, width = 1, height = 1),
                      xmin = 0.647, xmax = 1,
                      ymin = 0.5, ymax = 0.9333) +
    annotation_custom(rasterGrob(brain_img_6, width = 1, height = 1),
                      xmin = 0.647, xmax = 1,
                      ymin = 0.0667, ymax = 0.5) +
    coord_fixed(ratio = 0.5882)

  if (!is.null(legend_figure)){
    brain_map <- ggarrange(brain_figs, legend_figure, ncol = 2, nrow = 1, widths = c(0.85, 0.15))
  } else {
    brain_map <- brain_figs
  }

  return(brain_map)
}


# Function to get RGB values for Yale Brain Atlas parcels and ggplot legend based on inputted values along defined lower, middle, upper bound colors for a specific Neurosynth term
get_rgb_for_neurosynth_term <- function(functional_term, lower_bound_color = c(0, 0, 255), middle_bound_color = c(255, 255, 255), upper_bound_color = c(255, 0, 0), zeroed = TRUE, legend = FALSE, legend_title = NULL){

  stopifnot("functional_term must be a single character element containing a functional term found in Neurosynth. For a list of all 1334 functional terms, visit https://neurosynth.org/analyses/terms/" = (is.character(functional_term) & (length(functional_term) == 1) & (functional_term %in% colnames(parcelsynth_df))))
  stopifnot("parcel_values must be a vector containing numeric or integer values." = is.vector(values) & (all(is.numeric(values) | is.integer(values))))
  stopifnot("lower_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(lower_bound_color) & (all(is.numeric(lower_bound_color) | is.integer(lower_bound_color))) & (length(lower_bound_color) == 3) & all(lower_bound_color >= 0 & lower_bound_color <= 255))
  stopifnot("middle_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(middle_bound_color) & (all(is.numeric(middle_bound_color) | is.integer(middle_bound_color))) & (length(middle_bound_color) == 3) & all(middle_bound_color >= 0 & middle_bound_color <= 255))
  stopifnot("upper_bound_color must be a vector of length 3 (R, G, B) containing numeric or integer values on the RGB scale from 0 to 255." = is.vector(upper_bound_color) & (all(is.numeric(upper_bound_color) | is.integer(upper_bound_color))) & (length(upper_bound_color) == 3) & all(upper_bound_color >= 0 & upper_bound_color <= 255))
  stopifnot("zeroed must be either TRUE (1) or FALSE (0)." = zeroed == TRUE | zeroed == FALSE)
  stopifnot("legend must be either TRUE (1) or FALSE (0)." = legend == TRUE | legend == FALSE)
  stopifnot("legend_title must be a single character element." = is.null(legend_title) | (is.character(legend_title) & (length(legend_title) == 1)))

  # Get parcel_values based on Neurosynth values for that term across 690 YBA parcels
  parcel_values <- parcelsynth_df[,functional_term]

  # Function to interpolate between two colors based on a ratio
  interpolate_colors <- function(color1, color2, ratio) {
    result_color <- (1 - ratio) * color1 + ratio * color2
    round(result_color)
  }

  # Normalize input values to range [0, 1]
  if (zeroed == TRUE){
    max_value <- max(c(abs(min(parcel_values)), abs(max(parcel_values))))
    min_value <- -max_value
  } else {
    min_value <- min(parcel_values)
    max_value <- max(parcel_values)
  }
  normalized_values <- (parcel_values - min_value) / (max_value - min_value)

  # Initialize output data frame
  parcel_rgb_df <- data.frame("R" = numeric(0), "G" = numeric(0), "B" = numeric(0))

  # Interpolate colors for each value in the input vector
  for (value in normalized_values) {
    if (value <= 0.5) {
      # Interpolate between lower_bound_color and middle_bound_color
      ratio <- value / 0.5
      color <- interpolate_colors(lower_bound_color, middle_bound_color, ratio)
    } else {
      # Interpolate between middle_bound_color and upper_bound_color
      ratio <- (value - 0.5) / 0.5
      color <- interpolate_colors(middle_bound_color, upper_bound_color, ratio)
    }

    parcel_rgb_df <- rbind(parcel_rgb_df, color)
  }
  names(parcel_rgb_df) <- c("R", "G", "B")
  parcel_rgb_df$values <- parcel_values

  if (legend == TRUE & zeroed == FALSE){
    parcel_rgb_legend_fig <- ggplot(parcel_rgb_df,
                                    aes(x=values,
                                        y=values,
                                        color=values)) +
      geom_point() +
      scale_color_gradientn(limits = c(min(parcel_rgb_df$values), max(parcel_rgb_df$values)),
                            colors = c(rgb(lower_bound_color[1]/255, lower_bound_color[2]/255, lower_bound_color[3]/255), rgb(middle_bound_color[1]/255, middle_bound_color[2]/255, middle_bound_color[3]/255), rgb(upper_bound_color[1]/255, upper_bound_color[2]/255, upper_bound_color[3]/255)),
                            name = legend_title) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=15),
            legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(2, 'cm'))
    parcel_rgb_legend_data <- get_legend(parcel_rgb_legend_fig)
    parcel_rgb_legend <- as_ggplot(parcel_rgb_legend_data)
  }

  if (legend == TRUE & zeroed == TRUE){
    parcel_rgb_legend_fig <- ggplot(parcel_rgb_df,
                                    aes(x=values,
                                        y=values,
                                        color=values)) +
      geom_point() +
      scale_color_gradientn(limits = c(-pmax(abs(min(parcel_rgb_df$values)),
                                             abs(max(parcel_rgb_df$values))),
                                       pmax(abs(min(parcel_rgb_df$values)),
                                            abs(max(parcel_rgb_df$values)))),
                            colors = c(rgb(lower_bound_color[1]/255, lower_bound_color[2]/255, lower_bound_color[3]/255), rgb(middle_bound_color[1]/255, middle_bound_color[2]/255, middle_bound_color[3]/255), rgb(upper_bound_color[1]/255, upper_bound_color[2]/255, upper_bound_color[3]/255)),
                            name = legend_title) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=15),
            legend.key.size = unit(1, 'cm'),
            legend.key.height = unit(2, 'cm'))
    parcel_rgb_legend_data <- get_legend(parcel_rgb_legend_fig)
    parcel_rgb_legend <- as_ggplot(parcel_rgb_legend_data)
  }

  parcel_rgb_df_legend <- list(parcel_rgb_df, parcel_rgb_legend)

  return(parcel_rgb_df_legend)
}


# Function to generate 3D plot of Yale Brain Atlas colored according to Neurosynth term
plot_neurosynth_brain_3d_for_term <- function(functional_term, brain_section = "whole", hide_axes = TRUE) {

  stopifnot("functional_term must be a single character element containing a functional term found in Neurosynth. For a list of all 1334 functional terms, visit https://neurosynth.org/analyses/terms/" = (is.character(functional_term) & (length(functional_term) == 1) & (functional_term %in% colnames(parcelsynth_df))))
  stopifnot("brain_section must be either 'whole', 'LH', 'RH', 'deepLH', or 'deepRH' as a single character element." = (is.character(brain_section) & (brain_section %in% c('whole', 'LH', 'RH', 'deepLH', 'deepRH'))))
  stopifnot("hide_axes must be either TRUE (1) or FALSE (0)." = hide_axes == TRUE | hide_axes == FALSE)

  parcel_values <- parcelsynth_df[,functional_term]
  functional_term_rgb_df <- get_rgb_for_neurosynth_term(functional_term = functional_term, legend = FALSE)[[1]]
  R_values <- functional_term_rgb_df$R
  G_values <- functional_term_rgb_df$G
  B_values <- functional_term_rgb_df$B

  if (brain_section == "whole"){
    atlas_section_positions <- atlas_whole_positions
    atlas_section_indices <- atlas_whole_indices
  } else if (brain_section == "LH"){
    atlas_section_positions <- atlas_LH_positions
    atlas_section_indices <- atlas_LH_indices
  } else if (brain_section == "RH"){
    atlas_section_positions <- atlas_RH_positions
    atlas_section_indices <- atlas_RH_indices
  } else if (brain_section == "deepLH"){
    atlas_section_positions <- atlas_deepLH_positions
    atlas_section_indices <- atlas_deepLH_indices
  } else if (brain_section == "deepRH"){
    atlas_section_positions <- atlas_deepRH_positions
    atlas_section_indices <- atlas_deepRH_indices
  }

  brain_3d <- plot_ly(
    type = "mesh3d",
    x = atlas_section_positions$x,
    y = atlas_section_positions$y,
    z = atlas_section_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict[,brain_section] == 1)]/255, G_values[which(parcel_dict[,brain_section] == 1)]/255, B_values[which(parcel_dict[,brain_section] == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)]),
    opacity = 1,
    i = atlas_section_indices$i,
    j = atlas_section_indices$j,
    k = atlas_section_indices$k,
    hovertext = paste0(
      "Parcel: ",
      rep(parcel_dict$Name[which(parcel_dict[,brain_section] == 1)], vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)]),
      "<br>Z-score: ",
      rep(round(parcel_values[which(parcel_dict[,brain_section] == 1)], 3), vertex_number_by_parcel[which(parcel_dict[,brain_section] == 1)])
    )
  ) %>%
    layout(hoverlabel = list(
      bgcolor = "white",
      bordercolor = "darkblue",
      font = list(size = 15, family = "Verdana")
    )) %>%
    config(displaylogo = F, displayModeBar = T) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    )))

  if (hide_axes == TRUE){
    brain_3d <- brain_3d %>%
      layout(scene = list(xaxis = list(visible = FALSE),
                          yaxis = list(visible = FALSE),
                          zaxis = list(visible = FALSE)))
  }

  return(brain_3d)
}


# Function to generate a single image of six views of Yale Brain Atlas for a specific Neurosynth term
plot_neurosynth_brain_map_for_term <- function(functional_term, plot_title = NULL){

  stopifnot("functional_term must be a single character element containing a functional term found in Neurosynth. For a list of all 1334 functional terms, visit https://neurosynth.org/analyses/terms/" = (is.character(functional_term) & (length(functional_term) == 1) & (functional_term %in% colnames(parcelsynth_df))))
  stopifnot("plot_title must be a single character element." = is.null(plot_title) | (is.character(plot_title) & (length(plot_title) == 1)))


  functional_term_rgb_df <- get_rgb_for_neurosynth_term(functional_term = functional_term, legend = TRUE, legend_title = functional_term)[[1]]
  R_values <- functional_term_rgb_df$R
  G_values <- functional_term_rgb_df$G
  B_values <- functional_term_rgb_df$B
  functional_term_rgb_legend <- get_rgb_for_neurosynth_term(functional_term = functional_term, legend = TRUE, legend_title = functional_term)[[2]]

  brain_list <- list()
  brain_list[[1]] <- plot_ly(
    type = "mesh3d",
    x = atlas_whole_positions$x,
    y = atlas_whole_positions$y,
    z = atlas_whole_positions$z,
    vertexcolor = rep(
      rgb(R_values/255, G_values/255, B_values/255, 1),
      vertex_number_by_parcel
    ),
    i = atlas_whole_indices$i,
    j = atlas_whole_indices$j,
    k = atlas_whole_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 0, y = 0, z = 1.4),
      up = list(x = 0, y = 1, z = 0),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[2]] <- plot_ly(
    type = "mesh3d",
    x = atlas_whole_positions$x,
    y = atlas_whole_positions$y,
    z = atlas_whole_positions$z,
    vertexcolor = rep(
      rgb(R_values/255, G_values/255, B_values/255, 1),
      vertex_number_by_parcel
    ),
    i = atlas_whole_indices$i,
    j = atlas_whole_indices$j,
    k = atlas_whole_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 0, y = 0, z = -1.5),
      up = list(x = 0, y = 1, z = 0),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[3]] <- plot_ly(
    type = "mesh3d",
    x = atlas_LH_positions$x,
    y = atlas_LH_positions$y,
    z = atlas_LH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$LH == 1)]/255, G_values[which(parcel_dict$LH == 1)]/255, B_values[which(parcel_dict$LH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$LH == 1)]),
    opacity = 1,
    i = atlas_LH_indices$i,
    j = atlas_LH_indices$j,
    k = atlas_LH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = -1.3, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[4]] <- plot_ly(
    type = "mesh3d",
    x = atlas_LH_positions$x,
    y = atlas_LH_positions$y,
    z = atlas_LH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$LH == 1)]/255, G_values[which(parcel_dict$LH == 1)]/255, B_values[which(parcel_dict$LH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$LH == 1)]),
    opacity = 1,
    i = atlas_LH_indices$i,
    j = atlas_LH_indices$j,
    k = atlas_LH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[5]] <- plot_ly(
    type = "mesh3d",
    x = atlas_RH_positions$x,
    y = atlas_RH_positions$y,
    z = atlas_RH_positions$z,
    vertexcolor = rep(
      rgb(R_values[which(parcel_dict$RH == 1)]/255, G_values[which(parcel_dict$RH == 1)]/255, B_values[which(parcel_dict$RH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$RH == 1)]),
    i = atlas_RH_indices$i,
    j = atlas_RH_indices$j,
    k = atlas_RH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = 1.3, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  brain_list[[6]] <- plot_ly(
    type = "mesh3d",
    x = atlas_RH_positions$x,
    y = atlas_RH_positions$y,
    z = atlas_RH_positions$z,
    vertexcolor = rep(
      rgb((R_values[which(parcel_dict$RH == 1)])/255, G_values[which(parcel_dict$RH == 1)]/255, B_values[which(parcel_dict$RH == 1)]/255, 1), vertex_number_by_parcel[which(parcel_dict$RH == 1)]),
    i = atlas_RH_indices$i,
    j = atlas_RH_indices$j,
    k = atlas_RH_indices$k
  ) %>%
    config(displaylogo = F, displayModeBar = F) %>%
    layout(scene = list(camera = list(
      eye = list(x = -1.5, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1),
      center = list(x = 0, y = 0, z = 0)
    ))) %>%
    layout(scene = list(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        zaxis = list(visible = FALSE)))

  suppressWarnings(orca(brain_list[[1]], "../data/brain_map/brain1.png", width = 1500, height = 1500))
  suppressWarnings(orca(brain_list[[2]], "../data/brain_map/brain2.png", width = 1500, height = 1500))
  suppressWarnings(orca(brain_list[[3]], "../data/brain_map/brain3.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[4]], "../data/brain_map/brain4.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[5]], "../data/brain_map/brain5.png", width = 1800, height = 1300))
  suppressWarnings(orca(brain_list[[6]], "../data/brain_map/brain6.png", width = 1800, height = 1300))


  brain_img_1 <- readPNG("../data/brain_map/brain1.png")
  brain_img_2 <- readPNG("../data/brain_map/brain2.png")
  brain_img_3 <- readPNG("../data/brain_map/brain3.png")
  brain_img_4 <- readPNG("../data/brain_map/brain4.png")
  brain_img_5 <- readPNG("../data/brain_map/brain5.png")
  brain_img_6 <- readPNG("../data/brain_map/brain6.png")

  brain_figs <- ggplot() +
    theme_classic() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.line=element_blank()) +
    annotation_custom(rasterGrob(brain_img_1, width = 1, height = 1),
                      xmin = 0, xmax = 0.294,
                      ymin = 0.5, ymax = 1) +
    annotation_custom(rasterGrob(brain_img_2, width = 1, height = 1),
                      xmin = 0, xmax = 0.294,
                      ymin = 0, ymax = 0.5) +
    annotation_custom(rasterGrob(brain_img_3, width = 1, height = 1),
                      xmin = 0.294, xmax = 0.647,
                      ymin = 0.5, ymax = 0.9333) +
    annotation_custom(rasterGrob(brain_img_4, width = 1, height = 1),
                      xmin = 0.294, xmax = 0.647,
                      ymin = 0.0667, ymax = 0.5) +
    annotation_custom(rasterGrob(brain_img_5, width = 1, height = 1),
                      xmin = 0.647, xmax = 1,
                      ymin = 0.5, ymax = 0.9333) +
    annotation_custom(rasterGrob(brain_img_6, width = 1, height = 1),
                      xmin = 0.647, xmax = 1,
                      ymin = 0.0667, ymax = 0.5) +
    coord_fixed(ratio = 0.5882)

  brain_map <- ggarrange(brain_figs, functional_term_rgb_legend, ncol = 2, nrow = 1, widths = c(0.85, 0.15))

  return(brain_map)
}


# Function generate a wordcloud based on Neurosynth data for a specific Yale Brain Atlas parcel
plot_neurosynth_wordcloud_for_parcel <- function(parcel, term_subset = "all", scale_height = 3, scale_width = 0.25){
  stopifnot("parcel must be a single character element containing the parcel name of one parcel in Yale Brain Atlas." =  (length(parcel) == 1) & (parcel %in% parcel_dict$Name))
  stopifnot("term_subset must be either 'all' or 'functional'" = (length(term_subset) == 1) & (term_subset %in% c("all", "functional")))
  stopifnot("scale_height must be a single numeric element." = (length(scale_height) == 1) & (is.numeric(scale_height)))
  stopifnot("scale_width must be a single numeric element." = (length(scale_width) == 1) & (is.numeric(scale_width)))

  if (term_subset == "all"){
    parcelsynth_df_for_parcel <- parcelsynth_df[which(parcel_dict$Name == parcel), ]
  } else if (term_subset == "functional"){
    parcelsynth_df_for_parcel <- parcelsynth_functional_df[which(parcel_dict$Name == parcel), ]
  }

  num_nonzero_terms <- rowSums(parcelsynth_df_for_parcel != 0)
  if (num_nonzero_terms > 100){
    num_nonzero_terms <- 100
  }

  if (num_nonzero_terms != 0){
    wordcloud(
      words = colnames(parcelsynth_df_for_parcel)[sort(as.numeric(as.vector(parcelsynth_df_for_parcel)), decreasing = T, index.return = T)$ix[1:num_nonzero_terms]],
      freq = 100*round(sort(as.numeric(as.vector(parcelsynth_df_for_parcel)), decreasing = T)[1:num_nonzero_terms], 2),
      scale = c(scale_height, scale_width),
      colors = brewer.pal(8, "Dark2"),
      random.order = FALSE
    )
  } else {
    print("No localized Neurosynth terms.")
  }
}


# Function to generate a barplot based on Neurosynth data for a specific Yale Brain Atlas parcel
plot_neurosynth_barplot_for_parcel <- function(parcel, term_subset = "all", plot_title = NULL){
  stopifnot("parcel must be a single character element containing the parcel name of one parcel in Yale Brain Atlas." =  (length(parcel) == 1) & (parcel %in% parcel_dict$Name))
  stopifnot("term_subset must be either 'all' or 'functional'" = (length(term_subset) == 1) & (term_subset %in% c("all", "functional")))
  stopifnot("plot_title must be a single character element." = is.null(plot_title) | (is.character(plot_title) & (length(plot_title) == 1)))

  if (term_subset == "all"){
    parcelsynth_df_for_parcel <- parcelsynth_df[which(parcel_dict$Name == parcel), ]
  } else if (term_subset == "functional"){
    parcelsynth_df_for_parcel <- parcelsynth_functional_df[which(parcel_dict$Name == parcel), ]
  }

  num_nonzero_terms <- as.numeric(rowSums(parcelsynth_df_for_parcel != 0))
  if (num_nonzero_terms > 30){
    num_nonzero_terms <- 30
  }

  if (num_nonzero_terms != 0){
    parcelsynth_df_for_parcel_barplot <- data.frame("term" = colnames(parcelsynth_df_for_parcel)[sort(as.numeric(as.vector(parcelsynth_df_for_parcel)), decreasing = T, index.return = T)$ix[1:num_nonzero_terms]],
                                                    "zscore" = sort(as.numeric(as.vector(parcelsynth_df_for_parcel)), decreasing = T)[1:num_nonzero_terms])
    neurosynth_barplot <-
      ggplot(
        data = parcelsynth_df_for_parcel_barplot, # highest 30
        aes(y = reorder(term, zscore), x = zscore)
      ) +
      xlab(paste("Z-score for Parcel", parcel)) +
      ylab(paste("Top Terms for Parcel", parcel)) +
      geom_bar(position = "stack", stat = "identity", fill = "#1f77b4") +
      ggtitle(plot_title) +
      theme_classic2() +
      theme(plot.title = element_text(hjust = 0.5))

  } else {
    neurosynth_barplot <- "No localized Neurosynth terms."
  }
  return(neurosynth_barplot)
}
