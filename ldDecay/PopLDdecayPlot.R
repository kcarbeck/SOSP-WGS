## 28 August 2023
## katherine carbeck
## Plot LD decay output from PopLDdecay

####! step 3: better plots in R
###! using *.bin output from PopLddecay


library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)

setwd("/workdir/kcarbeck/ldDecay/out300kb/output")

# Set path
my_path <- "/workdir/kcarbeck/ldDecay/out300kb/output"  # update dir path if neccessary (local vs terimal plotting)

# get list of file names ending with "_ldDecay.stat.gz"
file_list <- list.files(path = my_path, pattern = ".bin$", full.names = TRUE)

# loop through each file, read data, create plot and then save it
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist = V1, mean_r2 = V2, sum_r2 = V4, n_pairs = V6) %>%
    dplyr::select(-V3, -V5) 

  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub(".bin", "", file_name)
  
  # Add a new column for file_name
  ld_data$file_name <- file_name

  # Create plot
  p <- ld_data %>%
    dplyr::filter(dist <= 5000) %>%
    ggplot(aes(x = dist, y = mean_r2)) +
    geom_line(size = 0.7) +
    scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
    facet_wrap(~ file_name, ncol = 4) +  # Create facets for each population
    theme_bw() +
    theme(axis.title.x = element_text(size = 11, color = 'black'), 
          axis.title.y = element_text(size = 11, color = 'black'), 
          axis.text.y = element_text(size = 10, color = 'black'), 
          axis.text.x = element_text(size = 10, color = 'black'), 
          legend.position = c(0.75, 0.85), panel.grid = element_blank(),
          legend.text = element_text(size = 10, color = 'black'),
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.08, "in")) +
    xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
  
  # Save
  ggsave(paste0("plots/", file_name, "_ldDecay_plot.png"), plot = p, width = 8, height = 7)
}


# imgcat /workdir/kcarbeck/ldDecay/out300kb/output/plots/adusta_MX_ldDecay_plot.png




###! combined plot

# Create a list to store data and population names
ld_data_list <- list()

# Loop through each file, read data, and store it in the list
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist=V1, mean_r2=V2, sum_r2=V4, n_pairs=V6) %>%
    dplyr::select(-V3,-V5) 
  
  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub(".bin", "", file_name)
  
  # Store data and population name in the list
  ld_data_list[[file_name]] <- ld_data
}

# Combine all data frames in the list into a single data frame
combined_data <- bind_rows(ld_data_list, .id = "Population")

##! LINE
# Create the plot with colored lines for each population
main_plot <- combined_data %>%
  dplyr::filter(dist <= 5000) %>%
  #dplyr::filter(bins <= 200000) %>%
  ggplot(aes(x = dist, y = mean_r2, color = Population)) +
  geom_line(size = 0.7) +
  # stat_smooth(method = "loess",se=F,size=0.7, span=0.8) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

ggsave("plots/combined_LdDecay_Plot.png", plot = main_plot, width = 9, height = 6)

# view in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/output/plots/combined_LdDecay_Plot.png

##! LOESS
# Create the plot with colored lines for each population
main_plot <- combined_data %>%
  dplyr::filter(dist <= 5000) %>%
  #dplyr::filter(bins <= 200000) %>%
  ggplot(aes(x = dist, y = mean_r2, color = Population)) +
  # geom_line(size = 0.7) +
  stat_smooth(method = "loess",se=F,size=0.7, span=0.75) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

ggsave("plots/combined_Loess_LdDecay_Plot.png", plot = main_plot, width = 9, height = 6)

# view in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/output/plots/combined_Loess_LdDecay_Plot.png



###! faceted plot
# Create a list to store data and population names
ld_data_list <- list()

# Loop through each file, read data, and store it in the list
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist=V1, mean_r2=V2, sum_r2=V4, n_pairs=V6) %>%
    dplyr::select(-V3,-V5) 
  
  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub(".bin", "", file_name)
  
  # Store data and population name in the list
  ld_data_list[[file_name]] <- ld_data
}

# Combine all data frames in the list into a single data frame
combined_data <- bind_rows(ld_data_list, .id = "Population")

# Create an empty ggplot object
p <- ggplot() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))

# Loop through each file, read data, and add to the ggplot object
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist = V1, mean_r2 = V2, sum_r2 = V4, n_pairs = V6) %>%
    dplyr::select(-V3, -V5)  %>%
    dplyr::filter(dist <= 7500)

  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub(".bin", "", file_name)
  
  # Add a new column for file_name
  ld_data$file_name <- file_name

  # Add to the existing ggplot object
  p <- p +
    # stat_smooth(data = ld_data, aes(x = dist / 1e3, y = mean_r2),
    stat_smooth(data = ld_data, aes(x = dist, y = mean_r2),
        method = "loess",se=F,size=0.7, span=0.5) +
    # geom_line(data = ld_data, aes(x = dist, y = mean_r2), size = 0.7) +
    facet_wrap(~ file_name, ncol = 4)
}

# Set y-axis limits and save the plot
p <- p + scale_y_continuous(limits = c(0, .7))
ggsave("plots/faceted_Loess_ldDecay_plot.png", plot = p, width = 9, height = 12)


# view in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/output/plots/faceted_Loess_ldDecay_plot.png



















#########!

#########!

#########!    with .stat.gz output

#########!

#########!


library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)

setwd("/workdir/kcarbeck/ldDecay/out300kb")

# Set path
my_path <- "/workdir/kcarbeck/ldDecay/out300kb"  # update dir path if neccessary (local vs terimal plotting)

# get list of file names ending with "_ldDecay.stat.gz"
file_list <- list.files(path = my_path, pattern = "_ldDecay.stat.gz$", full.names = TRUE)

# loop through each file, read data, create plot and then save it
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist = V1, mean_r2 = V2, sum_r2 = V4, n_pairs = V6) %>%
    dplyr::select(-V3, -V5) 

  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub("_ldDecay.stat.gz", "", file_name)
  
  # Add a new column for file_name
  ld_data$file_name <- file_name

  # Create plot
  p <- ld_data %>%
    dplyr::filter(dist <= 5000) %>%
    ggplot(aes(x = dist, y = mean_r2)) +
    geom_line(size = 0.7) +
    scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
    facet_wrap(~ file_name, ncol = 4) +  # Create facets for each population
    theme_bw() +
    theme(axis.title.x = element_text(size = 11, color = 'black'), 
          axis.title.y = element_text(size = 11, color = 'black'), 
          axis.text.y = element_text(size = 10, color = 'black'), 
          axis.text.x = element_text(size = 10, color = 'black'), 
          legend.position = c(0.75, 0.85), panel.grid = element_blank(),
          legend.text = element_text(size = 10, color = 'black'),
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.08, "in")) +
    xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
  
  # Save
  ggsave(paste0("stats.plots/", file_name, "_ldDecay_plot.png"), plot = p, width = 8, height = 7)
}

# view plots in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/adusta_MX_ldDecay_plot.png




###! combined plot

# Create a list to store data and population names
ld_data_list <- list()

# Loop through each file, read data, and store it in the list
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist=V1, mean_r2=V2, sum_r2=V4, n_pairs=V6) %>%
    dplyr::select(-V3,-V5) 
  
  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub("_ldDecay.stat.gz", "", file_name)
  
  # Store data and population name in the list
  ld_data_list[[file_name]] <- ld_data
}

# Combine all data frames in the list into a single data frame
combined_data <- bind_rows(ld_data_list, .id = "Population")




##! LINE
# Create the plot with colored lines for each population
main_plot <- combined_data %>%
  dplyr::filter(dist <= 3000) %>%
  #dplyr::filter(bins <= 200000) %>%
  ggplot(aes(x = dist, y = mean_r2, color = Population)) +
  geom_line(size = 0.7) +
  # stat_smooth(method = "loess",se=F,size=0.7, span=0.8) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

ggsave("stats.plots/combined_3000bp_LdDecay_Plot.png", plot = main_plot, width = 9, height = 6)

# view in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/combined_5000bp_LdDecay_Plot.png
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/combined_3000bp_LdDecay_Plot.png





##! LOESS
# Create the plot with colored lines for each population
main_plot <- combined_data %>%
  dplyr::filter(dist <= 5000) %>%
  #dplyr::filter(bins <= 200000) %>%
  ggplot(aes(x = dist, y = mean_r2, color = Population)) +
  # geom_line(size = 0.7) +
  stat_smooth(method = "loess",se=F,size=0.7, span=0.2) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

ggsave("stats.plots/combined_Loess_LdDecay_Plot.png", plot = main_plot, width = 9, height = 6)

# view in terminal
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/combined_Loess_LdDecay_Plot.png







###! faceted plot
# Create an empty ggplot object
p <- ggplot() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 11, color = 'black'), 
        axis.title.y = element_text(size = 11, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        panel.grid = element_blank()) +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))

# Loop through each file, read data, and add to the ggplot object
for (file in file_list) {
  ld_data <- read.table(file) %>%
    dplyr::rename(dist = V1, mean_r2 = V2, sum_r2 = V4, n_pairs = V6) %>%
    dplyr::select(-V3, -V5)  %>%
    dplyr::filter(dist <= 2000)

  # Extract file name w/o path and extension
  file_name <- gsub(".*/", "", file)
  file_name <- gsub("_ldDecay.stat.gz", "", file_name)
  
  # Add a new column for file_name
  ld_data$file_name <- file_name

  # Add to the existing ggplot object
  p <- p +
    # stat_smooth(data = ld_data, aes(x = dist / 1e3, y = mean_r2),
    # stat_smooth(data = ld_data, aes(x = dist, y = mean_r2),
    #             method = "loess",se=F,size=0.7, span=0.2) +
    geom_line(data = ld_data, aes(x = dist, y = mean_r2), size = 0.7) +
    facet_wrap(~ file_name, ncol = 4)
}

# Set y-axis limits and save the plot
p <- p + scale_y_continuous(limits = c(0, .8))
#ggsave("stats.plots/faceted_Loess_ldDecay_plot.png", plot = p, width = 9, height = 12)
# ggsave("stats.plots/faceted_Line_ldDecay_plot.png", plot = p, width = 9, height = 12)
ggsave("stats.plots/faceted_2000bpLine_ldDecay_plot.png", plot = p, width = 9, height = 12)


# view in terminal
# LOESS
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/faceted_Loess_ldDecay_plot.png
# LINE
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/faceted_Line_ldDecay_plot.png
# LINE 2000bp
imgcat /workdir/kcarbeck/ldDecay/out300kb/stats.plots/faceted_2000bpLine_ldDecay_plot.png


