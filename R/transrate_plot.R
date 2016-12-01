# R script to plot and tidy transrate output
# Author: Jeremias Brand, University of Basel
# 2016

# Transrate produces two main output file that should be plotted:
# assembles.csv gives stats on the transrate score
# contigs.csv provides per contig details on mapping and coverage
library(dplyr)
library(ggplot2)
assemblies <- read.table(file = "assemblies.csv", sep= ",",
                         header = T, stringsAsFactors = F)

filename_from_path <- function(path) {
  # Takes a unix path "path/to/filename" string as input
  # splits string and returns
  # last element in path (the filename) as a string
  tail( # tail gives the last elements from lists
    strsplit(path, "/")[[1]], n=1)
  }


# probably the best wouldbe to output this to a file again and then 
# concatenate all the filesfrom the assembly or feed it into the 
# db.
assemblies %>% mutate(assembly_name = filename_from_path(assembly))

contigs <- read.table(file = "contigs.csv", sep= ",",
                         header = T, stringsAsFactors = F)

contigs <- arrange(contigs, coverage)

contig_plot <- ggplot(data = contigs)

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10()


contigs %>% 
  select(coverage, contig_name) %>%
  transmute(cov = round(coverage/10, 0)) %>%
  group_by(cov) %>%
  summarise(N = n()) -> cov_count

cov_plot <- ggplot(data = cov_count, aes(cov, N, group = 1))

cov_plot + geom_point() + scale_x_log10()

ggplot(data= contigs) +
  geom_histogram(coverage)
opt <- options("scipen" = 20)

qplot(contigs$coverage, geom="histogram", log="yx", fill="blue") 

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10()  