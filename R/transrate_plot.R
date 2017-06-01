<<<<<<< HEAD
# R script to plot and tidy transrate output
# Author: Jeremias Brand, University of Basel
# 2016

# Transrate produces two main output file that should be plotted:
# assembles.csv gives stats on the transrate score
# contigs.csv provides per contig details on mapping and coverage

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, cowplot)

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

# cowplot provides better formatting
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
test_set <- contigs[1:1000,]

line <- ggplot(test_set) +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10() + scale_x_log10()

hist <- qplot(test_set$coverage, geom="histogram", log="yx", fill="blue") 

ggplot(test_set) +
  geom_histogram(aes(coverage), color = "blue", fill="red", binwidth = 1) +
  scale_x_log10() + scale_x_log10()

# cowplot also has a easier output function
save_plot("hist.test.pdf", hist)


plot_grid(hist, line)
plot_grid(hist, line, line, labels = c("F", "U"))

plotA <- plot_grid(hist, line, line, ncol=1, labels = c("F", "U"))
# This does not work, we need to tell it about the layout
save_plot("dd.pdf", plotA)


save_plot("dd.pdf", plotA,
          nrow = 3,
          ncol = 1)




contigs <- arrange(contigs, coverage)

contig_plot <- ggplot(data = contigs)

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10()


contigs %>% 
  select(coverage, contig_name) %>%
  transmute(cov = round(coverage/10, 0)*10) %>%
  group_by(cov) %>%
  summarise(N = n()) -> cov_count

cov_plot <- ggplot(data = cov_count, aes(cov, N, group = 1))

cov_plot + geom_point() + scale_x_log10()

cov_plot + geom_point() + scale_x_log10()



ggplot(data= contigs) +
  geom_histogram(round(coverage, 0))
opt <- options("scipen" = 20)

qplot(contigs$coverage, geom="histogram", log="yx", fill="blue") 

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
=======
# R script to plot and tidy transrate output
# Author: Jeremias Brand, University of Basel
# 2016

# Transrate produces two main output file that should be plotted:
# assembles.csv gives stats on the transrate score
# contigs.csv provides per contig details on mapping and coverage

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, cowplot)

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

# cowplot provides better formatting
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
test_set <- contigs[1:1000,]

line <- ggplot(test_set) +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10() + scale_x_log10()

hist <- qplot(test_set$coverage, geom="histogram", log="yx", fill="blue") 

ggplot(test_set) +
  geom_histogram(aes(coverage), color = "blue", fill="red", binwidth = 1) +
  scale_x_log10() + scale_x_log10()

# cowplot also has a easier output function
save_plot("hist.test.pdf", hist)


plot_grid(hist, line)
plot_grid(hist, line, line, labels = c("F", "U"))

plotA <- plot_grid(hist, line, line, ncol=1, labels = c("F", "U"))
# This does not work, we need to tell it about the layout
save_plot("dd.pdf", plotA)


save_plot("dd.pdf", plotA,
          nrow = 3,
          ncol = 1)




contigs <- arrange(contigs, coverage)

contig_plot <- ggplot(data = contigs)

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10()


contigs %>% 
  select(coverage, contig_name) %>%
  transmute(cov = round(coverage/10, 0)*10) %>%
  group_by(cov) %>%
  summarise(N = n()) -> cov_count

cov_plot <- ggplot(data = cov_count, aes(cov, N, group = 1))

cov_plot + geom_point() + scale_x_log10()

cov_plot + geom_point() + scale_x_log10()



ggplot(data= contigs) +
  geom_histogram(round(coverage, 0))
opt <- options("scipen" = 20)

qplot(contigs$coverage, geom="histogram", log="yx", fill="blue") 

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
>>>>>>> origin/master
  scale_y_log10()  