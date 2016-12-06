if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tydir, ggplot2, cowplot, seqTools, substring)

library(dplyr)
library(ggplot2)
library(cowplot)
library(seqTools)
library(tidyr)

R1_PATTERN <- "Mcla_R1_subsample.fastq.gz"
R2_PATTERN <- "Mcla_R2_subsample.fastq.gz"

fqdir<-paste0(getwd(), "/data")
fqq<-fastqq(file.path(fqdir,
                      c(R1_PATTERN, R2_PATTERN)),
            k=4,probeLabel=c("g4","g5"))
# Set this to some other location
basedir<-getwd()
# get the sample name to print in stats file
sample_name <- strsplit(R1_PATTERN, "_R1")[[1]][[1]]

# output to directory
sink(file = paste0("data/", sample_name, ".readdata.txt"))
sample_name
fqq
sink()


phred_plot <- plotPhredQuant(fqq, 1)
phred_score <- phred(fqq, 1) # gives a count for all prhed scroes 
# -> can be used to plot with ggplot


# convert the output of seqTools to dataframe
# tidy up data 
data.frame(phred_score) %>%
  mutate(phred = row.names(.)) %>% 
  gather(., key = nt.pos.tmp, value = count, -phred) %>%
  separate(., col = nt.pos.tmp, into = c("scratch", "nt.pos"), sep= "X") %>%
  select(-scratch) -> dd
  #mutate(nt.pos = strsplit(nt.pos.tmp, "X")[[1]][2])-> dd

substring(x, 1, 1)
# plotting quality scores
ggplot(dd) %>%
  geom_boxplot(aes(x= nt.pos, y = count))

assemblies <- read.table(file = "assemblies.csv", sep= ",",
                         header = T, stringsAsFactors = F)

filename_from_path <- function(path) {
  # Takes a unix path "path/to/filename" string as input
  # splits string and returns
  # last element in path (the filename) as a string
  tail( # tail gives the last elements from lists
    strsplit(path, "/")[[1]], n=1)
}


# plot contig stats ------------------------------------------------------------
contigs <- read.table(file = "contigs.csv", sep= ",",
                      header = T, stringsAsFactors = F)

# cowplot provides better formatting
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
test_set <- contigs[1:1000,]



line <- ggplot(test_set) +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10() 

contigs %>% 
  select(coverage, contig_name) %>%
  transmute(cov = round(coverage/10, 0)*10) %>%
  group_by(cov) %>%
  summarise(N = n()) -> cov_count

cov_plot <- ggplot(data = cov_count, aes(cov, N, group = 1))

opt <- options("scipen" = 20) # don't use scinotation
# x and y log, look strange
cov_plot + geom_point() + 
  scale_x_log10(breaks = c(0,10, 100, 1000, 10000, 100000, 1000000)) +
  scale_y_log10(breaks =c( 1, 5,10,100,1000,10000)) +
  xlab("coverage") + ylab("No. contigs with this coverage")

# x and y log, look strange
cov_plot + geom_point() + 
  scale_x_log10(breaks = c(0,10, 100, 1000, 10000, 100000, 1000000)) +
  scale_y_continuous(breaks =c(100, 1000, 5000, 10000)) +
  xlab("coverage") + ylab("No. contigs with this coverage")

# histogram showing the distribution of coverage
ggplot(data= contigs) +
  geom_histogram(aes(round(coverage + 1 , 0)),
                 bins = 100,
                 fill="#00B0F6" ) + 
  # I add 1 here because the log of 0 is -inf and then the function won't plot it
  scale_x_log10(breaks = c(0,10, 100, 1000, 10000, 100000, 1000000)) +
  ggtitle("coverage distribution")

ggplot(data= contigs) +
  geom_histogram(aes(length),
                 bins = 100,
                 fill="#00B0F6" ) + 
  scale_x_log10(breaks=c(50,500,5000))
  # I add 1 here because the log of 0 is -inf and then the function won't plot it
  scale_x_log10(breaks = c(0,10, 50, 100, 1000, 1500, 2000, 5000)) +
  ggtitle("coverage distribution")

cov_plot + geom_line() + scale_x_log10() + scale_y_log10() +


ggplot(data= cov_count) +
  geom_histogram(aes(N))


opt <- options("scipen" = 20)

qplot(contigs$coverage, geom="histogram", log="yx", fill="blue") 

contig_plot +
  geom_point(aes(x= reorder(factor(contig_name), coverage), y=coverage)) +
  scale_y_log10()  






hist <- qplot(test_set$coverage, geom="histogram", log="yx", fill="blue") 

ggplot(test_set) +
  geom_histogram(aes(coverage), color = "blue", fill="red", binwidth = 1) +
  scale_x_log10() + scale_x_log10()


library(gridExtra)
grid.arrange(hist, line, phred_plot, ncol=2, nrow =2)

grid.arrange(hist, line, ncol=2, nrow =2)


# Create the external graphical elements
# called a "grop" in Grid terminology
phred_grob = ggplotGrob(phred_plot)
p3_grob = ggplotGrob(p3)

# cowplot also has a easier output function
save_plot("hist.test.pdf", hist)


plot_grid(hist, line)
plot_grid(hist, line, phred_plot, ncol = 1)

plot_grid(hist, line, line, labels = c("F", "U"))

plotA <- plot_grid(hist, line, line, ncol=1, labels = c("F", "U"))
# This does not work, we need to tell it about the layout
save_plot("dd.pdf", plotA)


save_plot("dd.pdf", plotA,
          nrow = 3,
          ncol = 1)

# probably the best wouldbe to output this to a file again and then 
# concatenate all the filesfrom the assembly or feed it into the 
# db.
assemblies %>% mutate(assembly_name = filename_from_path(assembly))

contigs <- read.table(file = "contigs.csv", sep= ",",
                      header = T, stringsAsFactors = F)


