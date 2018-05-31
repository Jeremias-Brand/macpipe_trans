
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('transrate/summary/Macmini2_20170608.contigs.csv'),
    output = list('Macmini2_20170608.png'),
    params = list(),
    wildcards = list('Macmini2', "sample" = 'Macmini2'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("transrate" = '/home/jeremias/soft/transrate-1.0.2-linux-x86_64/transrate', "trinity" = '/home/jeremias/soft/trinityrnaseq-2.2.0/Trinity', "adapters_fasta" = '/home/jeremias/soft/Trimmomatic-0.36/adapters/Truseq_barcodes.fa', "BUSCO_p" = '/home/jeremias/soft/busco/BUSCO_plot.py', "trim_minlength" = '25', "BUSCO" = '/home/jeremias/soft/busco/BUSCO.py', "trim_trailing" = '30', "py3" = '/home/jeremias/.linuxbrew/var/pyenv/versions/3.5.2/bin/python', "trim_leading" = '30', "BuscoLib" = '/home/jeremias/soft/busco/eukaryota_odb9', "fastqdir" = '/home/jeremias/dev/macpipe_trans/data/', "transrateDir" = '/home/jeremias/dev/macpipe_trans/transrate/', "homedir" = '/home/jeremias/dev/macpipe_trans/'),
    rule = 'plot_transrate_contigs'
)
######## Original script #########
<<<<<<< HEAD
#!/home/jeremias/.linuxbrew/bin/R
# load and if necessary install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, cowplot)

opt <- options("scipen" = 20) # don't use sci notation


# prepare plots ----------------------------------------------------------------
contigs <- read.table(file = "contigs.csv", sep= ",",
                      header = T, stringsAsFactors = F)

contigs <- contigs[1:1000,]

# histogram showing the distribution of coverage
cov_dist_plot <- ggplot(data= contigs) +
  geom_histogram(aes(round(coverage + 1 , 0)),
                 bins = 100,
                 fill="#00B0F6" ) + 
# I add 1 here because the log of 0 is -inf and then the function won't plot it
 scale_x_log10() +
  xlab("coverage") + 
  ylab("No. contigs") +
  ggtitle("contig coverage") + annotation_logticks(sides = "b")

# histogram of the length distribution
len_dist_plot <- ggplot(data= contigs) +
  geom_histogram(aes(length),
                 bins = 300,
                 fill="#00B0F6" ) + 
  scale_x_log10() +  
  annotation_logticks(sides = "b") +
  xlab("contig length") +
  ylab("No. contigs") +
  ggtitle("contig length")

# histogram of the length for contigs under 1000nt
len1k_dist_plot <- ggplot(data= filter(contigs, length < 1000)) +
  geom_histogram(aes(length),
                 bins = 100,
                 fill="#00B0F6" ) + 
  xlab("contig length") +
  ylab("No. contigs") +
  ggtitle("dist contigs < 1000nt")


# output plots -----------------------------------------------------------------

#TODO integrate Nex90N50 into this plot



arranged.plot <- plot_grid(cov_dist_plot, len_dist_plot, NULL,len1k_dist_plot,
                           ncol = 2, nrow =2 )
# cowplot does not have a title function
# So we can just use a grid to place the title on top
title <- ggdraw() + draw_label("Sample assembly statistic", fontface='bold')
out.plot <- plot_grid(title, arranged.plot,
          ncol=1, rel_heights=c(0.1, 1)) # rel height controls the title margin

# This does not work, we need to tell it about the layout
save_plot("Sample.png", out.plot, base_width = 8, base_height = 8)




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
=======
#!/home/jeremias/.linuxbrew/bin/R
# load and if necessary install packages
# We need to expicitly tell R to use cairo for the server.
options(bitmapType='cairo')
X11.options(type="cairo")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, cowplot)


opt <- options("scipen" = 20) # don't use sci notation

# Snakemake provides the piplein imput with the snakemake@ command
# R is counting from 1!
# load the input file
input_file = snakemake@input[[1]]
# get the sample name from the input filename
filename <- tail(strsplit(input_file, "/")[[1]], 1)
sample_name <- strsplit(filename, ".contigs.csv")[[1]]


# prepare plots ----------------------------------------------------------------
contigs <- read.table(file = input_file, sep= ",",
                      header = T, stringsAsFactors = F)

# histogram showing the distribution of coverage
cov_dist_plot <- ggplot(data= contigs) +
  geom_histogram(aes(round(coverage + 1 , 0)),
                 bins = 100,
                 fill="#00B0F6" ) + 
# I add 1 here because the log of 0 is -inf and then the function won't plot it
 scale_x_log10() +
  xlab("coverage") + 
  ylab("No. contigs") +
  ggtitle("contig coverage") + annotation_logticks(sides = "b")

# histogram of the length distribution
len_dist_plot <- ggplot(data= contigs) +
  geom_histogram(aes(length),
                 bins = 300,
                 fill="#00B0F6" ) + 
  scale_x_log10() +  
  annotation_logticks(sides = "b") +
  xlab("contig length") +
  ylab("No. contigs") +
  ggtitle("contig length")

# histogram of the length for contigs under 1000nt
len1k_dist_plot <- ggplot(data= filter(contigs, length < 1000)) +
  geom_histogram(aes(length),
                 bins = 100,
                 fill="#00B0F6" ) + 
  xlab("contig length") +
  ylab("No. contigs") +
  ggtitle("dist contigs < 1000nt")


# output plots -----------------------------------------------------------------

#TODO integrate Nex90N50 into this plot
arranged.plot <- plot_grid(cov_dist_plot, len_dist_plot, NULL,len1k_dist_plot,
                           ncol = 2, nrow =2 )
# cowplot does not have a title function
# So we can just use a grid to place the title on top
title <- ggdraw() + draw_label(paste0(sample_name, "assembly statistics"), fontface='bold')
out.plot <- plot_grid(title, arranged.plot,
          ncol=1, rel_heights=c(0.1, 1)) # rel height controls the title margin

# This does not work, we need to tell it about the layout
save_plot(paste0("transrate/summary/", sample_name, ".png"), out.plot, base_width = 8, base_height = 8)
>>>>>>> origin/master
