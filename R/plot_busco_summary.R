# Read only the relevant line of busco summary file
# setwd(snakemakefile)
#options(bitmapType='cairo')
#X11.options(type="cairo")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, scales)


# input_file = 
busco_folder = snakemake@input[[1]] # path to the place where the busco files are stored
output_folder = busco_folder
TIMESTAMP = snakemake@params[[1]]
plot_basename = paste0("run_", TIMESTAMP, "_busco")

###############################################################################
# choose a pattern that fits the summaries you want to plot
PATTERN = "short_summary_Mac[A-Za-z0-9]*.+.busco.txt" # most inclusive
# PATTERN = "short_summary_Mac[A-Za-z0-9]*_*[0-9]*.busco.txt" # only with ID
summaries <- list.files(path = busco_folder, pattern = PATTERN)
# all the file names in the folder
sink("debug.file_plot_busco_summary")
sink()
my_output <- paste("busco/summary/","busco_figure.png",sep="/") 
my_width <- 20
my_height <- 15
my_unit <- "cm"
my_bar_height <- 0.75

pdf_width = 15
pdf_height = 10 + length(summaries) * 0.25
# Colors
my_colors <- c("#F04442", "#F0E442", "#3492C7", "#56B4E9")

# Bar height ratio

# Legend
my_title <- "BUSCO assessment
 results of SmartSeqv4 Trinity assemblies from run P428 \n
 used eukaryote db with 303 genes"

# Font
my_family <- "sans"
my_size_ratio <- 1

###############################################################################
# Setup busco results table template
n_runs <- length(summaries)
cols <- 7
busco_table <- data.frame(matrix(nrow = n_runs, ncol = cols))
names(busco_table) <- c("assembly", "complete",
                        "complete_single", "complete_dublicated", "fragmented",
                        "missing", "n")
###############################################################################
# function to get results from a file
extract_busco_results <- function(filename) {
  # input: busco short summary file from BUSCO 2.0 run
  # output: vector with assembly name followed by complete, single, duplicated,
  # fragemented, missing, n
  con <- file(filename)
  lines <- readLines(con = con, n = -1L)
  close(con)
  
  adline   <- lines[[5]]
  path     <- strsplit(adline, split = " ")[[1]][[9]]
  assembly <- tail(strsplit(path, split = "/")[[1]], n = 1)
  
  resline  <- lines[[8]]
  resline  <- gsub(pattern = "[^0-9 .]", replacement = " " , x = resline)
  resline  <- unlist(strsplit(resline, split = " +"))
  resline[1]  <- assembly
  return(resline)
}


###############################################################################
# gather all the results
for (i in 1 : length(summaries)){
  busco_res <- extract_busco_results(paste0(busco_folder, summaries[[i]]))
  for (res in 1 :length(busco_res)) {
    busco_table[i, res] = busco_res[res]
  }
}
tidy_busco_table <- gather(busco_table, key = category, value = percentage, -assembly)
tidy_busco_table <- filter(tidy_busco_table, category != "n", category != "complete")
tidy_busco_table$category <- factor(tidy_busco_table$category, 
                                    levels = c('missing', 'fragmented', 'complete_dublicated', 'complete_single'))
tidy_busco_table <- mutate(tidy_busco_table, assembly_name = gsub(".Trinity.fasta$","", assembly))


busco_plot <- ggplot(tidy_busco_table, aes(y = as.numeric(percentage), x = assembly_name, fill = category)) +
  geom_bar(position = "fill", stat="identity") +
  scale_y_continuous(labels = percent_format()) + 
  coord_flip() + 
  theme_gray(base_size = 8) + 
  scale_fill_manual(values = my_colors,labels =c(" Missing (M)",
                                                 " Fragmented (F)  ",
                                                 " Complete (C) and duplicated (D)",
                                                 " Complete (C) and single-copy (S)  ")) +
  ggtitle(my_title) + 
  xlab("") + 
  ylab("\n%BUSCOs") + 
  theme(plot.title = element_text(family=my_family, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
  theme(legend.position="top",legend.title = element_blank()) + 
  theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) + 
  theme(axis.ticks.length = unit(.85, "cm")) + 
  theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
  theme(axis.ticks.x = element_line(colour="#222222")) + 
  theme(axis.ticks.length = unit(0.4, "cm")) + 
  theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

pdf(file = paste(output_folder, paste0(plot_basename, ".pdf"), sep="/"),
    width = pdf_width, height = pdf_height)
busco_plot
dev.off()


# this part also adds transrate plot 

# Now I also want to add the results from the transrate run.
# I generated this by running a mini snakemake script called quick_gater.snake
# it allows to summarise all the assembly files in a folder
# transrate_assembly_summary <- read.csv("transrate/run_2017045_assemblies.csv")
# 
# # I clean the name to extract just the identifier of a assembly
# transrate_assembly_summary <- transrate_assembly_summary %>%
#   mutate(tmp_name = gsub(".Trinity.fasta$","", assembly)) %>%
#   mutate(tmp_name2 = gsub("^/.+/.+/.+/.+/","", tmp_name)) %>%
#   mutate(assembly_name = as.factor(tmp_name2)) %>%
#   select(-tmp_name, -tmp_name2) 
# 
# # I do the same thing for the busco data to have a clean dataset
# tidy_busco_table_wide <- busco_table %>%
#   mutate(tmp_name = gsub(".Trinity.fasta$","", assembly)) %>%
#   mutate(assembly_name = as.factor(tmp_name)) %>%
#   select(-tmp_name) 
# 
# # Join the two keeping only the ones for which I have a transrate score
# joined_busco_table <- full_join(transrate_assembly_summary, tidy_busco_table_wide, by = "assembly_name")
# 
# # Combine the name with the scores to plot them on the busco table
# joined_busco_table <- joined_busco_table %>%
#   # roundsocre for readability
#   mutate(score = round(score,3)) %>%
#   # merge columns, carefull first argument is name of new column
#   unite(assembly_name_stats, assembly_name, n_seqs, score, sep= " ") %>%
#   select(assembly_name_stats, complete_single,
#          complete_dublicated, fragmented, missing) %>%
#   # put in long form
#   gather(key = category, value = percentage, -assembly_name_stats) %>%
#   mutate(category = factor(category, 
#                 levels = c('missing', 'fragmented', 'complete_dublicated', 'complete_single')))
#          
# # Plotting again
# 
# busco_plot_score <- ggplot(joined_busco_table, aes(y = as.numeric(percentage), x = assembly_name_stats, fill = category)) +
#   geom_bar(position = "fill", stat="identity") +
#   scale_y_continuous(labels = percent_format()) + 
#   coord_flip() + 
#   theme_gray(base_size = 8) + 
#   scale_fill_manual(values = my_colors,labels =c(" Missing (M)",
#                                                  " Fragmented (F)  ",
#                                                  " Complete (C) and duplicated (D)",
#                                                  " Complete (C) and single-copy (S)  ")) +
#   ggtitle(my_title) + 
#   xlab("") + 
#   ylab("\n%BUSCOs") + 
#   theme(plot.title = element_text(family=my_family, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
#   theme(legend.position="top",legend.title = element_blank()) + 
#   theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
#   theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
#   theme(panel.grid.minor = element_blank()) + 
#   theme(panel.grid.major = element_blank()) +
#   theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
#   theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
#   theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) + 
#   theme(axis.ticks.length = unit(.85, "cm")) + 
#   theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
#   theme(axis.ticks.x = element_line(colour="#222222")) + 
#   theme(axis.ticks.length = unit(0.4, "cm")) + 
#   theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
#   guides(fill = guide_legend(override.aes = list(colour = NULL))) +
#   guides(fill=guide_legend(nrow=2,byrow=TRUE))
# 
# pdf("busco_single_pooled_transrate_score.pdf", width = 15, height = 35)
# busco_plot_score
# dev.off()
# 
