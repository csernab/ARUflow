## ---------------------------
##
## Script name: fast_analysis.R
## Purpose of script:
## Author: Carlos Serna
## Date Created: 28/07/2020
## Email: carlsern@ucm.es
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# --- Packages --- #
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

# ========== #
# Input data #
# ========== #

# Resfinder report
resfinder_report <- read.csv(snakemake@input[["out_report_resfinder"]],sep = "\t",stringsAsFactors = FALSE, na.strings='.')
#resfinder_report <- read.csv("~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/extensive_test/summary_abricate_resfinder.tsv", sep = "\t",stringsAsFactors = FALSE,na.strings='.')
# Plasmidfinder report
plasmidfinder_report <- read.csv(snakemake@input[["out_report_plasmidfinder"]],sep = "\t",stringsAsFactors = FALSE, na.strings='.')
#plasmidfinder_report <- read.csv("~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/extensive_test/summary_abricate_plasmidfinder.tsv", sep = "\t",stringsAsFactors = FALSE,na.strings='.')


# ---- Clean data ----- #

# Lo voy a intentar meter en una función
clean_report <- function(input_report) {
  # Convert as character name
  input_report$X.FILE <- as.character(input_report$X.FILE)
  # Get basename (samples name)
  input_report$X.FILE <- sub('\\..*$', '', basename(input_report$X.FILE))
  ## Gene content
  input_report <- input_report[-2]
  # Prepare to plot
  melt.report <- melt(input_report, id.vars="X.FILE", variable_name="resistance_gene")
  colnames(melt.report) <- c("sample","resistance_gene","value")
  # Clean coverage values
  melt.report[] <- lapply(melt.report, function(x) gsub("100.00", "100", x))
  melt.report[] <- lapply(melt.report, function(x) gsub(";.*", "\\1", x))
  melt.report <- transform(melt.report, value = as.numeric(value))
  
  # Return new clean report
  return(melt.report)
}

# ---- Final reports ----- #

final_resfinder_report <- clean_report(resfinder_report)
final_plasmidfinder_report <- clean_report(plasmidfinder_report)


# ---- Resfinder analysis plot ----- #

ggplot(final_resfinder_report, aes(x = resistance_gene, y = sample, fill = value)) + geom_tile(colour = "black",width = 0.94,size = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 8, face = "italic", colour = "black", hjust = 0.1, vjust = 0.2),axis.line = element_line(color="white"), 
        axis.line.x = element_line(color="white"), axis.text.y = element_text(size = 10 , face = "bold", colour = "black")) +
  scale_x_discrete(position = "top") + scale_y_discrete(limits = unique(rev(final_resfinder_report$sample))) +
  ylab("") + xlab("") +
  scale_fill_continuous(na.value="#F5F5F5",low="thistle2", high="darkred", 
                        guide="colorbar", name = "Coverage")

ggsave(snakemake@output[["plot_resfinder"]], dpi=300, dev='png', height=15.5, width=16.5, units="in")


# ---- Plasmidfinder analysis plot ----- #

ggplot(final_plasmidfinder_report, aes(x = resistance_gene, y = sample, fill = value)) + geom_tile(colour = "black",width = 0.94,size = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 8, face = "italic", colour = "black", hjust = 0.1, vjust = 0.2),axis.line = element_line(color="white"), 
        axis.line.x = element_line(color="white"), axis.text.y = element_text(size = 10 , face = "bold", colour = "black")) +
  scale_x_discrete(position = "top") + scale_y_discrete(limits = unique(rev(final_plasmidfinder_report$sample))) +
  ylab("") + xlab("") +
  scale_fill_continuous(na.value="#F5F5F5",low="#8FBC8F", high="#228B22", 
                        guide="colorbar", name = "Coverage")

ggsave(snakemake@output[["plot_plasmidfinder"]], dpi=300, dev='png', height=15.5, width=16.5, units="in")


# ---- Gen and Plasmid count analysis plot ----- #

count_report <- function(input_report) {
  # Convert as character name
  input_report$X.FILE <- as.character(input_report$X.FILE)
  # Get basename (samples name)
  input_report$X.FILE <- sub('\\..*$', '', basename(input_report$X.FILE))
  ## Gene count
  gene_count <-  data.frame("sample" = input_report$X.FILE, "count" = input_report$NUM_FOUND)
  # Return count report
  return(gene_count)
}


# Resistance genes
gene_number <- count_report(resfinder_report)
# Plasmids
plasmid_number <- count_report(plasmidfinder_report)


p1 <- ggplot(gene_number, aes(x = sample, y = count)) +
  geom_bar(stat = "identity", width = 0.66, colour = "black", fill = "darkred") +
  theme_classic() + coord_flip() +
  scale_x_discrete(limits = unique(rev(gene_number$sample))) + 
  ggtitle("Resistance genes count")

p2 <- ggplot(plasmid_number, aes(x = sample, y = count)) +
  geom_bar(stat = "identity", width = 0.66, colour = "black", fill = "#228B22") +
  theme_classic() + coord_flip() +
  scale_x_discrete(limits = unique(rev(plasmid_number$sample))) +
  ggtitle("Plasmids count")

g <- grid.arrange(p1, p2, nrow = 1)

ggsave(snakemake@output[["plot_count"]], g,dpi=300, dev='png', height=10, width=12, units="in")
