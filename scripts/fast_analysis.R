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
library(reshape)
library(ggplot2)


# ========== #
# Input data #
# ========== #

# Resfinder report
resfinder_report <- read.csv(snakemake@input[["out_report_resfinder"]])
# Plasmidfinder report
plasmidfinder_report <- read.csv(snakemake@input[["out_report_plasmidfinder"]])


# ---- Clean data ----- #

# Lo voy a intentar meter en una función
clean_report <- function(input_report) {
  # Convert as character name
  input_report$name <- as.character(input_report$name)
  # Get basename (samples name)
  input_report$name <- sub('\\..*$', '', basename(input_report$name))
  # Sample name vector
  sample <- input_report$name
  # Absence / presence vector
  report_content <- input_report[, !c(TRUE,FALSE)]
  # References
  references_content <- input_report[, !c(FALSE,TRUE)][-1]
  # empy vector
  references_vector <- character(0)
  # Get the reference names and accesions
  for (i in 1:ncol(references_content)) {
    if (dim(table(references_content[,i])) == 1) {
      ref <- unique(as.character(references_content[,i][!is.na(references_content[,i])]))
      ref <- gsub("*\\.[0-9]*\\_.*$","",ref)
      references_vector[i] <- ref 
    } else {
      ref <- unique(as.character(references_content[,i][!is.na(references_content[,i])]))
      ref <- gsub("*\\_.*$","", ref)
      ref <- paste(ref, collapse = " / ")
      references_vector[i] <- ref 
    }
  }
  # Change the column names in report content
  colnames(report_content) <- references_vector
  # Remove "_" final
  names(report_content) <- sub("\\_$","", names(report_content))
  # Add sample name column
  report_content["sample"] <- sample
  # Move last column to the start
  report_final <- report_content %>%
    select("sample", everything())
  
  # Return new clean report
  return(report_final)
}


# ---- Final reports ----- #

try(final_resfinder_report <- clean_report(resfinder_report))
try(final_plasmidfinder_report <- clean_report(plasmidfinder_report))


# ---- Plot ---- #


# RESFINDER
try(melt.resfinder <- melt(final_resfinder_report, id.vars="sample", variable_name="resistance_gene"))

try(qplot(data=melt.resfinder,
      x=resistance_gene,
      y=sample,
      fill=factor(value)) +
  geom_tile(colour = "black",width = 0.94,size = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 9, face = "italic", colour = "black", hjust = 0.1, vjust = 0.2),axis.line = element_line(color="white"), 
        axis.line.x = element_line(color="white"), axis.text.y = element_text(size = 10 , face = "bold", colour = "black")) +
  scale_x_discrete(position = "top") + scale_y_discrete(limits = unique(rev(melt.resfinder$sample))) +
  scale_fill_manual(values=c("no"="#F5F5F5", "yes"="#DC143C", "interrupted" = "#FF7F50","fragmented" = "#DAA520", 
                             "yes_nonunique" = "#800000", "partial" = "#FF8C00"), name = "") +
  ylab("") + xlab(""))

ggsave(snakemake@output[["plot_resfinder"]], dpi=300, dev='png', height=15, width=20, units="in")

# PLASMIDFINDER
try(melt.plasmidfinder <- melt(final_plasmidfinder_report, id.vars="sample", variable_name="plasmid"))

try(qplot(data=melt.plasmidfinder,
          x=plasmid,
          y=sample,
          fill=factor(value)) +
      geom_tile(colour = "black",width = 0.94,size = 0.4) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, size = 10, face = "italic", colour = "black",hjust = 0.1, vjust = 0.2),axis.line = element_line(color="white"), 
            axis.line.x = element_line(color="white"), axis.text.y = element_text(size = 10, face = "bold", colour = "black")) +
      scale_x_discrete(position = "top") + scale_y_discrete(limits = unique(rev(melt.plasmidfinder$sample))) +
      scale_fill_manual(values=c("no"="#F5F5F5", "yes"="#228B22","partial" = "#8FBC8F", 
                                 "yes_nonunique" = "#556B2F"), name = "") +
      ylab("") + xlab(""))

ggsave(snakemake@output[["plot_plasmidfinder"]], dpi=300, dev='png', height=15, width=14, units="in")
