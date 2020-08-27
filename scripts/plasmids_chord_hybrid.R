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
library(circlize)
library(ggplot2)


#resfinder_folder <- "~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/abricate/resfinder/"
#plasmidfinder_folder <- "~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/abricate/plasmidfinder/"
#resfinder_folder <- "~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/extensive_hybrid_test/abricate/resfinder_hybrid/"
#plasmidfinder_folder <- "~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/extensive_hybrid_test/abricate/plasmidfinder_hybrid/"

# ========== #
# Input data #
# ========== #

resfinder_folder <- paste(snakemake@config[["project_dir"]],"/output_data/","/",snakemake@config[["project"]],"/abricate/resfinder_hybrid/", sep = "")
plasmidfinder_folder <- paste(snakemake@config[["project_dir"]],"/output_data/","/",snakemake@config[["project"]],"/abricate/plasmidfinder_hybrid/", sep = "")


#resfinder_folder <- paste(outdir_input,"/abricate/resfinder/", sep = "")
#plasmidfinder_folder <- paste(outdir_input,"/abricate/plasmidfinder/", sep = "")

# Get files 
files_resfinder <- list.files(path = resfinder_folder)
files_plasmidfinder <- list.files(path = plasmidfinder_folder)

# Chord diagrams
for (n_file in 1:length(files_resfinder)){
  
  resistance_genes <- read.csv(file = paste(resfinder_folder,files_resfinder[n_file],sep = ""), sep = "\t")
  plasmids <-  read.csv(file = paste(plasmidfinder_folder,files_plasmidfinder[n_file],sep = ""), sep = "\t")
  
  file_name <- sub('\\..*$', '', basename(files_plasmidfinder[n_file]))
  
  plasmid_name <- plasmids$GENE
  res_gene_name <- resistance_genes$GENE
  
  matrix_sample <- matrix(nrow = length(plasmid_name), ncol = length(res_gene_name))
  rownames(matrix_sample) <- plasmid_name
  colnames(matrix_sample) <- res_gene_name
  
  contig_intersect <- intersect(resistance_genes$SEQUENCE, plasmids$SEQUENCE)
  
  if (length(matrix_sample) == 0) {
    message("No plasmids or Resistance genes")
  } else {
    if (length(contig_intersect) == 0){
      message("No common contigs")
    } else {
      for (i in 1:length(resistance_genes$GENE)) {
        for (j in 1:length(plasmids$GENE)) {
          if (resistance_genes$SEQUENCE[i] == plasmids$SEQUENCE[j]) {
            matrix_sample[j,i] <- 1
          } else {
            matrix_sample[j,i] <- 0
          }
        }
      }
    }
  }
  
  if (length(matrix_sample) == 0) {
    message("No plasmids or Resistance genes")
  } else {
    if (length(contig_intersect) == 0){
      message("No common contigs")
    } else {
      png(paste(snakemake@output[["chord_folder"]],file_name,".png",sep = ""),units="in", width=5, height=5, res=300)
#      png(paste("~/Doctorado_laboratorio/TFM_bbc/Github_snakemake/results/extensive_hybrid_test/",file_name,".png",sep = ""),units="in", width=6, height=6, res=300)
      chordDiagram(matrix_sample, annotationTrack = "grid", preAllocateTracks = 1)
      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + .34, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
        circos.axis(h = "top", labels = TRUE,labels.cex = 0.5, major.tick.percentage = 0.1, 
                    sector.index = sector.name, track.index = 2, minor.ticks = 0)
      }, bg.border = NA)
                   
      dev.off()
    }
  }
}
