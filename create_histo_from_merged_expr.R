#!/bin/R


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


###---
###--- LOCATE YOURSELF :
###---

work_dir <- "/projects/single_cell_skeletal_muscle/"
input_dir <- paste0(work_dir, "Input/")
output_dir <- paste0(work_dir,"R_graphs/")


###---
###--- GIVE INPUT VARIABLES
###---

args <- list()
args$gene_list <- paste0(input_dir,"my_genes.tsv" )

args$cell_type <- c("Anti-inflammatory_macrophages","B_T_NK_cells",
                    "Endothelial","FAPs","Mature_skeletal_muscle","Monocytes_Macrophages_Platelets",
                    "MuSCs_and_progenitors","Neural_Glial_Schwann_cells","Pro-inflammatory_macrophages",
                    "Resident_Macrophages_APCs","Smooth_muscle_cells","Tenocytes")
##> terms should be as given from metadata file of published S.C. data

args$output_dir <- paste0(output_dir,"my_GoI_graphs/")

merged_expr <-  read_delim(paste0(work_dir,"Output/Merged/Merged_MeanExpr+sem_All_CellTypes.tsv"),
                           delim="\t")


### READ INPUT

dir.create(args$output_dir, showWarnings = FALSE)
gene_table <- read.table(args$gene_list,sep="\t",header=TRUE,as.is=T)
gene_name <- gene_table$mm_gene_name
cell_type <- args$cell_type

###
### PLOT in facets per Cell-type
###

for(g in gene_name){
  if (g %in% levels(as.factor(merged_expr$GeneName))) {
      tab_to_plot <- merged_expr[intersect(which(merged_expr$GeneName == g ),
                                           which(merged_expr$CellType %in% cell_type)), ]

      png(paste0(args$output_dir,"Graph_of_Expr_",g,"_SCdata.png"))
          p <- ggplot(tab_to_plot, aes(x=TimePoint, y=MeanExpr, fill=TimePoint)) +
            geom_bar(stat="identity", width=0.9, position="dodge") +
            geom_errorbar(aes(ymin=MeanExpr-SEM, ymax=MeanExpr+SEM,color=TimePoint),
                          position="dodge",width=0.9,size=0.3)
          p2 <- p + facet_wrap(~CellType) + ggtitle(g) + theme_minimal()
          print(p2)
      dev.off()
  }
}
