#!/bin/R

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)


###---
###--- LOCATE YOURSELF :
###---

work_dir <- "~/PROJECTS/Single_Cell_Analysis/"
input_dir <- paste0(work_dir, "Input/")
output_dir <- paste0(work_dir,"R_create_graph_per_gene/")

###---
###--- PARSE INPUT ARGUMENTS :
###---

parser <- ArgumentParser(epilog="type --info to see options for Input")

parser$add_argument("-i", "--info",action='store_true',default=FALSE,help="read the gene.names and cell.types accepted for analysis.")
parser$add_argument("-g", "--gene_list",action="store",default=paste0(input_dir,""),
                    help="Choose one GeneName from those given in list to visualise across cell-types.Default: all genes given in previous step will be used.")
parser$add_argument("-c","--cell_type", action="store", default=FALSE, help="Give cell-type(s)  for which to see expression of the gene")
parser$add_argument("-o","--output_dir", action="store", default=FALSE, help="Give the path of directory to save graphs.")


###---
###--- RETRIEVE VARIABLES :
###---

args <- parser$parse_args()
print(args)


###---
###--- TEST MY OWN VARIABLES FOR LAUNCHING SCRIPT
###---

input_gene_list <- paste0(input_dir,"DM1_geneset.tsv" )
# input_gene_list <- paste0(input_dir,"DM1_geneset1_DF.tsv" )
# input_gene_list <- paste0(input_dir,"Geneset1_CP.tsv" )

args$cell_type <- "Anti-inflammatory_macrophages,MuSCs_and_progenitors,Smooth_muscle_cells,Tenocytes"

args$output_dir <- paste0(work_dir,"R_create_graph_per_gene/Graphs_D0-D2-D7_per_gene")

dir.create(args$output_dir, showWarnings = FALSE)


###---
###--- VERIFICATIONS OF INPUT
###---


if ( args$cell_type == FALSE && args$info == FALSE  ){ ##no argument given:
  ## gene_list will have the file_name by default
  parser$print_help()
  stop("EXIT!")

} else if (args$info) {

    merged_expr <-  readRDS("Merged_GeneExpr_per_CellType.rds")

    cat("\n>>The gene_names from which to choose: "); print(levels(as.factor(merged_expr$GeneName)))
    cat("\n>>The cell-types from which to choose: "); print(levels(as.factor(merged_expr$CellType)))

    stop("EXIT!")

} else if ((! any(args$gene_name %in% merged_expr$GeneName)) ||  (length(args$cell_type) == 0 )) {

    stop("Wrong Input: \n
          Please fill in at least one gene-name but no >1 ,\n
          and at least one cell-type\n")

} else if (args$gene_name != FALSE && args$cell_type != FALSE ){

    gene_name <- args$gene_name
    cell_type <- args$cell_type


    tab_to_plot <- merged_expr[intersect(which(merged_expr$GeneName %in% gene_name),
                                         which(merged_expr$CellType %in% cell_type)) , ]


    pdf(paste0(output_dir,"Graph_of_Expr_",gene_name,"_SCdata.pdf" ))
      ggplot(tab_to_plot, aes(x =TimePoint, y=MeanExpr, fill=CellType)) +
            geom_bar(stat="identity", width=0.5,position="dodge") +
            ggtitle(gene_name) + theme_minimal()
            ##+geom_errorbar(aes(x =TimePoint, ymin=MeanExpr-STD, ymax=MeanExpr+STD),position="dodge",colour="grey",size=0.4)+
            ## error bars are HUGE !!
    dev.off()

    #.....
    q("yes")
}
