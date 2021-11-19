#!/bin/R

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)


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
###--- VERIFICATIONS OF INPUT
###---


if ( args$cell_type == FALSE && args$info == FALSE  ){ ##no argument given:
  ## gene_list will have the file_name by default
  parser$print_help()
  stop("EXIT!")

} else if (args$info) {

    merged_expr <-  read_delim("Merged_GeneExpr_per_CellType.tsv",delim="\t")

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

    merged_expr <-  read_delim("Merged_GeneExpr_per_CellType.tsv",delim="\t") #>> to change with name of "Merged" table of expr.values created with "analyse_sc.py"

    for (g in gene_name){
        tab_to_plot <- merged_expr[intersect(which(merged_expr$GeneName %in% gene_name),
                                             which(merged_expr$CellType %in% cell_type)) , ]

        tab_to_plot$TimePoint <- factor(tab_to_plot$TimePoint)
        tab_to_plot$TimePoint <- fct_relevel(tab_to_plot$TimePoint, c("D0","D2","D5","D7"))

        png(paste0(args$output_dir,"/Graph_of_Expr_",gene_name,"_SCdata.png" ))
            p <- ggplot(tab_to_plot, aes(x=TimePoint, y=MeanExpr, fill= TimePoint)) +
                  geom_bar(stat="identity", width=0.9, position="dodge") +
                  theme(axis.text.x=element_blank()) +
                  geom_errorbar(aes(ymin=MeanExpr-SEM, ymax=MeanExpr+SEM,color=TimePoint),
                                position="dodge",width=0.9,size=0.1)
            ## width of errorbar and geom_bar must be same,otherwise no overlap !
            ## width= {0.0, .., 1}
            p2 <- p + facet_wrap(~CellType) + ggtitle(g) + theme_minimal() + theme(axis.text.x=element_text(angle=45))
            print(p2)
        dev.off()

    #.....
    q("yes")
}
