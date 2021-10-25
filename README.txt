## TO RUN SCRIPT:

install R/4.0.2
Launch and install package :
#R>
> install.packages(c("argparse","ggplot2"))
> quit() --> y

## In terminal, in the directory where the tables-data and the script are:

$ create_histo_by_user_choice.R --info (or -i)  #to see the gene-names and cell-types among which you can choose.

#Then copy/paste in the following command the chosen names :
$ create_histo_by_user_choice.R -g  <..>   -c  <..>
