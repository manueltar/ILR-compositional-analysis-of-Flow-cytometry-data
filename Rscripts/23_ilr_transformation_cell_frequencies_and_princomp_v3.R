
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
# suppressMessages(library("ggforce", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)


ilr_transformation = function(option_list)
{
  suppressMessages(library("compositions", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Master_path_analysis----
  
  Master_path_analysis = opt$Master_path_analysis
  
  cat("Master_path_analysis\n")
  cat(sprintf(as.character(Master_path_analysis)))
  cat("\n")
 
 
  #### Read FlowCyt_results_subset----
  
  FlowCyt_results_subset<-readRDS(file=opt$FlowCyt_results_subset)
  
  cat("FlowCyt_results_subset_0\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FlowCyt_results_subset$treatment))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FlowCyt_results_subset$treatment)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FlowCyt_results_subset$Sample_label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FlowCyt_results_subset$Sample_label)))))
  cat("\n")
  
  
  
  ##### create matrix samples by cell frequencies ----
  
  indx.int<-c(which(colnames(FlowCyt_results_subset) == 'Double_neg'),which(colnames(FlowCyt_results_subset) == 'CD235_single'),which(colnames(FlowCyt_results_subset) == 'Double_pos'),which(colnames(FlowCyt_results_subset) == 'CD41_single'))
  
  FlowCyt_results_subset_matrix<-as.matrix(FlowCyt_results_subset[,indx.int])
  row.names(FlowCyt_results_subset_matrix)<-FlowCyt_results_subset$Sample_label
  
  cat("FlowCyt_results_subset_matrix_0\n")
  cat(str(FlowCyt_results_subset_matrix))
  cat("\n")
  
  ilr_object<-as.data.frame(ilr(FlowCyt_results_subset_matrix),stringsAsFactors=F)
  ilr_object$Sample_label<-row.names(ilr_object)
  row.names(ilr_object)<-NULL
  
  indx.ilr<-which(colnames(ilr_object) != 'Sample_label')
  
  colnames(ilr_object)[indx.ilr]<-paste(colnames(ilr_object)[indx.ilr],'_ilr_value',sep='')
  
  
  cat("ilr_object_0\n")
  cat(str(ilr_object))
  cat("\n")
  
  ##### merge ----
  
  FlowCyt_results_subset<-merge(FlowCyt_results_subset,
                                ilr_object,
                                by='Sample_label')
  
  cat("FlowCyt_results_subset_1\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
 
  
  ################    SAVE   #######################
  
  setwd(out)
  
  
  write.table(FlowCyt_results_subset, file='FlowCyt_global_adapted_plus_ilr.csv', sep=",", quote=F, row.names = F)
  saveRDS(FlowCyt_results_subset, file='FlowCyt_global_adapted_plus_ilr.rds')
  
}

princomp_operation = function(option_list)
{
  suppressMessages(library("compositions", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Master_path_analysis----
  
  Master_path_analysis = opt$Master_path_analysis
  
  cat("Master_path_analysis\n")
  cat(sprintf(as.character(Master_path_analysis)))
  cat("\n")
  
  
  #### Read FlowCyt_results_subset----
  
  setwd(Master_path_analysis)
  
  FlowCyt_results_subset<-readRDS(file="FlowCyt_global_adapted_plus_ilr.rds")
  
  
  
  cat("FlowCyt_results_subset_0\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FlowCyt_results_subset$treatment))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FlowCyt_results_subset$treatment)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FlowCyt_results_subset$Sample_label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FlowCyt_results_subset$Sample_label)))))
  cat("\n")
  
 
  
  ##### calculate principal components ------
  
  DEBUG <- 0
  
  
  Matrix_1<-as.matrix(FlowCyt_results_subset[,c(which(colnames(FlowCyt_results_subset) == 'CD41_MFI'),
                                                  which(colnames(FlowCyt_results_subset) == 'time'))])
  
  row.names(Matrix_1)<-FlowCyt_results_subset$Sample_label
  
  if(DEBUG ==1)
  {
    cat("Matrix_1_0\n")
    cat(str(Matrix_1))
    cat("\n")
  }
  
  
  PC_list<-princomp(Matrix_1)
  
  if(DEBUG ==1)
  {
    cat("PC_list_0\n")
    cat(str(PC_list))
    cat("\n")
  }
  
  PC_scores<-PC_list$scores
  
  if(DEBUG ==1)
  {
    cat("PC_scores_0\n")
    cat(str(PC_scores))
    cat("\n")
  }
  
  df_COMP<-as.data.frame(PC_scores, stringsAsFactors=F)
  colnames(df_COMP)<-paste(colnames(df_COMP),"CD41_MFI",sep="_")
  row.names(df_COMP)<-NULL
  df_COMP$Sample_label<-row.names(PC_scores)
  
  if(DEBUG ==1)
  {
    cat("df_COMP_0\n")
    cat(str(df_COMP))
    cat("\n")
  }
  
  FlowCyt_results_subset<-merge(FlowCyt_results_subset,
                                df_COMP,
                                by="Sample_label")
  
  cat("FlowCyt_results_subset_1\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
  
  Matrix_2<-as.matrix(FlowCyt_results_subset[,c(which(colnames(FlowCyt_results_subset) == 'CD41_GeoMFI'),
                                                which(colnames(FlowCyt_results_subset) == 'time'))])
  
  row.names(Matrix_2)<-FlowCyt_results_subset$Sample_label
  
  if(DEBUG ==1)
  {
    cat("Matrix_2_0\n")
    cat(str(Matrix_2))
    cat("\n")
  }
  
  
  PC_list<-princomp(Matrix_2)
  
  if(DEBUG ==1)
  {
    cat("PC_list_0\n")
    cat(str(PC_list))
    cat("\n")
  }
  
  PC_scores<-PC_list$scores
  
  if(DEBUG ==1)
  {
    cat("PC_scores_0\n")
    cat(str(PC_scores))
    cat("\n")
  }
  
  df_COMP<-as.data.frame(PC_scores, stringsAsFactors=F)
  colnames(df_COMP)<-paste(colnames(df_COMP),"CD41_GeoMFI",sep="_")
  row.names(df_COMP)<-NULL
  df_COMP$Sample_label<-row.names(PC_scores)
  
  if(DEBUG ==1)
  {
    cat("df_COMP_0\n")
    cat(str(df_COMP))
    cat("\n")
  }
  
  FlowCyt_results_subset<-merge(FlowCyt_results_subset,
                                df_COMP,
                                by="Sample_label")
  
  cat("FlowCyt_results_subset_2\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
  
  ################    SAVE   #######################
  
  setwd(out)
  
  
  write.table(FlowCyt_results_subset, file='FlowCyt_global_adapted_plus_ilr_and_princomp.csv', sep=",", quote=F, row.names = F)
  saveRDS(FlowCyt_results_subset, file='FlowCyt_global_adapted_plus_ilr_and_princomp.rds')
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--FlowCyt_results_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_path_analysis"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  ilr_transformation(opt)
  princomp_operation(opt)
  
  
  
}


###########################################################################

system.time( main() )