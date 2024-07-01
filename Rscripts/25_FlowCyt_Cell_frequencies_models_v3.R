
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
suppressMessages(library("tidyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)


model_princomp = function(option_list)
{
  
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
  
  #### LOOP TO PRINT ----
  
  path_models<-paste(out,'models','/',sep='')
  
  if (file.exists(path_models)){
    
  }else{
    
    dir.create(file.path(path_models))
    
  }#path_models
  
  path_models<-paste(out,'models','/','princomp_model','/',sep='')
  
  if (file.exists(path_models)){
    
  }else{
    
    dir.create(file.path(path_models))
    
  }#path_models
  
  
  
  #### LOOP to select specific genotype comparisons----
  
  array_comparisons<-c("wt","homALT","wt","Del80","homALT","Del80")
  
  DEBUG <-0
  
  START<-1
  
  result_DEF<-data.frame()
  
  for(i in seq(from=START, to=length(array_comparisons), by=2))
  {
    array_comparisons_sel<-c(array_comparisons[i],array_comparisons[i+1])
    
    
    cat("----------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(array_comparisons_sel)))
    cat("\n")
    
    FlowCyt_results_subset_G_sel<-droplevels(FlowCyt_results_subset[which(FlowCyt_results_subset$Genotype%in%array_comparisons_sel),])
    
    if(DEBUG ==1)
    {
      cat("FlowCyt_results_subset_G_sel_0\n")
      cat(str(FlowCyt_results_subset_G_sel))
      cat("\n")
    }
    
    #### linear model of Comp.1 against genotype and cell type ----
    
    model_1<-lm(Comp.1_CD41_MFI ~ time + Genotype + Genotype:time, data=FlowCyt_results_subset_G_sel)
    
    # if(DEBUG ==1)
    # {
    #   cat("model_1_0\n")
    #   cat(str(model_1))
    #   cat("\n")
    # }
    
    summary_model_1<-summary(model_1)
    
    # if(DEBUG ==1)
    # {
    #   cat("summary_model_1_0\n")
    #   cat(str(summary_model_1))
    #   cat("\n")
    # }
    
    coefficient_matrix<-as.data.frame(summary_model_1$coefficients)
    coefficient_matrix$term.id<-row.names(coefficient_matrix)
    
    if(DEBUG ==1)
    {
      cat("coefficient_matrix_0\n")
      cat(str(coefficient_matrix))
      cat("\n")
    }
    
    
    coefficient_matrix_subset<-coefficient_matrix[,c(which(colnames(coefficient_matrix)=='term.id'),
                                                     4)]
    coefficient_matrix_subset$feature<-paste('Comp.1_CD41_MFI')
    coefficient_matrix_subset$comparison<-paste(array_comparisons_sel, collapse = "_vs_")
    
    result_DEF<-rbind(result_DEF,coefficient_matrix_subset)
    
  }# i in seq(from=START, to=length(array_comparisons), by=2
  
  
  if(dim(result_DEF)[1] >0)
  {
    if(DEBUG ==1)
    {
      cat("result_DEF_0\n")
      cat(str(result_DEF))
      cat("\n")
    }
    
    setwd(path_models)
    
    write.table(result_DEF, file="Princomp_model_result.tsv", sep="\t", quote=F, row.names = F)
    
   
    
  }# dim(result_DEF)[1] >0
  
  
  
  
  
  
  
  
 
}

model_compositional = function(option_list)
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
  
  #### melt ----
  
  DEBUG <-1
  
  
  
  #### LOOP TO PRINT ----
  
  path_models<-paste(out,'models','/',sep='')
  
  if (file.exists(path_models)){
    
  }else{
    
    dir.create(file.path(path_models))
    
  }#path_models
  
  path_models<-paste(out,'models','/','compositional','/',sep='')
  
  if (file.exists(path_models)){
    
  }else{
    
    dir.create(file.path(path_models))
    
  }#path_models
  
  
  
  #### LOOP to select specific genotype comparisons----
  
  array_comparisons<-c("wt","homALT","wt","Del80","homALT","Del80")
  
  
  
  START<-1
  
  result_DEF<-data.frame()
  
  for(i in seq(from=START, to=length(array_comparisons), by=2))
  {
    array_comparisons_sel<-c(array_comparisons[i],array_comparisons[i+1])
    
    
    cat("-------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(array_comparisons_sel)))
    cat("\n")
    
    FlowCyt_results_subset_G_sel<-droplevels(FlowCyt_results_subset[which(FlowCyt_results_subset$Genotype%in%array_comparisons_sel),])
    
    if(DEBUG ==1)
    {
      cat("FlowCyt_results_subset_G_sel_0\n")
      cat(str(FlowCyt_results_subset_G_sel))
      cat("\n")
    }
    
    model_df<-FlowCyt_results_subset_G_sel
    
    if(DEBUG ==1)
    {
      cat("model_df_0\n")
      cat(str(model_df))
      cat("\n")
    }
    
    
    reduced_lm<-lm(cbind(V1_ilr_value,V2_ilr_value,V3_ilr_value) ~ time, data=model_df)
    intermediate_lm<-lm(cbind(V1_ilr_value,V2_ilr_value,V3_ilr_value) ~ time +Genotype, data=model_df)
    
    
    model_comparison<-as.data.frame(anova(reduced_lm,intermediate_lm), stringsAsFactors=F)
    
    if(DEBUG ==1)
    {
      cat("model_comparison_0\n")
      cat(str(model_comparison))
      cat("\n")
    }
    
    pval<-model_comparison[,which(colnames(model_comparison) == 'Pr(>F)')][!is.na(model_comparison[,which(colnames(model_comparison) == 'Pr(>F)')])]
    
    
    
    if(DEBUG ==1)
    {
      cat("pval_0\n")
      cat(str(pval))
      cat("\n")
    }
    
    
    indx.int<-c(which(colnames(model_df) == 'Double_neg'),which(colnames(model_df) == 'CD235_single'),which(colnames(model_df) == 'Double_pos'),which(colnames(model_df) == 'CD41_single'))
    
    Repeat_matrix<-as.matrix(model_df[,indx.int])
    row.names(Repeat_matrix)<-model_df$Sample_label
    
    if(DEBUG ==1)
    {
      cat("Repeat_matrix_0\n")
      cat(str(Repeat_matrix))
      cat("\n")
    }
    
    coefs<-as.data.frame(ilrInv(coef(intermediate_lm),orig=Repeat_matrix), stringAsFactors=F)
    
    colnames(coefs)<-colnames(Repeat_matrix)
    coefs$term.id<-row.names(coefs)
    
    
    if(DEBUG ==1)
    {
      cat("coefs_0\n")
      cat(str(coefs))
      cat("\n")
    }
    
    coefs$pval<-as.numeric(pval)
    coefs$comparison<-paste(array_comparisons_sel, collapse = "_vs_")
    coefs$Minuslogpval<--1*log10(coefs$pval)
   
    if(DEBUG ==1)
    {
      cat("coefs_1\n")
      cat(str(coefs))
      cat("\n")
    }
    
    result_DEF<-rbind(result_DEF,coefs)
    
  }# i in seq(from=START, to=length(array_comparisons), by=2
  
  
  if(dim(result_DEF)[1] >0)
  {
    if(DEBUG ==1)
    {
      cat("result_DEF_0\n")
      cat(str(result_DEF))
      cat("\n")
    }
    
    setwd(path_models)
    
    write.table(result_DEF, file="compositional_ilr_model_result.tsv", sep="\t", quote=F, row.names = F)
    
    
    
  }# dim(result_DEF)[1] >0
  
  
  
  
  
  
  
  
  
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
  
 
  model_princomp(opt)
  model_compositional(opt)


  
}


###########################################################################

system.time( main() )