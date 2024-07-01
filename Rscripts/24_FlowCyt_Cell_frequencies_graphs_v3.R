
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

Graph_panoramic_scatter_plot = function(option_list)
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
  
  
  #### melt ----
  
  FlowCyt_results_subset.m<-melt(FlowCyt_results_subset, id.vars=c("Sample_label","sample","Genotype","time_point","treatment",'time'), variable.name='Feature', value.name='feature_magnitude')
  
  cat("FlowCyt_results_subset.m_0\n")
  cat(str(FlowCyt_results_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset.m$Feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset.m$Feature))))
  cat("\n")
  
  #### LOOP TO PRINT ----
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  path_graphs<-paste(out,'graphs','/','Feature_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  feature_array<-unique(FlowCyt_results_subset$Feature)
  
  cat("feature_array_0\n")
  cat(str(feature_array))
  cat("\n")
  
  
  test_index<-which(feature_array == 'Mean_FSC')
  
  
    
  feature_array<-c("CD41_MFI","CD41_GeoMFI")
  
  
  START<-test_index
  START<-1
  
  DEBUG<-0
  
  list_DEF<-list()
  
  for(i in seq(from=START, to=length(feature_array), by=1))
  {
    feature_array_sel<-feature_array[i]
    
   
    cat("----------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(feature_array_sel)))
    cat("\n")
      
   
    
    FlowCyt_results_subset.m_sel<-unique(droplevels(FlowCyt_results_subset.m[which(FlowCyt_results_subset.m$Feature %in% feature_array_sel),]))
    
    if(DEBUG == 1)
    {
      cat("FlowCyt_results_subset.m_sel_0\n")
      cat(str(FlowCyt_results_subset.m_sel))
      cat("\n")
    }
    
    REP<-FlowCyt_results_subset.m_sel
    
    ### feature_magnitude 1 ----
    
    
    ind.feature_1<-which(colnames(REP) == 'feature_magnitude')

    A<-round(summary(REP[,ind.feature_1][!is.na(REP[,ind.feature_1])]),2)

    step<-abs(A[6]-A[1])/10

    if(step == 0)
    {

      step<-1
    }

    if(DEBUG ==1)
    {
      cat("Summary_feature_magnitude\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }



    breaks.feature_1<-unique(sort(c(seq(from= A[1], to=A[6]+step,by=step))))
    labels.feature_1<-as.character(round(breaks.feature_1,1))


    if(DEBUG == 1)
    {
      cat("breaks.feature_1:\t")
      cat(sprintf(as.character(breaks.feature_1)))
      cat("\n")

      cat("labels.feature_1:\t")
      cat(sprintf(as.character(labels.feature_1)))
      cat("\n")

    }
    
    ### feature_magnitude 2 ----
    
    
    ind.feature_2<-which(colnames(REP) == 'time')
    
    # A<-round(summary(REP[,ind.feature_2][!is.na(REP[,ind.feature_2])]),2)
    # 
    # step<-abs(A[6]-A[1])/10
    # 
    # if(step == 0)
    # {
    #   
    #   step<-1
    # }
    # 
    # if(DEBUG ==1)
    # {
    #   cat("Summary_feature_magnitude\n")
    #   cat(sprintf(as.character(names(A))))
    #   cat("\n")
    #   cat(sprintf(as.character(A)))
    #   cat("\n")
    #   cat(sprintf(as.character(step)))
    #   cat("\n")
    # }
    # 
    
    
    breaks.feature_2<-c(0,16,24,48,72)
    labels.feature_2<-as.character(round(breaks.feature_2))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_2:\t")
      cat(sprintf(as.character(breaks.feature_2)))
      cat("\n")
      
      cat("labels.feature_2:\t")
      cat(sprintf(as.character(labels.feature_2)))
      cat("\n")
      
    }
    
    jitter_pos <- position_jitter(width=0.01, height = 0.01, seed = 1)


    graph_feature_magnitude<-ggplot() +
      geom_point(data=REP,
                 aes(x=REP[,ind.feature_2], 
                     y=REP[,ind.feature_1], 
                     color=sample,
                     shape=Genotype),
                 position=jitter_pos,size=4)+
      theme_bw()+
      scale_y_continuous(name=feature_array[i],breaks=breaks.feature_1,labels=labels.feature_1, limits=c(breaks.feature_1[1]-0.01,breaks.feature_1[length(breaks.feature_1)]+0.01))+
      scale_x_continuous(name='Time (hrs)',breaks=breaks.feature_2,labels=labels.feature_2, limits=c(breaks.feature_2[1]-0.01,breaks.feature_2[length(breaks.feature_2)]+0.01))+
      scale_color_brewer(palette = "Set3", drop=F)+
      scale_shape_manual(values=c(16,17,15))+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=14,color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "bottom")+
      guides(colour = guide_legend(nrow = 1))+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Part_I:\t")

    }
    

    setwd(path_graphs)

    svgname<-paste(paste(paste(feature_array_sel, collapse="__"),"Model_plot", sep='_'),".svg",sep='')
    makesvg = TRUE

    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph_feature_magnitude,
             device="svg",
             height=10, width=12)
    }

    if(DEBUG == 1)
    {
      cat("THE_END:------------------------------------------------------------------------------------------------\n")

    }
    
    
    # setwd(Master_path_analysis)
    # write.table(REP, file='test.tsv', sep="\t",quote=F, row.names = F)
    # 
    # quit(status = 1)
   
  }#i in 1:length(feature_array)
}

Graph_panoramic_dotplot_cell_frequencies = function(option_list)
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
  
  #### Read model_results----
  
  model_results<-as.data.frame(fread(file=opt$model_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("model_results_0\n")
  cat(str(model_results))
  cat("\n")
  
  #### melt ----
  
  FlowCyt_results_subset.m<-melt(FlowCyt_results_subset, id.vars=c("Sample_label","sample","Genotype","time_point","treatment",'time'), variable.name='Feature', value.name='feature_magnitude')
  
  cat("FlowCyt_results_subset.m_0\n")
  cat(str(FlowCyt_results_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset.m$Feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset.m$Feature))))
  cat("\n")
  
  #### LOOP TO PRINT ----
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  path_graphs<-paste(out,'graphs','/','Feature_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  feature_array<-unique(FlowCyt_results_subset$Feature)
  
  cat("feature_array_0\n")
  cat(str(feature_array))
  cat("\n")
  
  
  test_index<-which(feature_array == 'Mean_FSC')
  
  
  
  feature_array<-c("Double_neg","CD235_single","Double_pos","CD41_single")
  
  
  DEBUG<-0
  
  list_DEF<-list()
  

    
    
    FlowCyt_results_subset.m_sel<-unique(droplevels(FlowCyt_results_subset.m[which(FlowCyt_results_subset.m$Feature %in% feature_array),]))
    
    
    FlowCyt_results_subset.m_sel$Feature<-factor(FlowCyt_results_subset.m_sel$Feature,
                                                  levels=c("Double_neg","CD235_single","Double_pos","CD41_single"),
                                                  ordered=T)
    
    if(DEBUG == 1)
    {
      cat("FlowCyt_results_subset.m_sel_0\n")
      cat(str(FlowCyt_results_subset.m_sel))
      cat("\n")
    }
    
    REP_1<-FlowCyt_results_subset.m_sel
    
    if(DEBUG == 1)
    {
      cat("REP_1_0\n")
      cat(str(REP_1))
      cat("\n")
    }
    
    ### feature_magnitude 1 ----
    
    
    ind.feature_1<-which(colnames(REP_1) == 'feature_magnitude')
    
    # A<-round(summary(REP_1[,ind.feature_1][!is.na(REP_1[,ind.feature_1])]),2)
    # 
    # step<-abs(A[6]-A[1])/10
    # 
    # if(step == 0)
    # {
    #   
    #   step<-1
    # }
    # 
    # if(DEBUG ==1)
    # {
    #   cat("Summary_feature_magnitude\n")
    #   cat(sprintf(as.character(names(A))))
    #   cat("\n")
    #   cat(sprintf(as.character(A)))
    #   cat("\n")
    #   cat(sprintf(as.character(step)))
    #   cat("\n")
    # }
    
    
    
    breaks.feature_1<-seq(from= 0, to=100,by=10)
    labels.feature_1<-as.character(round(breaks.feature_1,0))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_1:\t")
      cat(sprintf(as.character(breaks.feature_1)))
      cat("\n")
      
      cat("labels.feature_1:\t")
      cat(sprintf(as.character(labels.feature_1)))
      cat("\n")
      
    }
    
    ### feature_magnitude 2 ----
    
    
    ind.feature_2<-which(colnames(REP_1) == 'time')
    
    # A<-round(summary(REP_1[,ind.feature_2][!is.na(REP_1[,ind.feature_2])]),2)
    # 
    # step<-abs(A[6]-A[1])/10
    # 
    # if(step == 0)
    # {
    #   
    #   step<-1
    # }
    # 
    # if(DEBUG ==1)
    # {
    #   cat("Summary_feature_magnitude\n")
    #   cat(sprintf(as.character(names(A))))
    #   cat("\n")
    #   cat(sprintf(as.character(A)))
    #   cat("\n")
    #   cat(sprintf(as.character(step)))
    #   cat("\n")
    # }
    # 
    
    
    breaks.feature_2<-c(0,16,24,48,72)
    labels.feature_2<-as.character(round(breaks.feature_2))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_2:\t")
      cat(sprintf(as.character(breaks.feature_2)))
      cat("\n")
      
      cat("labels.feature_2:\t")
      cat(sprintf(as.character(labels.feature_2)))
      cat("\n")
      
    }
    
    jitter_pos <- position_jitter(width=0, height = 0.2, seed = 1)
    
    
    graph_1<-ggplot() +
      geom_point(data=REP_1,
                 aes(x=REP_1[,ind.feature_2], 
                     y=REP_1[,ind.feature_1],
                     color=sample,
                     shape=Genotype),
                 position=jitter_pos,size=4)+
      theme_bw()+
      scale_y_continuous(name='Cell_frequencies',breaks=breaks.feature_1,labels=labels.feature_1, limits=c(breaks.feature_1[1]-0.02,breaks.feature_1[length(breaks.feature_1)]+0.02))+
      scale_x_continuous(name='Time (hrs)',breaks=breaks.feature_2,labels=labels.feature_2, limits=c(breaks.feature_2[1]-0.2,breaks.feature_2[length(breaks.feature_2)]+0.2))+
      scale_color_brewer(palette = "Set3", drop=F)+
      scale_shape_manual(values=c(16,17,15))
    
    graph_1<-graph_1+
      facet_grid(. ~ REP_1$Feature, scales='free_x', space='free_x') +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_blank(),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "bottom")+
      guides(colour = guide_legend(nrow = 1))+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("End_graph_1:\n----------------------------------------------------")
      
    }
    
    
    
    ####################### _ilr_value -------
    
    feature_array<-c("V1_ilr_value","V2_ilr_value","V3_ilr_value")
    
    
    
    FlowCyt_results_subset.m_sel<-unique(droplevels(FlowCyt_results_subset.m[which(FlowCyt_results_subset.m$Feature %in% feature_array),]))
    
    FlowCyt_results_subset.m_sel$Feature<-factor(FlowCyt_results_subset.m_sel$Feature,
                                                 levels=c("V1_ilr_value","V2_ilr_value","V3_ilr_value"),
                                                 ordered=T)
    
    if(DEBUG == 1)
    {
      cat("FlowCyt_results_subset.m_sel_0\n")
      cat(str(FlowCyt_results_subset.m_sel))
      cat("\n")
    }
    
    REP_2<-FlowCyt_results_subset.m_sel
    
    if(DEBUG == 1)
    {
      cat("REP_2_0\n")
      cat(str(REP_2))
      cat("\n")
    }
    
    ### feature_magnitude 1 ----
    
    
    ind.feature_1<-which(colnames(REP_2) == 'feature_magnitude')
    
    A<-round(summary(REP_2[,ind.feature_1][!is.na(REP_2[,ind.feature_1])]),2)

    step<-abs(A[6]-A[1])/10

    if(step == 0)
    {

      step<-1
    }

    if(DEBUG ==1)
    {
      cat("Summary_feature_magnitude\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }
    
    
    
    breaks.feature_1<-unique(sort(c(seq(from= A[1], to=A[6]+step,by=step))))
    labels.feature_1<-as.character(round(breaks.feature_1,1))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_1:\t")
      cat(sprintf(as.character(breaks.feature_1)))
      cat("\n")
      
      cat("labels.feature_1:\t")
      cat(sprintf(as.character(labels.feature_1)))
      cat("\n")
      
    }
    
    ### feature_magnitude 2 ----
    
    
    ind.feature_2<-which(colnames(REP_2) == 'time')
    
    # A<-round(summary(REP_2[,ind.feature_2][!is.na(REP_2[,ind.feature_2])]),2)
    # 
    # step<-abs(A[6]-A[1])/10
    # 
    # if(step == 0)
    # {
    #   
    #   step<-1
    # }
    # 
    # if(DEBUG ==1)
    # {
    #   cat("Summary_feature_magnitude\n")
    #   cat(sprintf(as.character(names(A))))
    #   cat("\n")
    #   cat(sprintf(as.character(A)))
    #   cat("\n")
    #   cat(sprintf(as.character(step)))
    #   cat("\n")
    # }
    # 
    
    
    breaks.feature_2<-c(0,16,24,48,72)
    labels.feature_2<-as.character(round(breaks.feature_2))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_2:\t")
      cat(sprintf(as.character(breaks.feature_2)))
      cat("\n")
      
      cat("labels.feature_2:\t")
      cat(sprintf(as.character(labels.feature_2)))
      cat("\n")
      
    }
    
    jitter_pos <- position_jitter(width=0, height = 0.2, seed = 1)
    
    
    graph_2<-ggplot() +
      geom_point(data=REP_2,
                 aes(x=REP_2[,ind.feature_2], 
                     y=REP_2[,ind.feature_1],
                     color=sample,
                     shape=Genotype),
                 position=jitter_pos,size=4)+
      theme_bw()+
      scale_y_continuous(name='Cell_frequencies ilr values',breaks=breaks.feature_1,labels=labels.feature_1, limits=c(breaks.feature_1[1]-0.01,breaks.feature_1[length(breaks.feature_1)]+0.01))+
      scale_x_continuous(name='Time (hrs)',breaks=breaks.feature_2,labels=labels.feature_2, limits=c(breaks.feature_2[1]-0.2,breaks.feature_2[length(breaks.feature_2)]+0.2))+
      scale_color_brewer(palette = "Set3", drop=F)+
      scale_shape_manual(values=c(16,17,15))
    
    graph_2<-graph_2+
      facet_grid(. ~ REP_2$Feature, scales='free_x', space='free_x') +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=14,color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "bottom")+
      guides(colour = guide_legend(nrow = 1))+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("End_graph_2:\n----------------------------------------------------")
      
    }
    
    
    graph_DEF<-plot_grid(graph_1,graph_2,
                         ncol = 1,
                         nrow=2,
                         rel_heights=c(1,1))
    
    
    
    
    setwd(path_graphs)
    
    svgname<-paste(paste("Panoramic_plot","Cell_frequencies", sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph_DEF,
             device="svg",
             height=10, width=12)
    }
    
    
    # setwd(out)
    # write.table(file="test.tsv",REP_2,quote=F,row.names = F)
    
    
    if(DEBUG == 1)
    {
      cat("THE_END:------------------------------------------------------------------------------------------------\n")
      
    }

    

}

Graph_selected_scatter_plot = function(option_list)
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
  
  
  
  #### melt ----
  
  FlowCyt_results_subset.m<-melt(FlowCyt_results_subset, id.vars=c("Sample_label","sample","Genotype","time_point","treatment",'time'), variable.name='Feature', value.name='feature_magnitude')
  
  cat("FlowCyt_results_subset.m_0\n")
  cat(str(FlowCyt_results_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset.m$Feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset.m$Feature))))
  cat("\n")
  
  #### LOOP TO PRINT ----
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  path_graphs<-paste(out,'graphs','/','Feature_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  feature_array<-unique(FlowCyt_results_subset$Feature)
  
  cat("feature_array_0\n")
  cat(str(feature_array))
  cat("\n")
  
  
  test_index<-which(feature_array == 'Mean_FSC')
  
  
  
  feature_array<-c("CD41_MFI","CD41_GeoMFI")
  
  
  START<-test_index
  START<-1
  
  DEBUG<-0
  
  list_DEF<-list()
  
  for(i in seq(from=START, to=length(feature_array), by=1))
  {
    feature_array_sel<-feature_array[i]
    
    
    cat("----------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(feature_array_sel)))
    cat("\n")
    
    
    
    FlowCyt_results_subset.m_sel<-unique(droplevels(FlowCyt_results_subset.m[which(FlowCyt_results_subset.m$Feature %in% feature_array_sel),]))
    
    if(DEBUG == 1)
    {
      cat("FlowCyt_results_subset.m_sel_0\n")
      cat(str(FlowCyt_results_subset.m_sel))
      cat("\n")
    }
    
    REP<-FlowCyt_results_subset.m_sel
    
    ### feature_magnitude 1 ----
    
    
    ind.feature_1<-which(colnames(REP) == 'feature_magnitude')
    
    A<-round(summary(REP[,ind.feature_1][!is.na(REP[,ind.feature_1])]),2)
    
    step<-abs(A[6]-A[1])/10
    
    if(step == 0)
    {
      
      step<-1
    }
    
    if(DEBUG ==1)
    {
      cat("Summary_feature_magnitude\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }
    
    
    
    breaks.feature_1<-unique(sort(c(seq(from= 0, to=45000,by=5000))))
    labels.feature_1<-as.character(round(breaks.feature_1))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_1:\t")
      cat(sprintf(as.character(breaks.feature_1)))
      cat("\n")
      
      cat("labels.feature_1:\t")
      cat(sprintf(as.character(labels.feature_1)))
      cat("\n")
      
    }
    
    ### feature_magnitude 2 ----
    
    
    ind.feature_2<-which(colnames(REP) == 'time')
    
    # A<-round(summary(REP[,ind.feature_2][!is.na(REP[,ind.feature_2])]),2)
    # 
    # step<-abs(A[6]-A[1])/10
    # 
    # if(step == 0)
    # {
    #   
    #   step<-1
    # }
    # 
    # if(DEBUG ==1)
    # {
    #   cat("Summary_feature_magnitude\n")
    #   cat(sprintf(as.character(names(A))))
    #   cat("\n")
    #   cat(sprintf(as.character(A)))
    #   cat("\n")
    #   cat(sprintf(as.character(step)))
    #   cat("\n")
    # }
    # 
    
    
    breaks.feature_2<-c(0,16,24,48,72)
    labels.feature_2<-as.character(round(breaks.feature_2))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.feature_2:\t")
      cat(sprintf(as.character(breaks.feature_2)))
      cat("\n")
      
      cat("labels.feature_2:\t")
      cat(sprintf(as.character(labels.feature_2)))
      cat("\n")
      
    }
    
    jitter_pos <- position_jitter(width=0.01, height = 0.01, seed = 1)
    
    
    graph_feature_magnitude<-ggplot() +
      geom_point(data=REP,
                 aes(x=REP[,ind.feature_2], 
                     y=REP[,ind.feature_1], 
                     shape=Genotype),
                 position=jitter_pos,size=5, color="black")+
      theme_bw()+
      scale_y_continuous(name=feature_array[i],breaks=breaks.feature_1,labels=labels.feature_1, limits=c(breaks.feature_1[1]-0.02,breaks.feature_1[length(breaks.feature_1)]+0.02))+
      scale_x_continuous(name='Time (hrs)',breaks=breaks.feature_2,labels=labels.feature_2, limits=c(breaks.feature_2[1],breaks.feature_2[length(breaks.feature_2)]+0.02))+
      scale_color_brewer(palette = "Set3", drop=F)+
      scale_shape_manual(values=c(16,17,15))+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=14,color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "bottom")+
      guides(colour = guide_legend(nrow = 1))+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Part_I:\t")
      
    }
    
    
    setwd(path_graphs)
    
    svgname<-paste(paste(paste(feature_array_sel, collapse="__"),"Selected_plot", sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph_feature_magnitude,
             device="svg",
             height=10, width=12)
    }
    
    if(DEBUG == 1)
    {
      cat("THE_END:------------------------------------------------------------------------------------------------\n")
      
    }
    
    
    # setwd(Master_path_analysis)
    # write.table(REP, file='test.tsv', sep="\t",quote=F, row.names = F)
    # 
    # quit(status = 1)
    
  }#i in 1:length(feature_array)
}

Graph_selected_cell_frequencies = function(option_list)
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
  
  #### Read model_results----
  
  model_results<-as.data.frame(fread(file=opt$model_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("model_results_0\n")
  cat(str(model_results))
  cat("\n")
  
  
  #### Read FlowCyt_results_subset----
  
  FlowCyt_results_subset<-readRDS(file=opt$FlowCyt_results_subset)
  
  cat("FlowCyt_results_subset_0\n")
  cat(str(FlowCyt_results_subset))
  cat("\n")
  
  
  #### melt ----
  
  FlowCyt_results_subset.m<-melt(FlowCyt_results_subset, id.vars=c("Sample_label","sample","Genotype","time_point","treatment",'time'), variable.name='Feature', value.name='feature_magnitude')
  
  cat("FlowCyt_results_subset.m_0\n")
  cat(str(FlowCyt_results_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(FlowCyt_results_subset.m$Feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(FlowCyt_results_subset.m$Feature))))
  cat("\n")
  
  #### LOOP TO PRINT ----
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  path_graphs<-paste(out,'graphs','/','Feature_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  feature_array<-unique(FlowCyt_results_subset$Feature)
  
  cat("feature_array_0\n")
  cat(str(feature_array))
  cat("\n")
  
  
  test_index<-which(feature_array == 'Mean_FSC')
  
  
  
  feature_array<-c("Double_neg","CD235_single","Double_pos","CD41_single")
  
  
  DEBUG<-1
  
  list_DEF<-list()
  
  
  
  
  FlowCyt_results_subset.m_sel<-unique(droplevels(FlowCyt_results_subset.m[which(FlowCyt_results_subset.m$Feature %in% feature_array),]))
  
  
  FlowCyt_results_subset.m_sel$Feature<-factor(FlowCyt_results_subset.m_sel$Feature,
                                               levels=c("Double_neg","CD235_single","Double_pos","CD41_single"),
                                               ordered=T)
  
  if(DEBUG == 1)
  {
    cat("FlowCyt_results_subset.m_sel_0\n")
    cat(str(FlowCyt_results_subset.m_sel))
    cat("\n")
  }
  
 
  FlowCyt_results_subset.m_sel.dt<-data.table(FlowCyt_results_subset.m_sel, key=c("Genotype","time_point","treatment",'time','Feature'))
  
  
  
  Freq.table_mean<-as.data.frame(FlowCyt_results_subset.m_sel.dt[,.(mean_feature_magnitude =mean(feature_magnitude)), by=key(FlowCyt_results_subset.m_sel.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Freq.table_mean_0\n")
    cat(str(Freq.table_mean))
    cat("\n")
  }
  
  
  # quit(status = 1)
  
  #### Graph
  
  breaks.Rank<-(seq(0,100,by=25))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors<-brewer.pal(8,"Reds")[c(1,3,5,7)]
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
  
  
  stacked_barplot<-ggplot(data=Freq.table_mean,
                              aes(x=Genotype, y=mean_feature_magnitude, fill=Feature)) +
    geom_bar(stat="identity",colour='black')+
    scale_y_continuous(name=paste("Percentage of cells",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors,drop=F)
  
 
  

  stacked_barplot<-stacked_barplot+
    facet_grid(. ~ time_point, scales='free_x', space='free_x', switch="y", drop=F)+
    theme_cowplot(font_size = 4)+
    theme( strip.background = element_blank(),
           strip.placement = "outside",
           strip.text = element_text(size=6),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="white",size=0,5),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme_classic()+
    theme(axis.title.y=element_text(size=8, color="black", family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=6, color="black", family="sans"),
          axis.text.x=element_text(angle=45, hjust=1,size=6, color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4))+
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.25, 'cm'), #change legend key size
          legend.key.height = unit(0.25, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.position="right")+
    ggeasy::easy_center_title()
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Stacked_barplot_mean","Cell_frequencies", sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 5, height = 3)
  print(stacked_barplot)
  dev.off()
  
  ### graph parameters_Minuslogpval ----
  
  model_results_subset<-unique(model_results[,c(which(colnames(model_results)=='comparison'),
                                         which(colnames(model_results)=='Minuslogpval'))])
  
  model_results_subset$COORD<-"1"
  
  model_results_subset$comparison<-factor(model_results_subset$comparison,
                                     levels=rev(c('wt_vs_homALT','wt_vs_Del80','homALT_vs_Del80')),
                                     ordered=T)
  
  cat("model_results_subset_0\n")
  cat(str(model_results_subset))
  cat("\n")
  
  
  indx_Minuslogpval<-which(colnames(model_results_subset) == 'Minuslogpval')
  
  A_Minuslogpval<-summary(model_results_subset[,indx_Minuslogpval])
  
  
  if(DEBUG == 1)
  {
    cat("A_Minuslogpval\n")
    cat(sprintf(as.character(names(A_Minuslogpval))))
    cat("\n")
    cat(sprintf(as.character(A_Minuslogpval)))
    cat("\n")
  }
  
  max_value<-A_Minuslogpval[6]
  min_value<-A_Minuslogpval[1]
  
  
  step<-round(abs(max_value-min_value)/3,1)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks_Minuslogpval<-unique(sort(unique(c(1,seq(min_value,max_value+step, by=step)))))
  labels_Minuslogpval<-as.character(round(breaks_Minuslogpval),1)
  
  if(DEBUG == 1)
  {
    cat("step_Minuslogpval\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("labels_Minuslogpval\n")
    cat(sprintf(as.character(labels_Minuslogpval)))
    cat("\n")
  }
  
  vector_colors<-brewer.pal(3, "Dark2")[c(3,2,1)]
  
  if(DEBUG == 1)
  {
    cat("vector_colors\n")
    cat(str(vector_colors))
    cat("\n")
  }
  
  #### dotplot
  
  if(DEBUG == 1)
  {
    cat("Gene_set_dotplot_START:\n")
    
  }
  
  
  Model_dotplot<-ggplot(data=model_results_subset,
                           aes(y=as.numeric(comparison),
                               x=COORD))+
    geom_point(aes(size=Minuslogpval, fill=comparison), stroke=1, shape=21)+
    scale_size(range = c(0,10), name='-log10pval',
               breaks=breaks_Minuslogpval, labels=labels_Minuslogpval, limits=c(breaks_Minuslogpval[1]-1,breaks_Minuslogpval[length(breaks_Minuslogpval)]+1))+
    scale_y_discrete(name=NULL, breaks=NULL)+
    scale_x_discrete(name=NULL, breaks=NULL)+
    scale_fill_manual(values=vector_colors, drop=F)+
    theme_classic()+
    theme(axis.title.y=element_text(size=8, color="black", family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=6, color="black", family="sans"),
          axis.text.x=element_text(angle=45, vjust=1,size=6, color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4))+
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.25, 'cm'), #change legend key size
          legend.key.height = unit(0.25, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.position="right")+
    ggeasy::easy_center_title()
  
  if(DEBUG == 1)
  {
    cat("Gene_set_dotplot_END:\n")
  }
  
  ##### plot model coefficients ----
  
  
  model_results.m<-melt(model_results, id.vars=c('term.id',"pval","comparison","Minuslogpval"), variable.name='Cell_type', value.name='Coefficient')
  
  cat("model_results.m_0\n")
  cat(str(model_results.m))
  cat("\n")
  
  model_results.m$Cell_type<-factor(model_results.m$Cell_type,
                                    levels=c("Double_neg","CD235_single","Double_pos","CD41_single"),
                                    ordered=T)
  
  cat("model_results.m_1\n")
  cat(str(model_results.m))
  cat("\n")
  cat(sprintf(as.character(unique(model_results.m$comparison))))
  cat("\n")
  
  levels_comp<-levels(model_results.m$comparison)
  
  model_results.m$comparison<-factor(model_results.m$comparison,
                                     levels=c('wt_vs_homALT','wt_vs_Del80','homALT_vs_Del80'),
                                     ordered=T)
  
  
  
  cat("model_results.m_2\n")
  cat(str(model_results.m))
  cat("\n")
  
 

  A_Coefficient<-summary(model_results.m$Coefficient)
  
  
  if(DEBUG == 1)
  {
    cat("A_Coefficient\n")
    cat(sprintf(as.character(names(A_Coefficient))))
    cat("\n")
    cat(sprintf(as.character(A_Coefficient)))
    cat("\n")
  }
  
  max_value<-A_Coefficient[6]
  min_value<-A_Coefficient[1]
  
  
  step<-round(abs(max_value-min_value)/3,1)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks_Coefficient<-unique(sort(unique(c(seq(min_value,max_value+step, by=step)))))
  labels_Coefficient<-as.character(breaks_Coefficient)
  
  breaks_Coefficient<-seq(from=0, to=1, by=0.2)
  labels_Coefficient<-as.character(breaks_Coefficient)
  
  if(DEBUG == 1)
  {
    cat("step_Coefficient\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("breaks_Coefficient\n")
    cat(sprintf(as.character(breaks_Coefficient)))
    cat("\n")
    cat("labels_Coefficient\n")
    cat(sprintf(as.character(labels_Coefficient)))
    cat("\n")
  }
  
  
  vector_colors<-brewer.pal(8,"Reds")[c(1,3,5,7)]
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
  
  

  
  
  
 
  
  
  stacked_barplot<-ggplot(data=model_results.m,
                          aes(x=term.id, y=Coefficient, fill=Cell_type)) +
    geom_bar(stat="identity",colour='black')+
    scale_y_continuous(name=paste("Coefficient in model",sep=" "),breaks=breaks_Coefficient,labels=labels_Coefficient,
                       limits=c(breaks_Coefficient[1]-0.001,breaks_Coefficient[length(breaks_Coefficient)]+0.001))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors,drop=F)
  
  
  
  stacked_barplot<-stacked_barplot +
    facet_grid(. ~ comparison, scales='free_x', space='free_x', switch="y", drop=F)+
    theme_cowplot(font_size = 4)+
    theme( strip.background = element_blank(),
           strip.placement = "outside",
           strip.text = element_text(size=6),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="white",size=0,5),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme_classic()+
    theme(axis.title.y=element_text(size=8, color="black", family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=6, color="black", family="sans"),
          axis.text.x=element_text(angle=45, hjust=1,size=6, color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4))+
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.25, 'cm'), #change legend key size
          legend.key.height = unit(0.25, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.position="right")+
    ggeasy::easy_center_title()
  
  setwd(path_graphs)
  
  graph_DEF<-plot_grid(Model_dotplot,stacked_barplot,
                       ncol = 2,
                       nrow=1,
                       rel_widths=c(0.5,1))
  
  svgname<-paste(paste("Model_plots","Cell_frequencies", sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 6, height = 4)
  print(graph_DEF)
  dev.off()
  
  setwd(out)
  write.table(file='test.tsv',model_results_subset, sep="\t", quote=F, row.names = F)
  write.table(file='test2.tsv',model_results.m, sep="\t", quote=F, row.names = F)
  
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
    make_option(c("--model_results"), type="character", default=NULL, 
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
  
 
  Graph_panoramic_scatter_plot(opt)
  Graph_panoramic_dotplot_cell_frequencies(opt)
  Graph_selected_scatter_plot(opt)
  Graph_selected_cell_frequencies(opt)


  
}


###########################################################################

system.time( main() )