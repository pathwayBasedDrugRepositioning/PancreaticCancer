###Running this code generates Fig2 of the paper and saves it in the figs directory


#Generates individual sub panels of figure 2 by taking parts of the original data, visualizing them and saving figures in the figs directory
MainPlot<-function(){
  library("feather")
  library("ggplot2")
  library("ComplexHeatmap")
  library("stringr")
  library("circlize")  
  
  path_ACP=read.csv("./data/ACP Pathways.txt",sep="\t",header = TRUE)
  path_CP=read.csv("./data/CP pathways.txt",sep="\t",header = TRUE)
  
  ACP_IDs=path_ACP[,"ID"];
  path_ACP_names=path_ACP[,"Pathway.Name"];
  names(path_ACP_names)=ACP_IDs;
  
  CP_ID=path_CP[,"ID"];
  path_IDs=c(ACP_IDs,CP_ID);
  path_IDs=unique(path_IDs);
  rr=which(is.na(path_IDs))
  if(length(rr)!=0){ path_IDs=path_IDs[-rr];}
  path_IDs=as.character(path_IDs);
  path_IDs=c("InstanceID",path_IDs);
  
  Matrix=read_feather("./data/LINCS_pathway_subset.feather",columns = path_IDs)
  MatIDs=Matrix[,1];
  MatIDs=as.character(MatIDs[[1]])
  
  disease=read.csv("./data/disease.txt",sep="\t",header = TRUE)
  compounds=read.csv("./data/Compounds.txt",sep="\t",header = TRUE)
  
  instanceIDs=compounds[2]
  instanceIDs=as.character(instanceIDs[[1]]);
  compoundNames=compounds[1][[1]]
  
  
  ii=which(MatIDs %in%  instanceIDs)
  
  Matrix=Matrix[ii,]
  MatIDs=Matrix[,1];
  MatIDs=as.character(MatIDs[[1]])
  
   draw_disease(Matrix,disease,path_CP)
  
  ht_acp=draw_acp_plot(compounds,path_ACP,Matrix,MatIDs,disease)
  ht_Score1=draw_ht_score(ScoreType="Score1",bshow_legend=TRUE,compounds,path_CP,Matrix,MatIDs,disease)
  ht_Score2=draw_ht_score(ScoreType="Score2",bshow_legend=FALSE,compounds,path_CP,Matrix,MatIDs,disease)
  ht_Res_score=draw_ht_score(ScoreType="Res-score",bshow_legend=FALSE,compounds,path_CP,Matrix,MatIDs,disease)
  ht_Res_score=draw_ht_score(ScoreType="Selected",bshow_legend=FALSE,compounds,path_CP,Matrix,MatIDs,disease)
  
}

#################################################

##Draws disease sub-panel figure and saves it
draw_disease<-function(Matrix,disease,path_CP){
  ii=which(path_CP[,"Score.Type"] %in% c("Score1","Score2"))
  path_CP_s=path_CP[ii,];
  pID=path_CP_s[,"ID"]
  pName=path_CP_s[,"Pathway.Name"];
  pName=str_sub(pName,1,70)
  names(pName)=pID;
  
  ii=which(disease[,"ID"] %in% pID)
  
  dis2=disease[ii,1:5];
  rownames(dis2)=pName[as.character(disease[ii,"ID"])]
  
  col1_fun = colorRamp2(c(-8,-5, 0,5, 8), c(rgb(0,32,96,maxColorValue =255),rgb(0,176,240,maxColorValue =255),
                                            "white",rgb(255,72,72,maxColorValue =255), rgb(192,0,0,maxColorValue =255))) 
  w=ncol(dis2) * 0.313 *1.313
  h= nrow(dis2) * 0.313 *1.313

  ### limit rowName length
 nameLength=26
 rnames2=str_sub(rownames(dis2),1,nameLength) 
 ii=which(nchar(rownames(dis2))>nameLength)
 added=rep("",nrow(dis2))
 added[ii]="..."
 rnames2=paste0(rnames2,added)
 rownames(dis2)=rnames2
 ##

  ht=Heatmap(dis2,name = "Enrichment"
             
             ,col = col1_fun
             #,heatmap_legend_param =list(legend_direction = "vertical" ,fontsize=12,fontface="bold")
             ,width = unit(w,"cm")
             ,height = unit(h,"cm") 
             ,row_names_gp = gpar(fontsize = 8,fontface="bold")
             ,column_names_gp = gpar(fontsize = 8,fontface="bold")
             
             ,cluster_rows = TRUE
             ,cluster_columns = FALSE
             ,column_title = "Pancreatic cancer"
             ,column_title_gp = gpar(fontsize=8,fontface="bold")
             ,show_heatmap_legend = FALSE
             ,column_names_rot=45
             ,row_names_rot=45
             
  ) 
  
  filename=paste0("./figs/disease.tiff")
  if(length(dev.list())!=0)dev.off()
  
  #tiff(filename = filename,width = (w+9),height = (h+8),units = "cm",res = 600)
  tiff(filename = filename,width = 7.5,height = 8.5,units = "cm",res = 600)
  
  draw(ht,heatmap_legend_side ="left",padding=unit(c(0,0,2,0),"cm"))
  dev.off()
  return(ht);
  
}

#################################################

#draws ACP plot and saves it
draw_acp_plot<-function(compounds,path_ACP,Matrix,MatIDs,disease){

  ii=which(abs(path_ACP[,"Panc1"])>5 & abs(path_ACP[,"Gemcitabine1"])>5)
  path_ACP2=path_ACP[ii,];
  rem=c("105801","730317","105818","160942","105807","105765")
  ii=which( path_ACP2[,"ID"] %in% rem)
  path_ACP2=path_ACP2[-ii,]
  
  path_ACP_ID=path_ACP2[,"ID"];
  path_ACP_pathName=path_ACP2[,"Pathway.Name"];
  names(path_ACP_pathName)=path_ACP_ID;
  
  ii=which(duplicated(compounds[,2]))
  compounds=compounds[-ii,]
  
  compound_IDs = compounds[,2];
  compound_names = compounds[,1];
  names(compound_names)=compound_IDs;
  
  Matrix_ACP=Matrix[,as.character(path_ACP_ID)]
  colnames(Matrix_ACP)= path_ACP_pathName[colnames(Matrix_ACP)];
  
  names=compound_names[as.character(Matrix[,1][[1]])]
  rownames(Matrix_ACP)=names;
  #########
  ii=which(disease[,"ID"] %in% path_ACP_ID)
  Panc1 = as.vector(unlist(disease[ii,"Panc1"]));
  
  dfPanc1=as.data.frame(Panc1);
  rownames(dfPanc1)=path_ACP_pathName[as.character(disease[ii,"ID"])];
  names(Panc1)=path_ACP_pathName[as.character(disease[ii,"ID"])];
 
  ### add results
  a=read.table("./data/results.txt",sep="\t",header = TRUE)
  com_res=(a[,"Compound.Name"])
  
  ii=which(com_res %in% rownames(Matrix_ACP))
  
  loewe=as.data.frame(a[ii,"Loewe.Score"])
  bliss=as.data.frame(a[ii,"Bliss.Score"])
  com_res=com_res[ii]
  rownames(bliss)=com_res
  rownames(loewe)=com_res
  jj=which(!(rownames(Matrix_ACP) %in% com_res));
  nn=length(jj);
  emp=as.data.frame(rep(0,nn))
  rownames(emp)=rownames(Matrix_ACP)[jj]
  bliss2=(c(as.vector(unlist(bliss)),as.vector(unlist(emp))))
  names(bliss2)=c(rownames(bliss),rownames(emp))
  bliss2=bliss2[rownames(Matrix_ACP)]
  
  loewe2=(c(as.vector(unlist(loewe)),as.vector(unlist(emp))))
  names(loewe2)=c(rownames(loewe),rownames(emp))
  loewe2=loewe2[rownames(Matrix_ACP)]
  
  an_bliss=columnAnnotation(Bliss=anno_barplot(bliss2,gp = gpar(fill= "#002060",fontsize = 8),ylim=c(0,26.7),axis_param =list(gp=gpar(fontsize = 8)))
                            , Loewe=anno_barplot(loewe2,gp = gpar(fill= "#C00000",fontsize = 8),ylim=c(0,51.313),axis_param =list(gp=gpar(fontsize = 8)))
                            
                            ,show_legend = TRUE
                            ,annotation_name_gp  = gpar(fontsize = 8,fontface="bold")
                            ,annotation_name_rot=45
                            ,height =unit(1.313,"cm")
                            ,gp = gpar(fontsize = 8)
  )
  
  ###
  col1_fun = colorRamp2(c(-8,-5, 0,5, 8), c(rgb(0,32,96,maxColorValue =255),rgb(0,176,240,maxColorValue =255),
                                            "white",rgb(255,72,72,maxColorValue =255), rgb(192,0,0,maxColorValue =255))) 
  ra=rowAnnotation(PANC1=anno_simple(Panc1,col=col1_fun)
                       ,name="PANC-1",show_annotation_name = TRUE
                       ,show_legend = TRUE
                       ,annotation_name_gp  = gpar(fontsize = 8,fontface="bold")                      
  )  

  Matrix_ACP=t(Matrix_ACP)
  w=ncol(Matrix_ACP) * 0.313 *1.313
  h= nrow(Matrix_ACP) * 0.313 * 1.313

 ### limit rowName length
 nameLength=26
 rnames2=str_sub(rownames(Matrix_ACP),1,nameLength) 
 ii=which(nchar(rownames(Matrix_ACP))>nameLength)
 added=rep("",nrow(Matrix_ACP))
 added[ii]="..."
 rnames2=paste0(rnames2,added)
 rownames(Matrix_ACP)=rnames2
 ##
  ht=Heatmap(Matrix_ACP,name = "Entichment"
             ,right_annotation = ra
             ,col = col1_fun
             ,heatmap_legend_param =list(legend_direction = "vertical"
              , title_gp = gpar(fontsize = 8, fontface = "bold"),labels_gp = gpar(fontsize = 8))
             
             ,width=unit(w,"cm")
             ,height =  unit(h,"cm")
             ,row_names_max_width= unit(4,"cm")
             ,row_dend_width = unit(7, "mm")
             ,column_dend_height = unit(7, "mm")
             
             ,row_names_gp = gpar(fontsize = 8,fontface="bold")
             ,column_names_gp = gpar(fontsize = 8,fontface="bold")
             
             ,cluster_rows = TRUE
             ,cluster_columns = TRUE
             ,column_title = "ACP Pathways"
             ,column_title_gp = gpar(fontsize=8,fontface="bold")
             ,show_heatmap_legend = TRUE
             ,top_annotation = an_bliss
             ,column_names_rot=45
             ,row_names_rot=45
  )
  filename=paste0("./figs/ACP.tiff")
  if(length(dev.list())!=0)dev.off()
  
  #tiff(filename = filename,width = (w+7.4),height = (h+5.4),units = "cm",res = 600)
  tiff(filename = filename,width = (21),height = (12),units = "cm",res = 600)
  draw(ht,heatmap_legend_side ="left")
  dev.off()
  return(ht);  
}

#################################################

#draws CP scores and saves it in the figs directory
draw_ht_score<-function(ScoreType,bshow_legend,compounds,path_CP,Matrix,MatIDs,disease){
  ii=which(compounds[,"Score.Type"]==ScoreType)
  compound_IDs_s=compounds[ii,"InstanceID"];
  compound_name_s=compounds[ii,"Compound"]
  names(compound_name_s)=compound_IDs_s;
  jj=which(path_CP[,"Score.Type"]==ScoreType)
  path_IDs_s=path_CP[jj,"ID"]
  path_names_s=path_CP[jj,"Pathway.Name"]
  path_names_s=str_replace(path_names_s,"-","")
  names(path_names_s)=path_IDs_s;
  
  MatIDs_s=which(MatIDs %in% compound_IDs_s)
  Matrix_s=as.matrix(Matrix[MatIDs_s,c("InstanceID",path_IDs_s)])
  rownames(Matrix_s)=as.character(Matrix_s[,1]);
  Matrix_s=Matrix_s[compound_IDs_s,]
  Matrix_s=cbind(compound_name_s[compound_IDs_s],Matrix_s)
  colnames(Matrix_s)[1]="compoundName"
  
  Matrix_s_n=Matrix_s[,as.character(path_IDs_s)];
  
  Matrix_s_n=apply(Matrix_s_n,2,as.numeric)
  rownames(Matrix_s_n)=Matrix_s[,"compoundName"]
  colnames(Matrix_s_n)=str_sub(path_names_s[colnames(Matrix_s_n)],1,70)
  
  ##
  ii=which(disease[,"ID"] %in% path_IDs_s)
  Panc1 = as.vector(unlist(disease[ii,"Panc1"]));
  
  dfPanc1=as.data.frame(Panc1);
  rownames(dfPanc1)=path_names_s[as.character(disease[ii,"ID"])];
  names(Panc1)=path_names_s[as.character(disease[ii,"ID"])];
  
  ### add results
  a=read.table("./data/results.txt",sep="\t",header = TRUE)
  com_res=(a[,"Compound.Name"])
  
  ii=which(com_res %in% rownames(Matrix_s_n))
  
  loewe=as.data.frame(a[ii,"Loewe.Score"])
  bliss=as.data.frame(a[ii,"Bliss.Score"])
  com_res=com_res[ii]
  rownames(bliss)=com_res
  rownames(loewe)=com_res
  jj=which(!(rownames(Matrix_s_n) %in% com_res));
  nn=length(jj);
  emp=as.data.frame(rep(0,nn))
  rownames(emp)=rownames(Matrix_s_n)[jj]
  bliss2=(c(as.vector(unlist(bliss)),as.vector(unlist(emp))))
  names(bliss2)=c(rownames(bliss),rownames(emp))
  bliss2=bliss2[rownames(Matrix_s_n)]
  
  
  loewe2=(c(as.vector(unlist(loewe)),as.vector(unlist(emp))))
  names(loewe2)=c(rownames(loewe),rownames(emp))
  loewe2=loewe2[rownames(Matrix_s_n)]
  
  an_bliss=columnAnnotation(Bliss=anno_barplot(bliss2,gp = gpar(fill= "#002060",fontsize = 8),ylim=c(0,26.7),axis_param =list(gp=gpar(fontsize = 8)))
                         , Loewe=anno_barplot(loewe2,gp = gpar(fill= "#C00000",fontsize = 8),ylim=c(0,51.313),axis_param =list(gp=gpar(fontsize = 8)))
                         
                         ,show_legend = TRUE
                         ,annotation_name_gp  = gpar(fontsize = 8,fontface="bold")
                         ,annotation_name_rot=45
                         ,height =unit(1.313,"cm")
                         ,gp = gpar(fontsize = 8)
  )
  ##
  
  col1_fun = colorRamp2(c(-8,-5, 0,5, 8), c(rgb(0,32,96,maxColorValue =255),rgb(0,176,240,maxColorValue =255), "white",rgb(255,72,72,maxColorValue =255), rgb(192,0,0,maxColorValue =255))) 
  ra=rowAnnotation(PANC1=anno_simple(Panc1,col=col1_fun)                       
                       ,show_legend = TRUE
                       ,annotation_name_gp  = gpar(fontsize = 8,fontface="bold")
  )
  
  Matrix_s_n=Matrix_s_n[,rownames(dfPanc1)]
  
  rownames(Matrix_s_n)=as.vector(rownames(Matrix_s_n))
  
  Matrix_s_n=t(Matrix_s_n);
  w=ncol(Matrix_s_n) * 0.313 *1.313
  h= nrow(Matrix_s_n) * 0.313 * 1.313
  title=paste0("CP Pathways: ",ScoreType)
  if(ScoreType=="Selected") title="Selected Pathways"
  ##
 ### reduce pathway name lengths
 nameLength=26
 rnames2=str_sub(rownames(Matrix_s_n),1,nameLength) 
 ii=which(nchar(rownames(Matrix_s_n))>nameLength)
 added=rep("",nrow(Matrix_s_n))
 added[ii]="..."
 rnames2=paste0(rnames2,added)
 rownames(Matrix_s_n)=rnames2
 ##

  ht=Heatmap(Matrix_s_n,name = "Entichment"
             ,right_annotation = ra
             ,top_annotation = an_bliss
             ,col = col1_fun
             ,heatmap_legend_param =list(legend_direction = "vertical", title_position = "topcenter", title_gp = gpar(fontsize = 8, fontface = "bold"),labels_gp = gpar(fontsize = 8))
             
             ,width=unit(w,"cm")
             ,height =  unit(h,"cm")
             ,row_names_max_width= unit(4,"cm")
             ,row_dend_width = unit(7, "mm")
             ,column_dend_height = unit(7, "mm")
             
             ,row_names_gp = gpar(fontsize = 8,fontface="bold")
             ,column_names_gp = gpar(fontsize = 8,fontface="bold")
             ,column_names_rot=45
             ,row_names_rot=45
             ,cluster_rows = TRUE
             ,cluster_columns = TRUE
             ,column_title = title
             ,column_title_gp = gpar(fontsize=8,fontface="bold")
             ,show_heatmap_legend = bshow_legend
            
  )
  filename=paste0("./figs/",ScoreType,".tiff")
  if(length(dev.list())!=0)dev.off()
  if(bshow_legend)h=h+0.8

  if(ScoreType=="Score1"){
    tiff(filename = filename,width = 11,height = 8,units = "cm",res = 600)
    
  }else if(ScoreType== "Res-score"){
    tiff(filename = filename,width = 10,height = 8,units = "cm",res = 600)
  
   }else if(ScoreType== "Selected"){
    tiff(filename = filename,width = 8,height = 8,units = "cm",res = 600)
  }  
  else{
    tiff(filename = filename,width = 7,height = 8,units = "cm",res = 600)
  }
  draw(ht,heatmap_legend_side ="left")
  dev.off()
  return(ht);
}

#################################################

##Merge figures to generate final figure 2 and saves it in figs directory; 
merge_tifs<-function(){
  library("ggplot2")
  library("cowplot")
  library("magick")
  
  disPath="./figs/disease.tiff"
  Score1Path="./figs/Score1.tiff";
  Score2Path="./figs/Score2.tiff";
  ResScorePath="./figs/Res-score.tiff";
  selectedPath="./figs/Selected.tiff";
  ACP_Path="./figs/ACP.tiff";
  ##
  
  score1=ggdraw() +  draw_image(Score1Path);
  score2=ggdraw() +  draw_image(Score2Path);

  ResScore=ggdraw() +  draw_image(ResScorePath);
  selected=ggdraw() +  draw_image(selectedPath);
  
  ACP=ggdraw() +  draw_image(ACP_Path);
  disease=ggdraw() +  draw_image(disPath);
  
  ########## all figs
  pp=plot_grid(ACP, labels =  c("A"),nrow=1)
  pp2=plot_grid(score1,ResScore, labels =  c("B","D"),nrow=1,ncol=2)
  pp3=plot_grid(score2,selected,disease, labels =  c("C","E","F"),nrow=1)
  pp4=plot_grid(pp,pp2,pp3,selected,nrow=3,ncol=1,rel_heights=c(12,8,8)) 
  ###
ggsave(filename = "./figs/Fig2.tiff",plot = pp4,width = 21, height = 28, dpi = 600,units="cm", compression = "lzw")
 
}
##
##Check if all packages used in for this is installed and if not installs them. The packages required are: feather, ggplot2, ComplexHeatmap, stringr, circlize, cowplot, magick
check_Installed_Packages<-function(){
  packages=c("feather","ggplot2","ComplexHeatmap","stringr","circlize","cowplot","magick")
  for(pck in packages){
    if(!(pck %in% installed.packages())){
      if(pck=="ComplexHeatmap"){
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        
        BiocManager::install("ComplexHeatmap")
      }else{
        install.packages(pck);  
      }      
    }
  }
  ii=which(!(packages %in% installed.packages()))
  if(length(ii)==1){ ##not all packages are installed    
    print(paste0( packeges[ii] ," not installed properly! Try installing it again!" ))
    return(FALSE)    
  }else if(length(ii)>1){
     paste0(paste0(packeges[ii] ,collapse =" and " )," are not installed properly! Try installing them again!")
     return(FALSE)
  }
  return(TRUE)
}

#################################################
checked=check_Installed_Packages()
if(checked){
  MainPlot()
  merge_tifs()
}

