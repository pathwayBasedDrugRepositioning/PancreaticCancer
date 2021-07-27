
##
gemcitabineInstances<-function(){
    print(getwd())
    tbl= as.matrix(read.table("./data2/gemcitabine instances.txt",header=FALSE,sep=","))
    print(tbl)
    View(tbl)
    upTan=tanimato(tbl[2,(3:ncol(tbl))],tbl[4,(3:ncol(tbl))],"up")*100
    downTan=tanimato(tbl[1,(3:ncol(tbl))],tbl[3,(3:ncol(tbl))],"down")*100
    
    
}
##
tanimato<-function(vec1,vec2,dir){
    intersection= length(which(vec1 %in% vec2))
    union=length(unique(c(vec1,vec2)))
    
    tan=intersection/union

    print(paste0(dir," : intersection=",intersection, "  union=",union," tanimoto=",tan))
   
    return(tan)
}
##
gemcitabinePathwayInstances<-function(){
    library(feather)
    Matrix=as.matrix(read_feather("./data/LINCS_pathway_subset.feather"))    
    ii=grep("CPC006_A375_6H:BRD-K15108141-001-01-7:0.08_true",Matrix[,1])
    gem_inst1=as.numeric(Matrix[ii,(2:ncol(Matrix))])
    jj=grep("RAD001_MCF7_24H:BRD-K15108141-003-01-3:0.3704_false",Matrix[,1])
    gem_inst2=as.numeric(Matrix[jj,(2:ncol(Matrix))])
    View(gem_inst1)
    View(gem_inst2)
    View(Matrix)
    correlation=cor(gem_inst1,gem_inst2)
    print(correlation)
    xx=1
  
}
##
#gemcitabineGeneInstances()
gemcitabinePathwayInstances()