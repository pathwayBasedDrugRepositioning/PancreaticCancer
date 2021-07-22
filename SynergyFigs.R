  library("ggplot2")
  library("cowplot")
  library("magick")
  library(stringr)
figs4drug<-function(drugName,cell="PANC-1"){
  
  
  folder=paste0("./SynergyFinderOutput/",drugName,"/")
  
  blissH=paste0(folder,"HBliss.svg");
  HSAH=paste0(folder,"HHSA.svg")
  loeweH=paste0(folder,"HLoewe.svg")
  zipH = paste0(folder,"HZIP.svg")  

  bliss3d=paste0(folder,"Bliss.svg");
  HSA3d=paste0(folder,"HSA.svg")
  loewe3d=paste0(folder,"Loewe.svg")
  zip3d = paste0(folder,"ZIP.svg")

  dose=paste0(folder,"newplot.svg")
  bar=dir(folder,"SynergyFinder_bar_plot",full.names=TRUE)
  ##
  
  fig_bliss=ggdraw() +  draw_image(blissH);
  fig_HSA=ggdraw() +  draw_image(HSAH);
  fig_loewe=ggdraw() +  draw_image(loeweH);
  fig_zip=ggdraw() +  draw_image(zipH);
  
  fig_bliss3d=ggdraw() +  draw_image(bliss3d);
  fig_HSA3d=ggdraw() +  draw_image(HSA3d);
  fig_loewe3d=ggdraw() +  draw_image(loewe3d);
  fig_zip3d=ggdraw() +  draw_image(zip3d);
  
  dose=ggdraw() +  draw_image(dose);
  bar=ggdraw() +  draw_image(bar);
  drugName=str_split(drugName,"_")[[1]][1]
  label=paste0(drugName, " vs Gemcitabine on ",cell) 
  ########## all figs
  title <- ggdraw() + draw_label(label)
  pp=plot_grid(dose,fig_bliss,fig_HSA,fig_loewe,fig_zip,fig_loewe3d,nrow=1)
  pp2=plot_grid(fig_bliss3d,fig_loewe3d,fig_HSA3d,fig_zip3d,nrow=2,ncol=2)
  pp3=plot_grid(title,pp2,bar,nrow=3,rel_heights=c(0.1,1,1))
  ###
  ggsave(filename = paste0("./figs/supplementaryFig_",label,".tiff"),plot = pp3,width = 21, height = 28, dpi = 600,units="cm", compression = "lzw")
  
  ggsave(filename = paste0("./figs/supplementaryFig_",label,".pdf"),plot = pp3,width = 21, height = 28, dpi = 600,units="cm")
  return(pp)
}
gen_fig3<-function(){
  f1=figs4drug("Entinostat")
  f2=figs4drug("Loperamide")
  f3=figs4drug("Thioridazine")
  f4=figs4drug("Saracatinib")
  f5=figs4drug("Scriptaid")
  f6=figs4drug("Palbociclib")
  f7=figs4drug("Racecadotril")
  f8=figs4drug("STK525924")
  f9=figs4drug("BX795")
  f10=figs4drug("Semagacestat")
  
  f11=figs4drug("Entinostat_HPAF-II",cell="HPAF-II")
  f12=figs4drug("Entinostat_K8484",cell="K8484")
  f13=figs4drug("Entinostat_MIA PaCa-2",cell="MIA PaCa-2")
  f14=figs4drug("Entinostat_TB32048",cell="TB32048")
  
  pp=plot_grid(f1,f2,f3,f4,f5,ncol=1,labels=LETTERS[1:5])
  pp2=plot_grid(f6,f7,f8,f9,f10,ncol=1,labels=LETTERS[1:5])
  pp3=plot_grid(f11,f12,f13,f14,f1,ncol=1,labels=LETTERS[1:5])

  ggsave(filename = paste0("./figs/fig3.tiff"),plot = pp,width = 11, height = 8.3, dpi = 600,units="in", compression = "lzw")
  ggsave(filename = paste0("./figs/fig4.tiff"),plot = pp2,width = 11, height = 8.3, dpi = 600,units="in", compression = "lzw")
  ggsave(filename = paste0("./figs/SupplementaryFig_Entinostat.tiff"),plot = pp3,width = 11, height = 8.3, dpi = 600,units="in", compression = "lzw")
  
  ggsave(filename = paste0("./figs/fig3.pdf"),plot = pp,width = 11, height = 8.3, dpi = 600,units="in")
  ggsave(filename = paste0("./figs/fig4.pdf"),plot = pp2,width = 11, height = 8.3, dpi = 600,units="in")
  ggsave(filename = paste0("./figs/SupplementaryFig_Entinostat.pdf"),plot = pp3,width = 11, height = 8.3, dpi = 600,units="in")

}
gen_fig5<-function(){
 folder=paste0("./SynergyFinderOutput/All/")
  
  bliss=paste0(folder,"Bliss.svg");
  HSA=paste0(folder,"HSA.svg")
  loewe=paste0(folder,"Loewe.svg")
  zip = paste0(folder,"ZIP.svg")  

  fig_bliss=ggdraw() +  draw_image(bliss);
  fig_HSA=ggdraw() +  draw_image(HSA);
  fig_loewe=ggdraw() +  draw_image(loewe);
  fig_zip=ggdraw() +  draw_image(zip);
  
  pp=plot_grid(fig_bliss,fig_loewe,ncol=1,labels=LETTERS[1:2])
  pp2=plot_grid(fig_HSA,fig_zip,ncol=1,labels=LETTERS[1:2])
 
  ggsave(filename = paste0("./figs/fig5.tiff"),plot = pp,width = 8, height = 11, dpi = 600,units="in", compression = "lzw")
  ggsave(filename = paste0("./figs/SupplementaryFig_synergy.tiff"),plot = pp2,width = 8, height = 11, dpi = 600,units="in", compression = "lzw")
  
  ggsave(filename = paste0("./figs/fig5.pdf"),plot = pp,width = 8, height = 11, dpi = 600,units="in")
  ggsave(filename = paste0("./figs/SupplementaryFig_synergy.pdf"),plot = pp2,width = 8, height = 11, dpi = 600,units="in")

  
  xx=1
}

#figs4drug("Entinostat")
gen_fig3()
gen_fig5()
print("Happy")