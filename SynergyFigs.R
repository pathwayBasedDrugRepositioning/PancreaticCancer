  library("ggplot2")
  library("cowplot")
  library("magick")

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
  label=paste0(drugName, " vs Gemcitabine on ",cell) 
  ########## all figs
  title <- ggdraw() + draw_label(label)
  pp=plot_grid(dose,fig_bliss,fig_HSA,fig_loewe,fig_zip,fig_loewe3d,nrow=1)
  pp2=plot_grid(fig_bliss3d,fig_loewe3d,fig_HSA3d,fig_zip3d,nrow=2,ncol=2)
  pp3=plot_grid(title,pp2,bar,nrow=3,rel_heights=c(0.1,1,1))
  ###
  #ggsave(filename = paste0("./figs/heat_",drugName,".tiff"),plot = pp,width = 28, height = 21, dpi = 600,units="cm", compression = "lzw")
  ggsave(filename = paste0("./figs/supplementaryFig_",label,".tiff"),plot = pp3,width = 21, height = 28, dpi = 600,units="cm", compression = "lzw")
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
  
  pp=plot_grid(f1,f2,f3,f4,f5,ncol=1,labels=LETTERS[1:5])
  pp2=plot_grid(f6,f7,f8,f9,f10,ncol=1,labels=LETTERS[1:5])

  ggsave(filename = paste0("./figs/fig3.tiff"),plot = pp,width = 11, height = 8.3, dpi = 600,units="in", compression = "lzw")
  ggsave(filename = paste0("./figs/fig4.tiff"),plot = pp2,width = 11, height = 8.3, dpi = 600,units="in", compression = "lzw")

}
#figs4drug("Entinostat")
gen_fig3()