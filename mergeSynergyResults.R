mergeSynergyResults<-function(){
	synFinder=read.table("./data2/SynergyFinder_SynergyScores.csv",sep="\t",header=TRUE)
	combenefit=read.table("./data2/Combenefit_SynergyScores.csv",sep="\t",header=TRUE)

	merged=merge(combenefit,synFinder,by="compound1",all=TRUE)
	write.table(merged,file="./data2/merged_SynergyScores.csv",sep = "\t")
	
	xx=1
}

mergeSynergyResults()