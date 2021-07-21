library("readxl")
library(stringr)

mean2Mats<-function(mat1,mat2){
	mat=matrix(NA,nrow=nrow(mat1),ncol=ncol(mat1))
	
	for(i in 1:nrow(mat1)){
		for(j in 1:ncol(mat1)){
			mat[i,j]=mean(c(as.numeric(mat1[i,j]),as.numeric(mat2[i,j])))
		}
	}
	rownames(mat)=rownames(mat1)
	colnames(mat)=colnames(mat1)
	return(mat)
}
##
mean3Mats<-function(mat1,mat2,mat3){
	mat=matrix(NA,nrow=nrow(mat1),ncol=ncol(mat1))
	
	for(i in 1:nrow(mat1)){
		for(j in 1:ncol(mat1)){
			mat[i,j]=mean(c(as.numeric(mat1[i,j]),as.numeric(mat2[i,j]),as.numeric(mat3[i,j])))
		}
	}
	rownames(mat)=rownames(mat1)
	colnames(mat)=colnames(mat1)
	return(mat)
}
##
writeOutput4SynergyFinder<-function(file,matrix,drug1,drug2,cell,unit="nM"){
	header=c("Drug1:",drug1)
	header=rbind(header,c("Drug2:",drug2))
	header=rbind(header,c("Cell:",cell))
	header=rbind(header,c("ConcUnit:",unit))
	matrix=rbind(colnames(matrix),matrix)
	matrix=cbind(rownames(matrix),matrix)

	addedHeader=matrix(NA,nrow=nrow(header),ncol=(ncol(matrix)-ncol(header)))
	header=cbind(header,addedHeader)

	matrix=rbind(header,matrix)	
	write.table(matrix,file,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,na="")
	
}

##
reformatSynergy<-function(){
	dirs=dir("./rawSynergyData/")
	i=0
	for(directory in dirs){
		i=i+1
		files=dir(paste0("./rawSynergyData/",directory),full.names = TRUE)
		file=files[1]
		tbl=read_excel(file)
		tbl1=as.matrix(tbl)			
		ii=which(tbl1[,1]=="Agent 1")
		nrows1=ii-2
		ncols1=(ncol(tbl1))-1
		matrix1=tbl1[1:nrows1,2:ncols1]
		rownames(matrix1)=tbl1[1:(ii-2),1]
		drug1=tbl1[ii,2]
		drug2=tbl1[ii+1,2]
		title=tbl1[ii+4,2]
		print(directory)	
		if(drug1=="Entinostat")
		{
			cell=str_split(directory,"  ")[[1]][2]
		}else{
			cell="PANC-1"
		}
		
		file=files[2]
		tbl=read_excel(file)
		tbl1=as.matrix(tbl)			
		matrix2=tbl1[1:nrows1,2:ncols1]
		rownames(matrix2)=tbl1[1:(ii-2),1]
			
		if(length(files)==2){
			matrix=mean2Mats(matrix1,matrix2)
		}
		if(length(files)==3){
			file=files[3]
			tbl=read_excel(file)
			tbl1=as.matrix(tbl)			
			matrix3=tbl1[1:nrows1,2:ncols1]
			rownames(matrix3)=tbl1[1:(ii-2),1]
			matrix=mean3Mats(matrix1,matrix2,matrix3)
			}		
				
		writeOutput4SynergyFinder("data2/synergyAll.txt",matrix=matrix,drug1=drug1,drug2=drug2,cell = cell)
		xx=1
	}
	print("Happy coding")
}



reformatSynergy()