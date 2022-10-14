## >>>>>>>>>>>> Script for Figure 4 <<<<<<<<<<<

rm(list=ls())

library(ggplot2)

load("res_chr1A.Rdata")
load("res_chr2A.Rdata")
load("res_chr3A.Rdata")
load("res_chr4A.Rdata")
load("res_chr5A.Rdata")
load("res_chr6A.Rdata")
load("res_chr7A.Rdata")
load("res_chr1B.Rdata")
load("res_chr2B.Rdata")
load("res_chr3B.Rdata")
load("res_chr4B.Rdata")
load("res_chr5B.Rdata")
load("res_chr6B.Rdata")
load("res_chr7B.Rdata")
load("res_chr1D.Rdata")
load("res_chr2D.Rdata")
load("res_chr3D.Rdata")
load("res_chr4D.Rdata")
load("res_chr5D.Rdata")
load("res_chr6D.Rdata")
load("res_chr7D.Rdata")

data_chromosome <- read.csv2("carte_Axiom-TABW420k_WGAv1.csv",header=TRUE, sep = " ") #dataframe containing the position on the genome of each SNP
#data_chromosome$V1: chromosome, data_chromosome$V2: position on the chromosome, data_chromosome$name: name of the SNP

heading_QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]      #list of heading QTLs names
heading_QTL <- as.character(heading_QTL)

flowering_genes <- read.csv("GenesMajeursFloraison.csv",sep=";")[c(1:3,7:12),]   #dataframe containing the position on the genome of each flowering gene

## List of markers selected by SAEMVS on chr1A
res_chr1A

pos_select1A=data_chromosome$V2[which(data_chromosome$name%in%res_chr1A)] #position of selected SNPs 
pos_QTL1A=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr1A"] #position of heading QTL
abs(pos_select1A-pos_QTL1A)<=10^6               #Boolean determining whether QTLs of chr1A has been identified by SAEMVS                                  

## List of markers selected by SAEMVS on chr2A
res_chr2A

pos_select2A=data_chromosome$V2[which(data_chromosome$name%in%res_chr2A)] #position of selected SNPs
pos_gene2A=flowering_genes$Position[which(flowering_genes$Chromosome=="2A")] #position of flowering genes
abs(pos_select2A-pos_gene2A)<=10^6         #Boolean determining whether flowering genes of chr2A has been identified by SAEMVS

## LList of markers selected by SAEMVS on chr3A
res_chr3A

## List of markers selected by SAEMVS on chr4A
res_chr4A


## List of markers selected by SAEMVS on chr5A
res_chr5A

pos_select5A=data_chromosome$V2[which(data_chromosome$name%in%res_chr5A)]
pos_gene5A=flowering_genes$Position[which(flowering_genes$Chromosome=="5A")]
abs(pos_select5A-pos_gene5A)<=10^6

## List of markers selected by SAEMVS on chr6A
res_chr6A

pos_select6A=data_chromosome$V2[which(data_chromosome$name%in%res_chr6A)]
pos_QTL6A=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr6A"]
abs(pos_select6A-pos_QTL6A)<=10^6

## List of markers selected by SAEMVS on chr7A
res_chr7A

pos_select7A=data_chromosome$V2[which(data_chromosome$name%in%res_chr7A)]
pos_gene7A=flowering_genes$Position[which(flowering_genes$Chromosome=="7A")]
abs(pos_select7A-pos_gene7A)<=10^6

## List of markers selected by SAEMVS on chr1B
res_chr1B

## List of markers selected by SAEMVS on chr2B
res_chr2B

pos_select2B=data_chromosome$V2[which(data_chromosome$name%in%res_chr2B)]
pos_QTL2B=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr2B"]
abs(pos_select2B-pos_QTL2B)<=10^6

pos_gene2B=flowering_genes$Position[which(flowering_genes$Chromosome=="2B")]
abs(pos_select2B-pos_gene2B)<=10^6

## List of markers selected by SAEMVS on chr3B
res_chr3B

## List of markers selected by SAEMVS on chr4B
res_chr4B

## List of markers selected by SAEMVS on chr5B
res_chr5B

pos_select5B=data_chromosome$V2[which(data_chromosome$name%in%res_chr5B)]
pos_gene5B=flowering_genes$Position[which(flowering_genes$Chromosome=="5B")]
abs(pos_select5B-pos_gene5B)<=10^6

## List of markers selected by SAEMVS on chr6B
res_chr6B

## List of markers selected by SAEMVS on chr7B
res_chr7B

pos_select7B=data_chromosome$V2[which(data_chromosome$name%in%res_chr7B)]
pos_QTL7B=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr7B"]
abs(pos_select7B-pos_QTL7B)<=10^6

pos_gene7B=flowering_genes$Position[which(flowering_genes$Chromosome=="7B")]
abs(pos_select7B-pos_gene7B)<=10^6

## List of markers selected by SAEMVS on chr1D
res_chr1D

## List of markers selected by SAEMVS on chr2D
res_chr2D

pos_select2D=data_chromosome$V2[which(data_chromosome$name%in%res_chr2D)]
pos_QTL2D=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr2D"]
abs(pos_select2D[1]-pos_QTL2D)<=10^6
abs(pos_select2D[2]-pos_QTL2D)<=10^6

pos_gene2D=flowering_genes$Position[which(flowering_genes$Chromosome=="2D")]
abs(pos_select2D[1]-pos_gene2D)<=10^6
abs(pos_select2D[2]-pos_gene2D)<=10^6

## List of markers selected by SAEMVS on chr3D
res_chr3D

## List of markers selected by SAEMVS on chr4D
res_chr4D

## List of markers selected by SAEMVS on chr5D
res_chr5D

pos_select5D=data_chromosome$V2[which(data_chromosome$name%in%res_chr5D)]
pos_gene5D=flowering_genes$Position[which(flowering_genes$Chromosome=="5D")]
abs(pos_select5D-pos_gene5D)<=10^6

## List of markers selected by SAEMVS on chr6D
res_chr6D

## List of markers selected by SAEMVS on chr7D
res_chr7D

pos_select7D=data_chromosome$V2[which(data_chromosome$name%in%res_chr7D)]
pos_gene7D=flowering_genes$Position[which(flowering_genes$Chromosome=="7D")]
abs(pos_select7D-pos_gene7D)<=10^6


################## Script for Figure 4 ##################

data_chromosome <- read.csv2("carte_Axiom-TABW420k_WGAv1.csv",header=TRUE, sep = " ")

pos_select1A=data_chromosome$V2[which(data_chromosome$name%in%res_chr1A)]
pos_select1B=data_chromosome$V2[which(data_chromosome$name%in%res_chr1B)]
pos_select1D=data_chromosome$V2[which(data_chromosome$name%in%res_chr1D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL1A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr1A"]
pos_QTL1B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr1B"]
pos_QTL1D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr1D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf1A=genes_floraison$Position[which(genes_floraison$Chromosome=="1A")]
pos_gf1B=genes_floraison$Position[which(genes_floraison$Chromosome=="1B")]
pos_gf1D=genes_floraison$Position[which(genes_floraison$Chromosome=="1D")]

data_chr1 <- data.frame(chr=c(rep("chr1A",length(pos_select1A)+1),rep("chr1B",length(pos_select1B)+1),rep("chr1D",length(pos_select1D)+1)),pos_select=c(c(NA,pos_select1A),c(NA,pos_select1B),c(NA,pos_select1D)))

data_chr1_QTL <- data.frame(chr=c(rep("chr1A",length(pos_QTL1A)),rep("chr1B",length(pos_QTL1B)),rep("chr1D",length(pos_QTL1D))),pos=c(pos_QTL1A,pos_QTL1B,pos_QTL1D))

data_chr1_flor <- data.frame(chr=c(rep("chr1A",length(pos_gf1A)),rep("chr1B",length(pos_gf1B)),rep("chr1D",length(pos_gf1D))),pos=c(pos_gf1A,pos_gf1B,pos_gf1D))

pos_select2A=data_chromosome$V2[which(data_chromosome$name%in%res_chr2A)]
pos_select2B=data_chromosome$V2[which(data_chromosome$name%in%res_chr2B)]
pos_select2D=data_chromosome$V2[which(data_chromosome$name%in%res_chr2D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL2A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr2A"]
pos_QTL2B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr2B"]
pos_QTL2D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr2D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf2A=genes_floraison$Position[which(genes_floraison$Chromosome=="2A")]
pos_gf2B=genes_floraison$Position[which(genes_floraison$Chromosome=="2B")]
pos_gf2D=genes_floraison$Position[which(genes_floraison$Chromosome=="2D")]

data_chr2 <- data.frame(chr=c(rep("chr2A",length(pos_select2A)+1),rep("chr2B",length(pos_select2B)+1),rep("chr2D",length(pos_select2D)+1)),pos_select=c(c(NA,pos_select2A),c(NA,pos_select2B),c(NA,pos_select2D)))

data_chr2_QTL <- data.frame(chr=c(rep("chr2A",length(pos_QTL2A)),rep("chr2B",length(pos_QTL2B)),rep("chr2D",length(pos_QTL2D))),pos=c(pos_QTL2A,pos_QTL2B,pos_QTL2D))

data_chr2_flor <- data.frame(chr=c(rep("chr2A",length(pos_gf2A)),rep("chr2B",length(pos_gf2B)),rep("chr2D",length(pos_gf2D))),pos=c(pos_gf2A,pos_gf2B,pos_gf2D))

pos_select3A=data_chromosome$V2[which(data_chromosome$name%in%res_chr3A)]
pos_select3B=data_chromosome$V2[which(data_chromosome$name%in%res_chr3B)]
pos_select3D=data_chromosome$V2[which(data_chromosome$name%in%res_chr3D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL3A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr3A"]
pos_QTL3B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr3B"]
pos_QTL3D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr3D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf3A=genes_floraison$Position[which(genes_floraison$Chromosome=="3A")]
pos_gf3B=genes_floraison$Position[which(genes_floraison$Chromosome=="3B")]
pos_gf3D=genes_floraison$Position[which(genes_floraison$Chromosome=="3D")]

data_chr3 <- data.frame(chr=c(rep("chr3A",length(pos_select3A)+1),rep("chr3B",length(pos_select3B)+1),rep("chr3D",length(pos_select3D)+1)),pos_select=c(c(NA,pos_select3A),c(NA,pos_select3B),c(NA,pos_select3D)))

data_chr3_QTL <- data.frame(chr=c(rep("chr3A",length(pos_QTL3A)),rep("chr3B",length(pos_QTL3B)),rep("chr3D",length(pos_QTL3D))),pos=c(pos_QTL3A,pos_QTL3B,pos_QTL3D))

data_chr3_flor <- data.frame(chr=c(rep("chr3A",length(pos_gf3A)),rep("chr3B",length(pos_gf3B)),rep("chr3D",length(pos_gf3D))),pos=c(pos_gf3A,pos_gf3B,pos_gf3D))


pos_select4A=data_chromosome$V2[which(data_chromosome$name%in%res_chr4A)]
pos_select4B=data_chromosome$V2[which(data_chromosome$name%in%res_chr4B)]
pos_select4D=data_chromosome$V2[which(data_chromosome$name%in%res_chr4D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL4A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr4A"]
pos_QTL4B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr4B"]
pos_QTL4D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr4D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf4A=genes_floraison$Position[which(genes_floraison$Chromosome=="4A")]
pos_gf4B=genes_floraison$Position[which(genes_floraison$Chromosome=="4B")]
pos_gf4D=genes_floraison$Position[which(genes_floraison$Chromosome=="4D")]

data_chr4 <- data.frame(chr=c(rep("chr4A",length(pos_select4A)+1),rep("chr4B",length(pos_select4B)+1),rep("chr4D",length(pos_select4D)+1)),pos_select=c(c(NA,pos_select4A),c(NA,pos_select4B),c(NA,pos_select4D)))

data_chr4_QTL <- data.frame(chr=c(rep("chr4A",length(pos_QTL4A)),rep("chr4B",length(pos_QTL4B)),rep("chr4D",length(pos_QTL4D))),pos=c(pos_QTL4A,pos_QTL4B,pos_QTL4D))

data_chr4_flor <- data.frame(chr=c(rep("chr4A",length(pos_gf4A)),rep("chr4B",length(pos_gf4B)),rep("chr4D",length(pos_gf4D))),pos=c(pos_gf4A,pos_gf4B,pos_gf4D))


pos_select5A=data_chromosome$V2[which(data_chromosome$name%in%res_chr5A)]
pos_select5B=data_chromosome$V2[which(data_chromosome$name%in%res_chr5B)]
pos_select5D=data_chromosome$V2[which(data_chromosome$name%in%res_chr5D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL5A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr5A"]
pos_QTL5B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr5B"]
pos_QTL5D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr5D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf5A=genes_floraison$Position[which(genes_floraison$Chromosome=="5A")]
pos_gf5B=genes_floraison$Position[which(genes_floraison$Chromosome=="5B")]
pos_gf5D=genes_floraison$Position[which(genes_floraison$Chromosome=="5D")]

data_chr5 <- data.frame(chr=c(rep("chr5A",length(pos_select5A)+1),rep("chr5B",length(pos_select5B)+1),rep("chr5D",length(pos_select5D)+1)),pos_select=c(c(NA,pos_select5A),c(NA,pos_select5B),c(NA,pos_select5D)))

data_chr5_QTL <- data.frame(chr=c(rep("chr5A",length(pos_QTL5A)),rep("chr5B",length(pos_QTL5B)),rep("chr5D",length(pos_QTL5D))),pos=c(pos_QTL5A,pos_QTL5B,pos_QTL5D))

data_chr5_flor <- data.frame(chr=c(rep("chr5A",length(pos_gf5A)),rep("chr5B",length(pos_gf5B)),rep("chr5D",length(pos_gf5D))),pos=c(pos_gf5A,pos_gf5B,pos_gf5D))


pos_select6A=data_chromosome$V2[which(data_chromosome$name%in%res_chr6A)]
pos_select6B=data_chromosome$V2[which(data_chromosome$name%in%res_chr6B)]
pos_select6D=data_chromosome$V2[which(data_chromosome$name%in%res_chr6D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL6A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr6A"]
pos_QTL6B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr6B"]
pos_QTL6D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr6D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
genes_floraison=genes_floraison_taille[c(1:3,7:12),]
pos_gf6A=genes_floraison$Position[which(genes_floraison$Chromosome=="6A")]
pos_gf6B=genes_floraison$Position[which(genes_floraison$Chromosome=="6B")]
pos_gf6D=genes_floraison$Position[which(genes_floraison$Chromosome=="6D")]

data_chr6 <- data.frame(chr=c(rep("chr6A",length(pos_select6A)+1),rep("chr6B",length(pos_select6B)+1),rep("chr6D",length(pos_select6D)+1)),pos_select=c(c(NA,pos_select6A),c(NA,pos_select6B),c(NA,pos_select6D)))

data_chr6_QTL <- data.frame(chr=c(rep("chr6A",length(pos_QTL6A)),rep("chr6B",length(pos_QTL6B)),rep("chr6D",length(pos_QTL6D))),pos=c(pos_QTL6A,pos_QTL6B,pos_QTL6D))

data_chr6_flor <- data.frame(chr=c(rep("chr6A",length(pos_gf6A)),rep("chr6B",length(pos_gf6B)),rep("chr6D",length(pos_gf6D))),pos=c(pos_gf6A,pos_gf6B,pos_gf6D))


pos_select7A=data_chromosome$V2[which(data_chromosome$name%in%res_chr7A)]
pos_select7B=data_chromosome$V2[which(data_chromosome$name%in%res_chr7B)]
pos_select7D=data_chromosome$V2[which(data_chromosome$name%in%res_chr7D)]

QTL <- read.csv("marker_HD_INRmon13LN.csv",sep="")[,1]
QTL <- as.character(QTL)
pos_QTL7A=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr7A"]
pos_QTL7B=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr7B"]
pos_QTL7D=data_chromosome$V2[data_chromosome$name%in%QTL & data_chromosome$V1=="chr7D"]

genes_floraison_taille <- read.csv("GenesMajeursFloraison.csv",sep=";")
pos_gf7A=genes_floraison$Position[which(genes_floraison$Chromosome=="7A")]
pos_gf7B=genes_floraison$Position[which(genes_floraison$Chromosome=="7B")]
pos_gf7D=genes_floraison$Position[which(genes_floraison$Chromosome=="7D")]

data_chr7 <- data.frame(chr=c(rep("chr7A",length(pos_select7A)+1),rep("chr7B",length(pos_select7B)+1),rep("chr7D",length(pos_select7D)+1)),pos_select=c(c(NA,pos_select7A),c(NA,pos_select7B),c(NA,pos_select7D)))

data_chr7_QTL <- data.frame(chr=c(rep("chr7A",length(pos_QTL7A)),rep("chr7B",length(pos_QTL7B)),rep("chr7D",length(pos_QTL7D))),pos=c(pos_QTL7A,pos_QTL7B,pos_QTL7D))

data_chr7_flor <- data.frame(chr=c(rep("chr7A",length(pos_gf7A)),rep("chr7B",length(pos_gf7B)),rep("chr7D",length(pos_gf7D))),pos=c(pos_gf7A,pos_gf7B,pos_gf7D))


data_join_chr= rbind(data_chr1, data_chr2, data_chr3, data_chr4, data_chr5, data_chr6, data_chr7)
data_join_QTL=rbind(data_chr1_QTL, data_chr2_QTL, data_chr3_QTL, data_chr4_QTL, data_chr5_QTL, data_chr6_QTL, data_chr7_QTL)
data_join_flor=rbind(data_chr1_flor, data_chr2_flor, data_chr3_flor, data_chr4_flor, data_chr5_flor, data_chr6_flor, data_chr7_flor)

plot_chr_join=ggplot(data_join_chr,aes(x=chr,group=chr,y=pos_select,shape="Selected SNPs"))+ geom_point(size=6) + geom_point(data = data_join_QTL,aes(y=pos,group=chr,color="Heading QTLs"),pch=18,size=5,position = position_dodge(width = 0.9)) + geom_point(data = data_join_flor,aes(y=pos,group=chr,color="Major flowering genes"),pch=18,size=5,position = position_dodge(width = 0.9))+
  labs(x="Chromosome", y = "Position") + scale_shape_manual(name= " ", values = c("Selected SNPs"=4))+ scale_color_manual(name = " ",values = c("Heading QTLs"="red","Major flowering genes"="darkgreen")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"),legend.title = element_text(size=20,face="bold"),legend.text = element_text(size=15))+ theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "grey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  )

save(plot_chr_join,file="plot_chr_join.Rdata")
plot_chr_join
