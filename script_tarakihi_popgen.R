library(ape)
#read.dna
library(adegenet)
#DNAbin2genind
library(readxl)
#read excel

#####1.Read sequences####
read.dna("tar_alignment.fasta",format="fasta")->tar_dna 
tar_dna <- as.matrix(tar_dna)

read.dna("ktar_alignment.fasta",format="fasta")->ktar_dna 
ktar_dna <- as.matrix(ktar_dna)

read.dna("total_alignment.fasta",format="fasta")->tar_ktar_dna
tar_ktar_dna <- as.matrix(tar_ktar_dna)

####2.assign populations#####
tarind <- DNAbin2genind(tar_dna, polyThres=0)
ktarind <- DNAbin2genind(ktar_dna, polyThres=0)
tar_ktarind <- DNAbin2genind(tar_ktar_dna, polyThres=0)
#we retain only SNPs for which the second largest allele frequency is greater than 10%

variables_table <- read_excel("data matrix.xlsx")
variables_table_tar<-variables_table[-(c(129:143)),]
variables_table_ktar<-variables_table[(c(129:143)),]

#check if individuals match
variables_table_tar$Samples==indNames(tarind)
variables_table_ktar$Samples==indNames(ktarind)
variables_table$Samples==indNames(tar_ktarind)

pop(tarind)<-variables_table_tar$Location
tarpop <- genind2genpop(tarind)
pop(ktarind)<-variables_table_ktar$Location
ktarpop <- genind2genpop(ktarind)
pop(tar_ktarind)<-variables_table$Location
tar_ktarpop <- genind2genpop(tar_ktarind)

####3.rarefaction curve####
library(spider)
#haploAccum

haploAccum(tar_dna,permutations = 1000)->tar_accum #MS result was done with 10000
haploAccum(ktar_dna,permutations = 1000)->ktar_accum #MS results was done with 10000

plot(tar_accum)
plot(ktar_accum)

####4.fst bonferroni adjust for pairwise Fst results####
pvalues<-c(0.95996,0.41504,0.78223,0.10254,0.15332,0.05566,0.90332,0.96875,0.74121,0.13574,0.23340,0.16602,0.27832,0.00098,0.10840,0.76172,0.69531,0.67090,0.01953,0.56348,0.37500,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.69043,0.78516,0.68652,0.33691,0.90234,0.31738,0.72852,0.00000,0.42188,0.48730,0.22559,0.04395,0.45312,0.55078,0.40039,0.00000,0.69336,0.80566,0.95020,0.41602,0.01953,0.72461,0.25977,0.79004,0.00000,0.84473,0.58398,0.51562,0.82324,0.06055,0.22754,0.43359,0.01367,0.04102,0.00000,0.34375,0.22266,0.10547,0.80664,0.99609,0.30664,0.07324,0.89746,0.14160,0.37695,0.00000,0.68652,0.43848,0.79004,0.67773,0.58203,0.25586,0.56055,0.01855,0.62695,0.90625,0.85059,0.00000,0.72656,0.77148,0.65430,0.04297,0.32715,0.54297,0.90625,0.27930,0.21387,0.70605,0.05566,0.18750,0.00000,0.36523,0.38184,0.22852,0.74414,0.77930,0.12988)
length(which(pvalues<0.05))
which(pvalues<0.05)

p.adjust(pvalues,method="bonferroni")->pvalues_adjusted_bf
length(which(pvalues_adjusted_bf<0.05))
which(pvalues_adjusted_bf<0.05)

p.adjust(pvalues,method="holm")->pvalues_adjusted_holm
length(which(pvalues_adjusted_holm<0.05))
which(pvalues_adjusted_holm<0.05)

p.adjust(pvalues,method="fdr")->pvalues_adjusted_fdr
length(which(pvalues_adjusted_fdr<0.05))
which(pvalues_adjusted_fdr<0.05)
#we choose fdr because it is less stringent and more powerful

####5.1 PCA####

sum(is.na(tar_ktarind$tab)) #23 NAs corresponding to gaps
#remove NAs (=gaps)
tar_ktarind$tab[ , colSums(is.na(tar_ktarind$tab)) == 0]->tar_ktarind$tab
sum(is.na(tar_ktarind$tab))
#Allele presence absence data are extracted and NAs replaced using tab:
X <- tab(tar_ktarind, freq = TRUE)
class(X)
dim(X)
X[1:5,1:5]

variables_table <- read_excel("data matrix.xlsx")
#View(variables_table)

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 4)
barplot(pca1$eig[1:50], main = "PCA eigenvalues")
summary(pca1)

#calculate projected inertia#
axis1<-pca1$eig[1]/sum(pca1$eig)*100
axis1<-paste0("PC1: ",round(axis1,digits=2),"%")
axis2<-pca1$eig[2]/sum(pca1$eig)*100
axis2<-paste0("PC2: ",round(axis2,digits=2),"%")
axis3<-pca1$eig[3]/sum(pca1$eig)*100
axis3<-paste0("PC3: ",round(axis3,digits=2),"%")
axis4<-pca1$eig[4]/sum(pca1$eig)*100
axis4<-paste0("PC4: ",round(axis4,digits=2),"%")


library(ggplot2)
#### PCA FOR LOCATIONS

rgb(48,76,208,maxColorValue=255)->darkblue #FRDL
rgb(0,125,207,maxColorValue=255)->blue #SPWCSI
rgb(0,191,159,maxColorValue=255)->turquoise #NT
rgb(172,105,130,maxColorValue=255)->lightbrown #ENLD
rgb(182,84,0,maxColorValue=255)->brown #BPLE
rgb(125,69,0,maxColorValue=255)->darkbrown #SPEC
rgb(180,84,219,maxColorValue=255)->purple #HB
rgb(206,0,134,maxColorValue=255)->dark_pink #WAI
rgb(255,142,237,maxColorValue=255)->light_pink #WGTN
rgb(167,211,124,maxColorValue=255)->light_green #SPCC
rgb(192,198,2,maxColorValue=255)->olive #KAIK
rgb(0,187,82,maxColorValue=255)->green #CHCH
rgb(38,95,31,maxColorValue=255)->darkgreen #OTAG
rgb(255,135,78,maxColorValue=255)->orange #CHAT
rgb(190,0,70,maxColorValue=255)->red #KTAR

colors<-c("FRDL"=darkblue,"SPWCSI"=blue,
          "NT"=turquoise,"ENLD"=lightbrown,"BPLE"=brown,
          "SPEC"=darkbrown,"HB"=purple,"WAI"=dark_pink,
          "WGTN"=light_pink,"SPCC"=light_green,"KAIK"=olive,
          "CHCH"=green,"OTAG"=darkgreen,"CHAT"=orange,"KTAR"=red)

as.factor(variables_table$Location)->locations
plot <- ggplot()+ coord_fixed()+geom_hline(yintercept=0, col="darkgrey")+
  geom_vline(xintercept=0, col="darkgrey") 
PCA.dfs <- data.frame(pca1$li[,1:4],locations)

pca.plot12<-plot + geom_point(data=PCA.dfs[,1:2],shape=21,size=5,
                              color="black",aes(x=Axis1, y=Axis2,
                                                fill=locations))+
  stat_ellipse(aes(x=PCA.dfs$Axis1, y=PCA.dfs$Axis2,color=locations),
               type="norm",size=1)+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  labs(x=axis1, y=axis2)+
  ggtitle("Principal Component Analysis: TAR and King TAR")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black",fill=NA))
pca.plot12

#Axes 3:4
pca.plot34<-plot + geom_point(data=PCA.dfs[,3:4],shape=21,size=5,
                              color="black",aes(x=Axis3, y=Axis4,
                                                fill=locations))+
  stat_ellipse(aes(x=PCA.dfs$Axis3, y=PCA.dfs$Axis4,color=locations),
               type="norm",size=1)+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  labs(x=axis3, y=axis4)+
  ggtitle("Principal Component Analysis: TAR and King TAR")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), line = element_line(colour = "black"),panel.border=element_rect(colour = "black",fill=NA))
pca.plot34

##### 5.2 PCA WITHOUT KTAR ####

sum(is.na(tarind$tab)) #23 NAs corresponding to gaps
#remove NAs (=gaps)
tarind$tab[ , colSums(is.na(tarind$tab)) == 0]->tarind$tab
sum(is.na(tarind$tab))
#Allele presence absence data are extracted and NAs replaced using tab:
X <- tab(tarind, freq = TRUE)
class(X)
dim(X)
X[1:5,1:5]

pca2 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 4)
barplot(pca2$eig[1:50], main = "PCA eigenvalues")
summary(pca2)

#calculate projected inertia#
axis1<-pca2$eig[1]/sum(pca2$eig)*100
axis1<-paste0("PC1: ",round(axis1,digits=2),"%")
axis2<-pca2$eig[2]/sum(pca2$eig)*100
axis2<-paste0("PC2: ",round(axis2,digits=2),"%")
axis3<-pca2$eig[3]/sum(pca2$eig)*100
axis3<-paste0("PC3: ",round(axis3,digits=2),"%")
axis4<-pca2$eig[4]/sum(pca2$eig)*100
axis4<-paste0("PC4: ",round(axis4,digits=2),"%")

#### PCA FOR LOCATIONS

rgb(48,76,208,maxColorValue=255)->darkblue #FRDL
rgb(0,125,207,maxColorValue=255)->blue #SPWCSI
rgb(0,191,159,maxColorValue=255)->turquoise #NT
rgb(172,105,130,maxColorValue=255)->lightbrown #ENLD
rgb(182,84,0,maxColorValue=255)->brown #BPLE
rgb(125,69,0,maxColorValue=255)->darkbrown #SPEC
rgb(180,84,219,maxColorValue=255)->purple #HB
rgb(206,0,134,maxColorValue=255)->dark_pink #WAI
rgb(255,142,237,maxColorValue=255)->light_pink #WGTN
rgb(167,211,124,maxColorValue=255)->light_green #SPCC
rgb(192,198,2,maxColorValue=255)->olive #KAIK
rgb(0,187,82,maxColorValue=255)->green #CHCH
rgb(38,95,31,maxColorValue=255)->darkgreen #OTAG
rgb(255,135,78,maxColorValue=255)->orange #CHAT

colors<-c("FRDL"=darkblue,"SPWCSI"=blue,
          "NT"=turquoise,"ENLD"=lightbrown,"BPLE"=brown,
          "SPEC"=darkbrown,"HB"=purple,"WAI"=dark_pink,
          "WGTN"=light_pink,"SPCC"=light_green,"KAIK"=olive,
          "CHCH"=green,"OTAG"=darkgreen,"CHAT"=orange)

as.factor(variables_table_tar$Location)->locations
plot <- ggplot()+ coord_fixed()+geom_hline(yintercept=0, col="darkgrey")+
  geom_vline(xintercept=0, col="darkgrey") 
PCA.dfs <- data.frame(pca2$li[,1:4],locations)

pca.plot12<-plot + geom_point(data=PCA.dfs[,1:2],shape=21,size=5,
                              color="black",aes(x=Axis1, y=Axis2,
                                                fill=locations))+
  stat_ellipse(aes(x=PCA.dfs$Axis1, y=PCA.dfs$Axis2,color=locations),
               type="norm",size=1)+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  labs(x=axis1, y=axis2)+
  ggtitle("Principal Component Analysis: TAR")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black",fill=NA))
pca.plot12

#Axes 3:4
pca.plot34<-plot + geom_point(data=PCA.dfs[,3:4],shape=21,size=5,
                              color="black",aes(x=Axis3, y=Axis4,
                                                fill=locations))+
  stat_ellipse(aes(x=PCA.dfs$Axis3, y=PCA.dfs$Axis4,color=locations),
               type="norm",size=1)+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  labs(x=axis3, y=axis4)+
  ggtitle("Principal Component Analysis: TAR and King TAR")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), line = element_line(colour = "black"),panel.border=element_rect(colour = "black",fill=NA))
pca.plot34

#### 5.3 PCA PROJECT SIZES ####
variables_table_tar$`Fork length (mm)`->size

PCA.dfs <- data.frame(pca2$li[,1:4],size)

pca.plot12<-plot + geom_point(data=PCA.dfs[,1:2],shape=21,size=5,
                              color="black",aes(x=Axis1, y=Axis2,
                                                fill=size))+
  scale_fill_gradient(low="blue",high="red")+
  labs(x=axis1, y=axis2)+
  ggtitle("Principal Component Analysis: TAR size")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black",fill=NA))
pca.plot12

#Axes 3:4
pca.plot34<-plot + geom_point(data=PCA.dfs[,3:4],shape=21,size=5,
                              color="black",aes(x=Axis3, y=Axis4,
                                                fill=size))+
  scale_fill_gradient(low="blue",high="red")+
  labs(x=axis3, y=axis4)+
  ggtitle("Principal Component Analysis: TAR size")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), line = element_line(colour = "black"),panel.border=element_rect(colour = "black",fill=NA))
pca.plot34

#### 5.4 PCA PROJECT WEIGHT ####
variables_table_tar$`Total weigth (g)`->weight

PCA.dfs <- data.frame(pca2$li[,1:4],weight)

pca.plot12<-plot + geom_point(data=PCA.dfs[,1:2],shape=21,size=5,
                              color="black",aes(x=Axis1, y=Axis2,
                                                fill=weight))+
  scale_fill_gradient(low="blue",high="red")+
  labs(x=axis1, y=axis2)+
  ggtitle("Principal Component Analysis: TAR weight")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black",fill=NA))
pca.plot12

#Axes 3:4
pca.plot34<-plot + geom_point(data=PCA.dfs[,3:4],shape=21,size=5,
                              color="black",aes(x=Axis3, y=Axis4,
                                                fill=weight))+
  scale_fill_gradient(low="blue",high="red")+
  labs(x=axis3, y=axis4)+
  ggtitle("Principal Component Analysis: TAR weight")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), line = element_line(colour = "black"),panel.border=element_rect(colour = "black",fill=NA))
pca.plot34

#### 6 TEST DAPC ####
##assign populations again##
tarind <- DNAbin2genind(tar_dna, polyThres=0)
ktarind <- DNAbin2genind(ktar_dna, polyThres=0)
tar_ktarind <- DNAbin2genind(tar_ktar_dna, polyThres=0)
#we retain only SNPs for which the second largest allele frequency is greater than 10%
#check if individuals match
variables_table_tar$Samples==indNames(tarind)
variables_table_ktar$Samples==indNames(ktarind)
variables_table$Samples==indNames(tar_ktarind)
pop(tarind)<-variables_table_tar$Location
tarpop <- genind2genpop(tarind)
pop(ktarind)<-variables_table_ktar$Location
ktarpop <- genind2genpop(ktarind)
pop(tar_ktarind)<-variables_table$Location
tar_ktarpop <- genind2genpop(tar_ktarind)

#### 6.1 DAPC pre-defined populations KTAR####
dapc_tar_ktarind<-dapc(tar_ktarind, var.contrib = TRUE, scale = FALSE, 
                       n.pca = 30, n.da = nPop(tar_ktarind) - 1)

scatter(dapc_tar_ktarind, cell = 0, pch = 18:23, cstar = 0,
        mstree = TRUE, lwd = 2, lty = 2)


#### 6.2 DAPC pre-defined populations ####
dapc_tarind<-dapc(tarind, var.contrib = TRUE, scale = FALSE, n.pca = 30, 
                  n.da = nPop(tarind) - 1)

scatter(dapc_tarind, cell = 0, pch = 18:23, cstar = 0,
        mstree = TRUE, lwd = 2, lty = 2)

#### 6.3 DAPC not pre-defined populations ####
find.clusters(tarind,max.n.clust=40,n.pca=200,n.clust=15)->grp
#Keep all the PCs, there is no downside: 200 to put more
#we chose 15 clusters but could be more.
table(pop(tarind), grp$grp)

table.value(table(pop(tarind), grp$grp), col.lab=paste("inf", 1:15),
            row.lab=paste(levels(pop(tarind))))

dapc1 <- dapc(tarind, grp$grp,var.contrib = TRUE, scale = FALSE,
              n.pca = 50, n.da = nPop(tar_ktarind) - 1)
#We want to retain a minimum of PCs: keep 50
#Keep 14 eigenvalues

scatter(dapc1, cell = 0, pch = 18:23, cstar = 0,
        mstree = TRUE, lwd = 2, lty = 2)

#### 7 TEST ISOLATION BY DISTANCE ####
library("gdistance")
library("MASS")
##prepare shapefile and raster##
shapefile<-shapefile("shapefiles/NewZealand_Boundary.shp") 
plot(shapefile)
ext <- extent(-180, 180, -47.36973, -34.39339) 
#important to have extant -180 180 take chatham into account
gridsize <- 0.1
r <- raster(ext, res=gridsize)
rr <- rasterize(shapefile, r)
plot(rr)
plot(rr,ext=c(-178, -175, -47.36973, -34.39339)) #here is chatham!

dist.dna(tar_dna,model="T92")->dnadist
View(as.matrix(dnadist)[1:50,1:50])

#extract coordiantes and Put back the real Chatham coordinates
variables_table_tar[,10:11]->latlon
latlon$long[latlon$long==182.9865]<--177.014
View(latlon)

####7.1 crow flight####
pC <- as.matrix(latlon[c("long","lat")])

geoDist <- pointDistance(pC, longlat = TRUE,allpairs=T)
geoDist <- as.dist(geoDist)

as.matrix(geoDist)->m_geoDist
colnames(as.matrix(dnadist))->colnames(m_geoDist)
rownames(as.matrix(dnadist))->rownames(m_geoDist)

View(m_geoDist[1:100,1:100])
#distances are in meters! 

mantel.rtest(geoDist,dnadist,nrepet=999)
#Observation: 0.02752448 , #p-value: 0.153
#plot(geoDist,dnadist)
#abline(lm(dnadist~geoDist), col="red",lty=2)
dens <- kde2d(geoDist,dnadist, n=150)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

plot(geoDist, dnadist,pch=20, 
     xlab="geo distance (meters)", ylab="genetic distance (T92)",cex.lab=1.3)
image(dens, col=transp(myPal(150),.7), add=TRUE)
abline(lm(dnadist~geoDist),lwd=3)
title("Isolation by distance: \"crow flight\"")


####7.2 restricted to ocean####
#move FRDL coord a bit to West so they are not considered enclave in land
latlon$long[latlon$long==166.9871]<--166.7

#Create a raster of ocean
rr->roce
roce[is.na(roce[])] <- 2
roce[roce[]==1] <- NA
roce[roce[]==2] <- 1
plot(roce)
plot(roce,ext=c(160, 180, -47.36973, -34.39339))

#do the distance analysis

troce <- transition(roce, mean, directions = 8) #less than 1 min
troceC <- geoCorrection(troce, "c", scl = TRUE)
pC <- as.matrix(latlon[c("long","lat")])
cosDist <- costDistance(troceC, pC)

as.matrix(cosDist)->m_cosDist
colnames(as.matrix(dnadist))->colnames(m_cosDist)
rownames(as.matrix(dnadist))->rownames(m_cosDist)
View(m_cosDist[1:100,1:100])

mantel.rtest(cosDist,dnadist,nrepet=999)
#Observation: 0.0271162, #p-value: 0.201
dens <- kde2d(cosDist,dnadist, n=150)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

plot(cosDist, dnadist,pch=20, 
     xlab="geographic least-cost distance", ylab="genetic distance (T92)",cex.lab=1.3)
image(dens, col=transp(myPal(150),.7), add=TRUE)
abline(lm(dnadist~cosDist),lwd=3)
title("Isolation by distance: restricted to ocean")

#### 7.3 restricted to ocean without Chatham Islands####

pC[-c(24:38),]->pC_nochat
tar_dna[-c(24:38),]->tar_dna_nochat

dist.dna(tar_dna_nochat,model="T92")->dnadist_nochat
View(as.matrix(dnadist_nochat)[1:50,1:50])

cosDist_nochat <- costDistance(troceC, pC_nochat)

as.matrix(cosDist_nochat)->m_cosDist_nochat
colnames(as.matrix(dnadist_nochat))->colnames(m_cosDist_nochat)
rownames(as.matrix(dnadist_nochat))->rownames(m_cosDist_nochat)
View(m_cosDist_nochat[1:100,1:100])

mantel.rtest(cosDist_nochat,dnadist_nochat,nrepet=999)
#Observation: 0.02419971, #p-value: 0.229
dens <- kde2d(cosDist_nochat,dnadist_nochat, n=150)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(cosDist_nochat, dnadist_nochat,pch=20, 
     xlab="geo least-cost distance (meters)", ylab="genetic distance (T92)",cex.lab=1.3)
image(dens, col=transp(myPal(150),.7), add=TRUE)
abline(lm(dnadist_nochat~cosDist_nochat),lwd=3)
title("Isolation by distance: \nleast-cost restricted to ocean, no Chatham")

#### 8 Mismatch distribution ####
library(ggplot2)
tar_mismatch_r_table <- read_excel("TAR_mismatch_table.xlsx")

plot1<-ggplot() +
  ylab("Frequency") + xlab("Pairwise Differences") +xlim(0,41) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data=tar_mismatch_r_table,aes(x=Differences,y=Observed)) +
  geom_line(data=tar_mismatch_r_table,aes(x=Differences,y=Observed),size=1) +
  geom_line(data=tar_mismatch_r_table,aes(x=Differences,y=Constant_Freq_Exp))
plot1

ktar_mismatch_r_table <- read_excel("KTAR_mismatch_table.xlsx")
plot1<-ggplot() +
  ylab("Frequency") + xlab("Pairwise Differences") +xlim(0,20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data=ktar_mismatch_r_table,aes(x=Differences,y=Observed)) +
  geom_line(data=ktar_mismatch_r_table,aes(x=Differences,y=Observed),size=1) +
  geom_line(data=ktar_mismatch_r_table,aes(x=Differences,y=Constant_Freq_Exp))
plot1