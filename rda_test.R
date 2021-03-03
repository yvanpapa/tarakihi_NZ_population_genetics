### import adegenet script env

#Script to perform a Redundancy Analysis (RDA) using spatial data, obtained from the paper:
#Meirmans, P.G. (2015) Seven common mistakes in population genetics and how to avoid them; Molecular Ecology.

#need environment from popgen script to import all required objects

library(hierfstat)
library(vegan)
library(SoDA)
library(reshape)

library("gdistance")
library("MASS")

library(ape)
#read.dna
library(adegenet)
#DNAbin2genind
library(readxl)
#read excel

#We need to replace population names by numbers
as.factor(variables_table_tar$Location)->locs
un<-unique(locs)
c()->Loc_num

for (i in 1:length(un)) {
Loc_num[which(locs==un[i])]<-i
}

pop(tarind)<-Loc_num

genind2hierfstat(tarind,pop=NULL)->dat


#custom functions to calculate allele frequencies (modified from script by Jerome Goudet)
"getal" <- function(dats)
{
  x<-dim(dats)
  if (max(dats[,2],na.rm=T)<1000000) modulo=1000
  if (max(dats[,2],na.rm=T)<10000) modulo<-100
  if (max(dats[,2],na.rm=T)<100) modulo<-10
  firstal<-dats[,-1] %/% modulo
  secal<-dats[,-1] %% modulo
  ind<-vector(length=0)
  nbpop<-max(dats[,1])
  for (i in 1:nbpop) {
    dum<-1:sum(dats[,1]==i)
    if (i==1) ind<-dum else ind<-c(ind,dum)
  }
  ind<-rep(ind,2)
  if (x[2]==2) dats.al<-data.frame(pop=rep(dats[,1],2),ind=ind,al=c(firstal,secal))
  else dats.al<-data.frame(pop=rep(dats[,1],2),ind=ind,al=rbind(firstal,secal))
  return(dats.al)
}


"pop.freqs"<-function(dats,diploid=TRUE)
{
  nbpop<-max(dats[,1])
  if (diploid) dats<-getal(dats)[,-2]
  nbloc<-dim(dats)[2]-1
  for (i in 1:nbloc){
    eff<-apply(table(dats[,(i+1)],dats[,1]),2,sum,na.rm=TRUE)
    freq<-sweep(table(dats[,(i+1)],dats[,1]),2,eff,FUN="/")
    if (i==1 ) all.freq<-freq else all.freq<-rbind(all.freq,freq)
  }
  return(t(all.freq))
}

##read the data (here, the name is for one of the simulated datasets used for this paper; fill in your own dataset, which should be in Fstat format)
#filename = sprintf("Simulated_Gradient_mig0.01_rep001.dat")
#dat = read.fstat.data(filename)

#coordinates of the populations, usually these are X and Y (or lat & long) of every population and read from a file
#here, since we have simulated data from a linear metapopulation, the coordinates are simply 1 to 20
#coords = cbind(X=1:20)

#calculate euclidean distance matrix

#extract pop-level coordiantes and Put back the real Chatham coordinates
variables_table2 <- read_excel("R PCA matrix2.xlsx")
variables_table_tar2<-variables_table2[-(c(129:143)),]
variables_table_tar2[,c(11,13)]->latlon2
latlon2$long_approx[latlon2$long_approx==182.9865]<--177.014
#move FRDL coord a bit to West so they are not considered enclave in land
latlon2$long_approx[latlon2$long_approx==166.9871]<-166.7
View(latlon2)

unique(latlon2)->latlon2

pC <- as.matrix(latlon2[c("long_approx","lat")])
cosDist2 <- costDistance(troceC, pC)
as.matrix(cosDist2)->m_cosDist2
#colnames(as.matrix(dnadist))->colnames(m_cosDist)
#rownames(as.matrix(dnadist))->rownames(m_cosDist)
#View(m_cosDist[1:100,1:100])



geomat = cosDist2

#necessary or you get an error: https://github.com/jgx65/hierfstat/issues/22
dat[] <- lapply(dat, type.convert)

#pairwise fst, may take a while
fstmat = array(0, c(14,14))
for(a in 2:14){
  for(b in 1:(a-1)){
    subdat = dat[dat[,1] == a | dat[,1] == b,]
    bs = basic.stats(subdat)
    fstmat[a,b] = fstmat[b,a] = bs$overall[8]
  }	
}
#do a mantel test
mant = mantel(fstmat,geomat)
fstlin = fstmat / (1-fstmat)
mantlin = mantel(fstlin,geomat)
print(mantlin)

#calculate overall fst
bstats = basic.stats(dat)
global_fst = bstats$overall[8]
print(bstats$overall)

#We need to keep only biallelic sites:

filter<-c()
for (i in 2:length(colnames(dat))) {
  if   
  (length(unique(dat[,i]))!=2) {
    filter<-c(filter,i)
  }}

dat_biallelic<-dat[,-filter]

#calculate pop freqs using custom function, not using hierfstat
freqs = pop.freqs(dat_biallelic)

#remove redundant frequencies; for every SNP lcous, only one of the allele frequencies is informative (as the other is its complement)
#don't use the code below for microsatellites or for SNPs with more than two alleles
allvar = apply(freqs,2,var)
freqs=freqs[,allvar>0]
freqs = freqs[,seq(1,dim(freqs)[2],2)]

#Put back the transformed coordinates for Chatham
latlon2$long_approx[latlon2$long_approx==-177.014]<-182.9865
#trick to convert to numeric
latlon2 <- as.matrix(latlon2[c("long_approx","lat")])


#make polynomials of spatial coordinates
space = poly(latlon2, degree=3)
#give correct names to polynomials depending on whether one or two dimensions are given in the coordinates (here only one is used)
if(dim(latlon2)[2]==1){
  colnames(space) <- c("X","X2","X3")
}else{
  colnames(space) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")	
}

#perform forward selection of spatial variables
#this may give a lot of output to the console
ord.spa = rda(freqs ~ ., data.frame(space), scale= FALSE)
(R2.all.spa = RsquareAdj(ord.spa)$adj.r.squared)
stp.spa = ordistep(rda(freqs ~ 1, data.frame(space)), scope = formula(ord.spa), scale= FALSE, direction="forward", pstep = 1000)
(selected.spa = attributes(stp.spa$terms)$term.labels)

#subsample space to the selected variables
#you may be interested to pause here to see which ones they are
space = space[,selected.spa]

perc.rda = 0
perc.tot = 0

#if there are any variables selected, perform an rda
if(length(selected.spa) > 0){
  #do the actual rda and study the output
  rda_res = rda(freqs, space,scale=FALSE) #nicer for plotys
  print(rda_res1)
  
  #calculate variance components
  tot_inert = rda_res$tot.chi
  orig_eigs = rda_res$CCA$eig
  
  all_perc = orig_eigs/orig_inert
  perc.rda = sum(orig_eigs) / tot_inert
  perc.tot = global_fst * perc.rda	
  
  #test overall significance and per axis
  rda_test = rda(freqs~space,scale=FALSE) #better for testing
  overall.test = anova.cca(rda_test, step=1000)
  axis.test = anova.cca(rda_test, by="axis", step= 1000)
  
  #set some plotting variables
  main =sprintf("Allele frequencies ~ Space (%.1f %%; Fst = %.2f; p = %s)", 100*perc.rda, perc.tot, overall.test[[5]][1]);
  xlab=sprintf("RDA1 (%.0f%%, p = %s)", 100*all_perc[1], axis.test[[5]][1])
  ylab=sprintf("RDA2 (%.0f%%, p = %s)", 100*all_perc[2], axis.test[[5]][2])
  lim = c(-1.2,1.2)
  
  #now draw plot
  plot(rda_res, scaling=1, main = main, xlab=xlab, ylab=ylab, display=c("sp","cn"), type="none", xlim=lim, ylim=lim)
  text(rda_res, scaling=1,display=c("cn"), col = "grey40", lwd=2, head.arrow = 0.12)
  points(rda_res, scaling=1, display = ("sp"), pch = 1, col="firebrick3")
  
  
}

