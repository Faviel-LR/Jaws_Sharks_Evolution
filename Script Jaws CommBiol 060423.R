library(Morpho) # Version 2.8
library(geomorph) # Version 4.0.0
library(phytools) # Version 0.7-80
library(geiger) # Version 2.0.7
library(viridisLite) # Version 0.4.0
library(vegan) # Version 2.5-7
library(ratematrix) # Version 1.2.3
library(mvMORPH) # Version 1.1.4
library(landvR) # Version 0.3
library(castor) # Version 1.6.9
library(dispRity) # Version 1.6.0
library(Hmisc) # Version 4.5-0
library(BTRTools) # Version  0.0.0.9000
library(coda) # Version 0.19-4
library(ggplot2) # Version 3.3.3
library(ggpubr) # Version 0.4.0
library(ggrepel) # Version 0.9.2
library(ggbreak) # Version 0.1.1
library(plotrix) # Version 3.8-2
library(Rphylopars) # Version 0.3.9
library(dendextend) # Version 1.16.0
library(heatmaply) # Version 1.4.0

########################################################################
###### Patch 
########################################################################

#### load the folder with the pts files of each specimen
setwd("/set/a/directory/pts")

ptslist<-dir(path=".",pattern='.pts',recursive=T)
ptsarray<-array(dim=c(196,3,145))
for(i in 1:length(ptslist)){
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],
                                      skip=2,header=F,sep="",row.names=1))
}

#### folder with the ply surface files
setwd("/set/a/directory/ply")

dimnames(ptsarray)[3]<-list(substr(dir(".",pattern=".ply"),1,
                                   (nchar(dir(".",pattern=".ply"))-4)))
arraylm<-ptsarray
JRsamp <- arraylm[-c(7,16,26,36,37,46,56,57,66,76,86,96,106,
                     116,126,127,136,146,147,156,166,167,176,
                     186,196),,]

#### semiautomatic surface landmarks placement, load the files
setwd("/set/a/directory")

template.lm<-read.table("template4.pts",skip=2,header=F,sep="")
TempRes.lm <- template.lm[-c(7,16,26,36,37,46,56,57,66,76,86,96,106,
                             116,126,127,136,146,147,156,166,167,176,
                             186,196),,]

patch2<-TempRes.lm[c(172:257),]

template.mesh<-file2mesh(filename = "./ply/CarcharhinusacronotusMCLMaraMotta.ply")
shade3d(template.mesh, col="seashell3") ### view template

fix<-as.integer(c(1:6))
slidings<-as.integer(c(7:171))
surface<-as.integer(c(172:257))

SC1<-as.integer(c(1,7:32,2))
SC2<-as.integer(c(2,33:49,3))
SC3<-as.integer(c(4,50:112,5))
SC4<-as.integer(c(3,113:128,4))
SC5<-as.integer(c(3,129:145,4))
SC6<-as.integer(c(2,146:171,6))
curv<-list(SC1,SC2,SC3,SC4,SC5,SC6)

### set a path to the folder with the ply files
meshlist<-dir("./ply",pattern=".ply")

atlas<-createAtlas(template.mesh,as.matrix(TempRes.lm[c(1:171),]),
                   patch = as.matrix(patch2),
                   corrCurves= list(c(7:32),c(33:49),c(50:112),c(113:128),
                                    c(129:145),c(146:171)),
                   patchCurves = NULL, keep.fix = fix)

#### supplenmentary figure 1
plotAtlas(atlas, meshcol = "white", render = "s")


patched1<-placePatch(atlas, JRsamp,
                     path=".", prefix =NULL,
                     fileext = ".ply", ray = TRUE,
                     inflate = 3.5, tol = 5,
                     relax.patch = FALSE, keep.fix = fix,
                     silent = FALSE)
checkLM(patched1,path="./",suffix=".ply",atlas=atlas,point = "s", 
        pt.size = 0.5,render="s",alpha=1)


relaxedtemplate1<-array(dim=c(257,3,145))
for (i in 1:dim(patched1)[3])
{
  print(i) 
  relaxedtemplate1[,,i]<-relaxLM(patched1[,,i], rbind(atlas[["landmarks"]],atlas[["patch"]]),
                                 slidings, outlines = curv, surp = surface,
                                 sur.name = NULL, mesh = file2mesh(paste("./",meshlist[i],
                                                                         sep="")),
                                 tol = 1e-50, deselect = FALSE, inc.check = F, 
                                 iterations = 100, fixRepro = TRUE, stepsize=0.2)
}

dimnames(relaxedtemplate1)[3]<-dimnames(patched1)[3]

#### save the result
save(relaxedtemplate1,file="./Resampled Raw Jaws.R")
load(file="./Resampled Raw Jaws.R")

########################################################################
###### GPA
########################################################################

setwd("/set/your/working/directory")

### load this files from your working directory
RJL <- readland.tps("RawJaws.TPS",specID = "ID",warnmsg = TRUE)
Jcurves<-as.matrix(read.csv("SlidersRes.csv",header = T, sep = ","))
JSurf<-as.matrix(read.table("Surf-Jaws2.txt",header = F))

### GPA
jaws3d.gpa<-gpagen(RJL,curves = Jcurves,surfaces = JSurf,
                   ProcD = F)

########################################################################
###### LASEC Function from Watanabe, how many landmarks are needed
########################################################################
lasec <- function(coord.data, n.dim, iter=1000, show.progress=T) {
  
  # coord.data is a n x p matrix with n=specimens, p=shape variables (could be
  # unaligned).
  # n.dim is spatial dimensionality of coord.data (i.e., 2 or 3 for 2-D, 3-D)
  # show.progress is for displaying the nth iteration on console.
  
  # LIBRARY #
  require(geomorph)
  require(vegan)
  
  # INITIALIZATION #
  matrix.pss.cs <- NULL
  matrix.pss <- NULL
  matrix.chosen <- NULL
  lm.maxpss <- NULL
  lm.minpss <- NULL
  output <- NULL
  if(show.progress==T) {
    pb <- txtProgressBar(min=0, max=iter, char=" >", style=3)
  }
  
  # DATA #
  n.lm <- ncol(coord.data)/n.dim   # no. landmarks in data
  list.lm <- seq(1:n.lm)   # make a list of landmarks
  coord.data <- arrayspecs(coord.data, n.lm, n.dim)
  gpa <- gpagen(coord.data, print.progress=FALSE)
  shape.data <- two.d.array(gpa$coords)
  cs.data <- gpa$Csize
  coord.data <- two.d.array(coord.data)
  
  for(i in 1:iter) {
    if(show.progress==T) {
      setTxtProgressBar(pb, i)
    }
    new.lm.order <- sample(n.lm, replace=F)   # new sequence of landmarks
    list.pss <- NULL
    list.pss.cs <- NULL
    for(j in 3:n.lm) {
      subsample.lm <- new.lm.order[1:j]
      if(n.dim==3) {
        chosen.nvar <- sort(c(3*subsample.lm-2, 3*subsample.lm-1, 3*subsample.lm))
      }
      if(n.dim==2) {
        chosen.nvar <- sort(c(2*subsample.lm-1, 2*subsample.lm))
      }
      sampled.coord <- coord.data[,chosen.nvar]
      sampled.coord <- arrayspecs(sampled.coord, j, n.dim)
      sampled.gpa <- gpagen(sampled.coord, print.progress=FALSE)
      sampled.cs <- sampled.gpa$Csize
      sampled.shape <- two.d.array(sampled.gpa$coords)
      
      # CALCULATE 'FIT' BETWEEN FULL & SUBSAMPLED DATA #
      pss <- protest(shape.data, sampled.shape, permutations=0)$ss
      pss.cs <- protest(cs.data, sampled.cs, permutations=0)$ss
      
      # RECORD 'FIT' BETWEEN FULL & SUBSAMPLED DATA #
      list.pss <- append(list.pss, pss)
      list.pss.cs <- append(list.pss.cs, pss.cs)
    }
    matrix.chosen <- append(matrix.chosen, new.lm.order)
    matrix.pss <- append(matrix.pss, c(NA, NA, list.pss))
    matrix.pss.cs <- append(matrix.pss.cs, c(NA, NA, list.pss.cs))
  }
  
  matrix.chosen <- matrix(matrix.chosen, nrow=iter, byrow=T)
  matrix.pss <- 1-matrix(matrix.pss, nrow=iter, byrow=T)
  median.pss <- apply(matrix.pss, 2, median)
  matrix.pss.cs <- 1-matrix(matrix.pss.cs, nrow=iter, byrow=T)
  median.pss.cs <- apply(matrix.pss.cs, 2, median)
  
  # PLOT PROC SS CURVE #
  pdf("LaSEC_SamplingCurve_Shape.pdf")
  plot(x=NA, xlim=c(0, n.lm), ylim=c(0,1), xlab="No. landmarks", ylab="Fit")
  for(i in 1:iter) {
    par(new=T)
    plot(matrix.pss[i,], xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="grey", xlab="", ylab="", axes=F)
  }
  par(new=T)
  plot(median.pss, xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="black", lwd=3, xlab="", ylab="", axes=F)
  dev.off()
  
  # PLOT PROC SS CENTROID SIZE CURVE #
  pdf("LaSEC_SamplingCurve_CS.pdf")
  plot(x=NA, xlim=c(0, n.lm), ylim=c(0,1), xlab="No. landmarks", ylab="Fit")
  for(i in 1:iter) {
    par(new=T)
    plot(matrix.pss.cs[i,], xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="grey", xlab="", ylab="", axes=F)
  }
  par(new=T)
  plot(median.pss.cs, xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="black", lwd=3, xlab="", ylab="", axes=F)
  dev.off()
  
  # OUTPUT #
  output$fit <- matrix.pss
  output$median.fit <- median.pss
  output$maxfit.landmark <- table(lm.maxpss)
  output$minfit.landmark <- table(lm.minpss)
  output$fit.cs <- median.pss.cs
  return(output)
}


J3D.smth <- two.d.array(jaws3d.gpa$coords)

# lasec(coord.data=example.txt, n.dim=3, iter=1000)

lasec(J3D.smth,3,iter=1000,show.progress = TRUE)
### supplementary figure 2 obtained after this function is applied

########################################################################
###### Resample number of landmarks after LASEC function, and GPA
########################################################################

RJL <- readland.tps("RawJaws.TPS",specID = "ID",warnmsg = TRUE)

Jaws3d.lm <- RJL[c(1:6,9,12,15,18,21,24,27,30,35,38,41,44,47,53,57,
                   61,65,69,73,77,81,85,89,93,97,101,105,109,114,
                   116,118,121,123,125,127,130,132,134,136,138,140,
                   142,144,148,151,154,157,160,163,167,170,173,175,
                   177,179,181,183,185,187,188,190,193,195,196,198,
                   201,203,204,206,208,210,212,214,216,218,221,223,
                   225,226,228,230,233,235,237,238,240,242,244,246,
                   248,250,252,254,256),,]

### use these sliders and surface files for the new configuration
Jcurves<-as.matrix(read.csv("SlidersJ3D.csv",header = T, sep = ","))
JSurf<-as.matrix(read.table("Surf-J3D.txt",header = F))

jaws3d.gpa<-gpagen(Jaws3d.lm,curves = Jcurves,surfaces = JSurf,
                   ProcD = F)

save(jaws3d.gpa,file="./GPA3DJ.R")
load(file="./GPA3DJ.R")

########################################################################
###### Cluster of stomach content by species
########################################################################
DietJaws<-read.table(as.matrix("Supplementary table 2.csv"),row.names = 1,
                     header=T, sep = ",")

### Subset without outliers to check clustering differences
DietJawsNOUT <- DietJaws[-c(4,20,33,37,39,54,73),]  

### substitute the elements accordingly
DietJ.dist<-vegdist(DietJaws, method = "bray")
DietJ.clust <- hclust(DietJ.dist,  method = "average")
plot(DietJ.clust, horiz =T)
Diet.Den<-as.dendrogram(DietJ.clust)
par(cex=0.5,  mar=c(3, 1, 1, 9))
plot(Diet.Den, horiz=T)

den_k<-find_k(DietJ.clust)
plot(den_k)
plot(color_branches(DietJ.clust, k = den_k$nc))
plot(color_branches(DietJ.clust, k = 8))

#### supplementary figure 3 & 4
heatmaply(as.matrix(DietJaws),
          colors = viridis(n=256, 
                           option = "inferno"),
          k_col = 0, k_row = 8, fontsize_row = 5, dendrogram = "row",
          distfun = function(x) vegdist(x,method = "bray"),
          hclust_method = "average")


########################################################################
###### Read Classifiers table
########################################################################

class.Raw <- read.csv("Supplementary table 1.csv",header = TRUE, sep = ",",
                      stringsAsFactors = T)
gp2 <- as.factor(paste(class.Raw$Order))
col2.gp <- c("darkseagreen4","black","cyan3","darkorchid",
             "gold","darkorange","lightsteelblue",
             "deepskyblue4","firebrick3","navy")
names(col2.gp) <- levels(gp2)
col2.gp <- col2.gp[match(gp2,names(col2.gp))]

########################################################################
###### Geomorph
########################################################################
########################################################################
###### PCA All specimens
########################################################################

Jaws.pc.all<-gm.prcomp(jaws3d.gpa$coords)
plot(Jaws.pc.all, axis1 = 1, axis2 = 2, pch=21,
     cex = 1.5, bg=col2.gp) ### figure 1

summary(Jaws.pc.all)
getMeaningfulPCs(Jaws.pc.all$sdev^2, n = 144, sdev = FALSE)
#### 4 components are the selected ones to describe
######## scree plot
varExp = (100*Jaws.pc.all$sdev^2)/sum(Jaws.pc.all$sdev^2)
varDF = data.frame(Dimensions=1:length(varExp),
                   varExp=varExp)

#### supplementary figure 8
ggplot(varDF,aes(x=Dimensions,y=varExp)) + geom_point() + 
  geom_col(fill="steelblue") + geom_line() + 
  theme_bw() + scale_x_continuous(breaks=1:nrow(varDF)) + 
  ylim(c(0,50)) +
  xlim(c(0,50)) +
  ylab("% variance explained")

### Optional Add labels on each point 
df <- data.frame(Sp =  class.Raw$Species,
                 Ord = class.Raw$Order,
                 Sp2 = class.Raw$Species2,
                 Fam = class.Raw$Fam2,
                 PC1 = Jaws.pc.all$x[,1],
                 PC2 = Jaws.pc.all$x[,2],
                 PC3 = Jaws.pc.all$x[,3],
                 PC4 = Jaws.pc.all$x[,4])

#### supplementary figure S6 & S7
ggplot(df, aes(PC1,PC2, color = Ord)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_text_repel(mapping =  aes(label = Sp2, fontface = "italic"),
                  colour = "black", size = 2)+
  scale_color_manual(values=c("darkseagreen4","black","cyan3","darkorchid",
                              "gold","darkorange","lightsteelblue",
                              "deepskyblue4","firebrick3","navy"))+
  xlab("PC1 (40.89%)") + ylab("PC2 (15.74%)")+
  theme_minimal()

##########################################################################
#### Display landmark variation part of figure 1 and Suppl. S6 & S7
##########################################################################

procrustes_2d_array <- geomorph::two.d.array(jaws3d.gpa$coords)
ordination <- stats::prcomp(procrustes_2d_array)

#### extract the values for the corresponding PC, example with PC1
PC1_variation <- variation.range(jaws3d.gpa, axis = 1,
                                 ordination = ordination, type = "spherical")

hyp_P1min <- as.matrix(Jaws.pc.all$shapes$shapes.comp1$min[1:100,],
                       nrow = 100, ncol = 3)
hyp_P1max <- as.matrix(Jaws.pc.all$shapes$shapes.comp1$max[1:100,],
                       nrow = 100, ncol = 3)
procrustes.var.plot(hyp_P1min, hyp_P1min, col = heat.colors, 
                    magnitude = 1, pt.size = 2,
                    col.val = PC1_variation[, "radius"], labels = FALSE)

########################################################################
###### Average coordinates by species
########################################################################

Jaws.res<-two.d.array(jaws3d.gpa$coords)
means<-(aggregate(Jaws.res~class.Raw$TipLab, FUN = mean))[,-1]
rownames(means)<-levels(class.Raw$TipLab)
ElasmoG<-arrayspecs(means,100,3)

########################################################################
###### Read Tree(s)
########################################################################

set.seed(183)
Phy<-read.tree("JTree.tre") ### MCC tree
Phys<-read.nexus("output.nex") ### 1000 trees

########################################################################
##### Select a subset of trees to account for phylogenetic uncertainty
########################################################################

random.trees<-sample(Phys,size=100)

#### only sharks
Sel.trees <- drop.tip.multiPhylo(random.trees,
                                 tip = "Callorhinchus_milii")

#### only orders with more than three species
Ord.trees <- drop.tip.multiPhylo(random.trees,
                                 tip = c("Callorhinchus_milii","Echinorhinus_brucus",
                                         "Heterodontus_francisci","Heterodontus_japonicus",
                                         "Pristiophorus_japonicus"))

#### only diet groups with more than three species
Die.trees <- drop.tip.multiPhylo(random.trees, 
                                 tip = c("Callorhinchus_milii","Apristurus_laurussonii",
                                         "Atelomycterus_macleayi","Carcharhinus_hemiodon",
                                         "Centrophorus_seychellorum","Cetorhinus_maximus",
                                         "Chiloscyllium_arabicum","Eucrossorhinus_dasypogon",
                                         "Halaelurus_boesemani","Hemiscyllium_trispeculare",
                                         "Heterodontus_francisci","Heterodontus_japonicus",
                                         "Isogomphodon_oxyrhynchus","Isurus_paucus",
                                         "Orectolobus_japonicus","Oxynotus_centrina",
                                         "Pseudocarcharias_kamoharai","Scyliorhinus_boa",
                                         "Squatina_japonica","Stegostoma_fasciatum"))

#### All diet groups
Die2.trees <- drop.tip.multiPhylo(random.trees, 
                                  tip = c("Apristurus_laurussonii","Atelomycterus_macleayi",
                                          "Callorhinchus_milii","Carcharhinus_hemiodon",
                                          "Centrophorus_seychellorum","Chiloscyllium_arabicum",    
                                          "Eucrossorhinus_dasypogon","Halaelurus_boesemani",
                                          "Heterodontus_japonicus","Isogomphodon_oxyrhynchus",
                                          "Isurus_paucus","Orectolobus_japonicus",
                                          "Pseudocarcharias_kamoharai","Scyliorhinus_boa",
                                          "Squatina_japonica"))

########################################################################
###### Build Subsets for the trimmed trees
########################################################################

class.Tre <- read.csv("Supplementary table 3.csv",header = TRUE,
                      sep = ",",row.names = 1, stringsAsFactors = TRUE)

Sub.sel <- subset(class.Tre, GP1 == "Sel")
Sub.sel <- droplevels(Sub.sel)

Sub.die <- subset(class.Tre, GP3 == "SC")
Sub.die <- droplevels(Sub.die)

Sub.die2 <- subset(class.Tre, GP4 == "SC")
Sub.die2 <- droplevels(Sub.die2)

Sub.ord <- subset(class.Tre, GP2 == "Mord")
Sub.ord <- droplevels(Sub.ord)

######################
Hab<-as.factor(Sub.sel$Hab)
names(Hab)<-rownames(Sub.sel)

TrL<-as.factor(Sub.sel$TL_C)
names(TrL)<-rownames(Sub.sel)

Ord<-as.factor(Sub.ord$Order)
names(Ord)<-rownames(Sub.ord)

Diet<-as.factor(Sub.die2$Diet)
names(Diet)<-rownames(Sub.die2)

######################
### subset the MCC tree

species<-Sub.sel$TipLab
Phy.sel<-drop.tip(Phy,Phy$tip.label
                  [-match(species, Phy$tip.label)])

species<-Sub.ord$TipLab
Phy.ord<-drop.tip(Phy,Phy$tip.label
                  [-match(species, Phy$tip.label)])

species<-Sub.die2$TipLab
Phy.dd3<-drop.tip(Phy,Phy$tip.label
                  [-match(species, Phy$tip.label)])

#####################
#### Subset landmark data for mvMorph

Array2D<-two.d.array(ElasmoG[,,c(1:90)])
Jaws.land <- Array2D[Phy$tip.label,]

Jaws.sel <- Jaws.land[-90,] ## only sharks
Jaws.ord <- Jaws.land[-c(66,67,82,86,90),] ## orders with more than 3 species
Jaws.die <- Jaws.land[-c(12,46,90,26,69,10,57,53,13,62,
                         66,67,25,7,54,72,3,49,85,65),] ## No monospecific diet
Jaws.die2 <- Jaws.land[-c(3,7,12,13,25,26,46,49,53,54,57,67,69,85,90),] ## all diet groups

########################################################################
### color vectors for specific subset
########################################################################

n<-4
p4 <- inferno(n)
names(p4) <- levels(Sub.sel$Hab)
p4.gp <- p4[match(Sub.sel$Hab,names(p4))]

n<-3
p3 <- viridis(n)
names(p3) <- levels(Sub.sel$TL_C)
p3.gp <- p3[match(Sub.sel$TL_C,names(p3))]

die.col<-c("#eeeb0f","#249c20","#2e2067","#ab1aab","#3597df",
           "#de5f1b","#dc9d1f","#b28130")
names(die.col) <- levels(Sub.die2$Diet)
colDie.gp <- die.col[match(Sub.die2$Diet,names(die.col))]

gpOrd <- as.factor(paste(class.Tre$Order))
col2.gp <- c("darkseagreen4","black","cyan3","darkorchid",
             "gold","darkorange","lightsteelblue",
             "deepskyblue4","firebrick3","navy")
names(col2.gp) <- levels(gpOrd)
col2.gp <- col2.gp[match(gpOrd,names(col2.gp))]

########################################################################
### Phylogenetically aligned PCA
########################################################################

### all species
All.paca <- gm.prcomp(ElasmoG, phy = Phy,
                      align.to.phy = TRUE)
summary(All.paca)

#### attention here
getMeaningfulPCs(All.paca$sdev, n = 89, sdev = FALSE)

### Only sharks
Sel.land<-ElasmoG[,,-7]
Sel.paca <- gm.prcomp(Sel.land, phy = Phy.sel,
                      align.to.phy = TRUE)

PhySig<-physignal(Sel.land, Phy.sel, iter = 1000)
summary(PhySig)

### substitute arguments accordingly after observing the variation explained
### by each PaCA
#### figure 2 a-d substitute with the corresponding bg value
plot(All.paca, axis1 = 1, axis2 = 2, pch=21, cex = 1.5,
     xlab = "PaCA 1: 70.88 %", ylab = "PaCA 23.76 %",
     bg=col2.gp, phylo = TRUE, time.plot = FALSE,
     phylo.par = list(edge.color="grey50", anc.states =FALSE,
                      node.cex=0, node.labels = FALSE,
                      tip.labels = FALSE,
                      tip.txt.cex = 0.5))

########################################################################
### dispRity 
########################################################################
### example with orders, select from the data the corresponding groups

dispJaws2<- dispRity.per.group(Jaws.sel,
                               list(Hexa=c(87:89),
                                    Squa=c(68:81),
                                    Squt=c(83:85),
                                    Orec=c(52:65),
                                    Lamn=c(1:11),
                                    Carc=c(12:51)),
                               metric=c(sum, variances))
summary(dispJaws2)

test.dispRity(dispJaws2, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(dispJaws2, ylab = "Procrustes Variance")

dispHab2 <- dispRity.per.group(as.matrix(J.dis.sel),
                               list(Deep=c(3,19:22,31:35,37,45,48,49,55,65,71,81:83,85),
                                    Pelagic=c(1,2,9,11,18,24,39,51:54,66,68,69,75,76,79),
                                    Shelf=c(13:15,23,25,38,40,42,43,50,56:59,67,70,72,74,
                                            78,80,84,86,89),
                                    Reef=c(4:8,10,12,16,17,26:30,36,41,44,46,47,60:64,73,
                                           77,87,88)),
                               metric=c(sum, variances))
summary(dispHab2)

dispTL <- dispRity.per.group(as.matrix(J.dis.sel),
                             list(LP=c(24,27,28,44,46,47,56:59,62,65,67,70,72,78,85,87,
                                       89),
                                  MP=c(3:6,10,12,13,15,16,19,20,22,23,25,26,29:32,35:38,
                                       41,43,45,55,60,61,63,64,68,71,73:75,77,81,84,86,88),
                                  TP=c(1,2,7:9,11,14,17,18,21,33,34,39,40,42,48:54,66,69,
                                       76,79,80,82,83)),
                             metric=c(sum, variances))

summary(dispTL)


dispDie <- dispRity.per.group(as.matrix(J.dis.sel),
                              list(BP=c(8,14,17,18,39,76),
                                   CEPH=c(1,2,11,21,22,29,31,35,37,43,45,48,49,
                                          60,66,73,74,79),
                                   CR=c(28,34,41,56:59,67,70,77,78,83),
                                   FISH=c(6,7,9,10,12,15,16,19,25,32,33,38,40,51,
                                          53:55,61,63,64,69,75,81,82,84,86,88),
                                   GEN=c(5,23,27,30,72,80,89)),
                              metric=c(sum, variances))

summary(dispDie)

BP<-as.matrix(t(dispDie$disparity$BP[[2]]))
CEPH<-as.matrix(t(dispDie$disparity$CEPH[[2]]))
CR<-as.matrix(t(dispDie$disparity$CR[[2]]))
FISH<-as.matrix(t(dispDie$disparity$FISH[[2]]))
GEN<-as.matrix(t(dispDie$disparity$GEN[[2]]))

Carc<-as.matrix(t(dispJaws2$disparity$Carc[[2]]))
Hexa<-as.matrix(t(dispJaws2$disparity$Hexa[[2]]))
Lamn<-as.matrix(t(dispJaws2$disparity$Lamn[[2]]))
Orec<-as.matrix(t(dispJaws2$disparity$Orec[[2]]))
Squa<-as.matrix(t(dispJaws2$disparity$Squa[[2]]))
Squt<-as.matrix(t(dispJaws2$disparity$Squt[[2]]))

Deep<-as.matrix(t(dispHab2$disparity$Deep[[2]]))
Pelagic<-as.matrix(t(dispHab2$disparity$Pelagic[[2]]))
Shelf<-as.matrix(t(dispHab2$disparity$Shelf[[2]]))
Reef<-as.matrix(t(dispHab2$disparity$Reef[[2]]))

LP<-as.matrix(t(dispTL$disparity$LP[[2]]))
MP<-as.matrix(t(dispTL$disparity$MP[[2]]))
TP<-as.matrix(t(dispTL$disparity$TP[[2]]))

DispAll <- as.data.frame(cbind(Carc,Lamn,Orec,Squa,Squt,Hexa,BP,
                               CEPH,CR,FISH,GEN,LP,MP,TP,Deep,Pelagic,Shelf,Reef))
colnames(DispAll) <- c("Carc","Lamn","Orec","Squa","Squt","Hexa","BP","CEPH","CR",
                       "FISH","GEN","LP","MP","TP","Deep","Pelagic","Shelf","Reef")

#### Figure 3b
ggplot(stack(DispAll), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1, alpha=0.8,
              aes(colour = factor(ind))) + 
  scale_colour_manual(values = c("darkseagreen4","darkorange","lightsteelblue","firebrick3","navy","gold",
                                 "darkslateblue","forestgreen","yellow2","darkorchid", "deepskyblue2",
                                 "#440154FF", "#21908CFF", "#FDE725FF",
                                 "#000004FF", "#781C6DFF","#A8AA6D","#ED6925FF")) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  scale_x_discrete(limits=c("Hexa","Squa","Squt","Orec","Lamn","Carc",
                            "BP","CEPH","CR","FISH","GEN",
                            "LP","MP","TP",
                            "Deep","Pelagic","Shelf","Reef"))+
  theme_minimal()+
  ylab("Procrustes Variance")+
  stat_summary(fun.y=mean, geom="point", shape=21, 
               size=5, color="black", fill="white")+
  theme(legend.position = "none", axis.title.x = element_blank())


########################################################################
### Disparity through time
########################################################################

#### with MCC Tree
dispRity_dtt <- dtt.dispRity(data = Jaws.sel, metric = c(sum, variances),
                             tree = Phy.sel, nsim = 100)

plot(dispRity_dtt)

#### with the subset of trees
#### make two matrices to store the values
DTT_Times <- matrix(NA, ncol = 89, nrow=length(Sel.trees))
DTT_DTT <- matrix(NA, ncol = 89, nrow=length(Sel.trees))

for(i in 1:length(Sel.trees)) {
  tree <- Sel.trees[[i]]
  dttSel  <- dtt.dispRity(data = Jaws.sel, metric = c(sum, variances),
                          tree = tree, nsim = 1)
  DTT_Times[i, 1:89] <- dttSel$times
  DTT_DTT[i, 1:89]<- dttSel$dtt
}

dtt.df <- as.data.frame(DTT_Times)
dta.df <- as.data.frame(DTT_DTT)
dtt.sin <- data.frame(y=unlist(dtt.df))
dta.sin <- data.frame(y=unlist(dta.df))

x1<-rep(c("T01","T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10",
          "T11","T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20",
          "T21","T22", "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
          "T31","T32", "T33", "T34", "T35", "T36", "T37", "T38", "T39", "T40",
          "T41","T42", "T43", "T44", "T45", "T46", "T47", "T48", "T49", "T50",
          "T51","T52", "T53", "T54", "T55", "T56", "T57", "T58", "T59", "T60",
          "T61","T62", "T63", "T64", "T65", "T66", "T67", "T68", "T69", "T70",
          "T71","T72", "T73", "T74", "T75", "T76", "T77", "T78", "T79", "T80",
          "T81","T82", "T83", "T84", "T85", "T86", "T87", "T88", "T89", "T90",
          "T91","T92", "T93", "T94", "T95", "T96", "T97", "T98", "T99", "T100"),
        times=89)

dtt.com <- cbind(dtt.sin,dta.sin,x1)
colnames(dtt.com) <- c("Time","DTT","Tree")

#### figure 5c
ggplot(dtt.com, aes(x= Time, y = DTT )) + 
  geom_line(aes(group=Tree), color = "darkred") + 
  geom_point(alpha=.3, size = 0.5) +
  theme_bw()

####################################################################
################ Discrete traits evolutionary model, AIC
### from the ACE function in: Brazeau, M.D., Giles, S., Dearden, R.P. et al. 
### Endochondral bone in an Early Devonian ‘placoderm’ from Mongolia.
### Nat Ecol Evol 4, 1477–1484 (2020).
### https://doi.org/10.1038/s41559-020-01290-2
### Example with Habitat and ER model, modify accordingly for each model and factor

Hab_L_ER <- matrix(NA, ncol= 2, nrow=length(Sel.trees))

for(i in 1:length(Sel.trees)) {
  tree <- Sel.trees[[i]]
  ## Order the character to match the one from the input tree
  character <- Sub.sel[Sel.trees[[i]]$tip.label,]
  ## Select the trait from the data frame column
  tokens <- as.numeric(character[,10]) 
  
  fitER  <- fit_mk(trees = tree, Nstates = 4,
                   tip_states = tokens, 
                   rate_model="ER", Ntrials=2, 
                   Nthreads = 4, 
                   Nbootstraps = 100,
                   verbose = TRUE)
  Hab_L_ER[i,1] <- fitER$AIC
  Hab_L_ER[i,2] <- fitER$loglikelihood
}

### do the same with the remaining models and assign column names
colnames(Hab_L_ER)  <- c("AIC_ER_Hab" , "logLik_ER_Hab")
colnames(Hab_L_SYM) <- c("AIC_SYM_Hab", "logLik_SYM_Hab")
colnames(Hab_L_ARD) <- c("AIC_ARD_Hab", "logLik_ARD_Hab")

### concatenate the results from the models
Hab_Mod <- cbind(Hab_L_SYM, Hab_L_ARD, Hab_L_ER)
apply(Hab_Mod, 2, mean)
boxplot(Hab_Mod[,c(1,3,5)])

###############################################################
####### Prepare landmark data for mvMorph
###############################################################

Y <- Jaws.sel
data <- list(Y=Jaws.sel)
Z <- Jaws.ord
dataZ <- list(Z=Jaws.ord)
D <- Jaws.die2
dataD <- list(D=Jaws.die2)

#######
### list for Sel landmark subset
data.sel=list(shape=Jaws.sel,HabJ=Hab,TrLJ=TrL,CS=log(Sub.sel$CS))
### list for Ord landmark subset
data.ord=list(shape=Jaws.ord,OrdJ=Ord,log(CS=Sub.ord$CS))
### list for Diet landmark subset
data.die=list(shape=Jaws.die,DieJ=Diet,CS=log(Sub.die$CS))
### list for Diet landmark subset (75spp)
data.die2=list(shape=Jaws.die2,DieJ=Diet,CS=log(Sub.die2$CS))

#### make simm.map with subset of trees for each corresponding trait

### figure 4a & 4b
my_treesSYM_Hab<-make.simmap(Sel.trees,Hab,model="SYM")

#### 
save(my_treesSYM_Hab,file="./Test_trees_Script.R")
load(file="./Test_trees_Script.R")

###############################################################
##### test model for continuous trait, example with BM, modify accordingly 

simBM_Coords <- matrix(NA, ncol= 3, nrow=length(Sel.trees))

simBM <- for(i in 1:length(Sel.trees)) {
  Y <- Y[Sel.trees[[i]]$tip.label,]
  
  fitBM <- mvgls(Y~1,
                 tree = Sel.trees[[i]], penalty = "RidgeArch",
                 model = "BM", method = "H&L",
                 error = TRUE, nbcores = 4L)
  simBMgic <- GIC(fitBM)
  simBM_Coords[i,1] <- simBMgic$LogLikelihood
  simBM_Coords[i,2] <- simBMgic$GIC
}

colnames(simBM_Coords) <- c("logLik_BM","GIC_BM","bias_BM")
colnames(simEB_Coords) <- c("logLik_EB","GIC_EB","bias_EB")
colnames(simOU_Coords) <- c("logLik_OU","GIC_OU","bias_OU")
colnames(simLA_Coords) <- c("logLik_LA","GIC_LA","bias_LA")

ECoord_Mod <- cbind(simBMR, simEB_Coords, simOU_Coords, simLA_Coords)

apply(ECoord_Mod, 2, mean)

fit1  <- mvgls(Y~1, tree = Phy.sel, model="BM", penalty="RidgeArch",
               method = "PL-LOOCV", nbcores = 4L)
fit2  <- mvgls(Y~1, tree = Phy.sel, model="OU", penalty="RidgeArch",
               method = "PL-LOOCV", nbcores = 4L)
fit3  <- mvgls(Y~1, tree = Phy.sel, model="EB", penalty="RidgeArch",
               method = "PL-LOOCV" , nbcores = 4L)

GIC(fitBM); GIC(fitEB); GIC(fitOU)

signal <- mvgls(Y~1, tree = Phy.sel, model="lambda", penalty="RidgeArch",
                method = "H&L", nbcores = 4L)
GIC(signal)

###############################################################
#### evolutionary rates mvMorph, example with Habitat

simHabT <- sapply(1:length(my_treesSYM_Hab), function(x) {
  fitBM <- mvgls(Y~1, data=data,
                 tree = my_treesSYM_Hab[[x]],
                 model = "BMM", method = "H&L",
                 error = TRUE, nbcores = 4L)
  fitBM$param
})

rowMeans(simHabT)
boxplot(t(simHabT))

### save the result
save(simHabT,file="./choose_a_name.R")
load(file="./choose_a_name.R")

load(file="./Hab_EvR100921.R")
load(file="./TrL_EvR100921.R")
load(file="./Die70_EvR100921.R")
load(file="./Ord_EvR100921.R")

ER_all <- as.data.frame(cbind(Ord_EvR,Die70_EvR,TrL_EvR,Hab_EvR))

#### Figure 3a
ggplot(stack(ER_all), aes(x = ind, y = values))+
  geom_jitter(position=position_jitter(0.2), size=1, alpha=0.8,
              aes(colour = factor(ind))) + 
  scale_colour_manual(values = c("darkseagreen4","gold","darkorange","lightsteelblue","firebrick3","navy",
                                 "darkslateblue","forestgreen","yellow2","darkorchid", "deepskyblue2",
                                 "#440154FF", "#21908CFF", "#FDE725FF",
                                 "#000004FF", "#781C6DFF","#ED6925FF","#A8AA6D")) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3)+
  theme_minimal()+
  ylab("Evolutionary Rates")+
  stat_summary(fun.y=mean, geom="point", shape=21, 
               size=5, color="black", fill="white")+
  scale_y_break(c(3.5e-06,1.03e-05))+
  scale_x_discrete(limits=c("Hexa","Squa","Squt","Orec","Lamn","Carc",
                            "BP","CEPH","CR","FISH","GEN",
                            "LP","MP","TP",
                            "deep","pelagic","shelf","reef"))+
  theme(legend.position = "none", axis.title.x = element_blank())

###############################################################
###testing for jaw shape differences, example with habitat
###############################################################

fit.Hab<-mvgls(shape ~ CS*Hab, data=data.sel,
               tree = Phy.sel, model="lambda", method="H&L",
               penalty = "RidgeArch",
               verbose = TRUE,
               nbcores = 4L)

summary(fit.Hab)

###############################################################
#### Performing multivariate tests on generalized least squares
#### linear model, example with Habitat

mult_hab <- manova.gls(fit.Hab, nperm=999, type="II",
                       test="Pillai", verbose = TRUE,
                       nbcores = 4L)
mult_hab

###############################################################
### Evolution rate single landmark
### load coordinates and mesh to visualize, example with Lamniformes
###############################################################

land.gp<-as.factor(c(1:100))
### select species from the GPA 
Lamn<-ElasmoG[,,c(1,2,18,19,25,52:56,69)]
### read the raw landmarks
lamn.lm<-read.table("lamnlm.pts",skip=2,header=F,sep="")
### resample to 100 coordinates
Lamn.lm <- lamn.lm[-c(7,16,26,36,37,46,56,57,66,76,86,96,106,
                      116,126,127,136,146,147,156,166,167,176,
                      186,196),,]
lamn.mesh<-file2mesh(filename = "./LamnaditropisMCLERB0937.ply")

### select species
Lamn.sp<-c("Alopias_superciliosus","Alopias_vulpinus","Carcharias_taurus",
           "Carcharodon_carcharias","Cetorhinus_maximus","Isurus_oxyrinchus",
           "Isurus_paucus","Lamna_ditropis","Lamna_nasus",
           "Mitsukurina_owstoni","Pseudocarcharias_kamoharai")

### trim the tree
Phy.lam <- drop.tip(Phy,Phy$tip.label
                    [-match(Lamn.sp, Phy$tip.label)])
plotTree(Phy.lam,fsize=0.5)

### mean rates
EMR<-compare.multi.evol.rates(A=Lamn,gp=land.gp, 
                              Subset=FALSE, 
                              phy= Phy.lam,iter=999)
LM.emr <- as.matrix(EMR$sigma.d.gp)
LM.sca <- as.factor(as.numeric(cut2(LM.emr, g = 20)))

n<-20
p4 <- viridis(n)
names(p4) <- levels(as.factor(LM.sca))
col.land <- p4[match(LM.sca, names(p4))]
shade3d(lamn.mesh, col="seashell3")
#### supplementary figure S9
spheres3d(Lamn.lm, color = col.land, specular = 1, radius = 0.7)

###############################################################
############ Branch-Specific rate, BayesTraites ################
#### details in: https://github.com/anjgoswami/mammal_jaws/blob/main/Script_jaw_mammals_paper.r
#loading packages
library(BTRTools)
library(phytools)
library(coda)

tracePlots <- function(file, burnin=0, thinning=1, plot=TRUE, 
                       display=c("Lh", "Lh...Prior", "No.Pram", 
                                 "Alpha", "Sigma.2")){
  require(BTRTools)
  require(coda)
  
  rjout <- loadRJ(file, burnin = burnin, thinning = thinning)
  chain_out <- type.convert(rjout$rj_output)
  rownames(chain_out) = chain_out[,"It"]
  chain_out = chain_out[,-1]
  # Retrieve numerical
  index <- sapply(chain_out,function(x) is.numeric(x))
  chain <- mcmc(chain_out[,index])
  
  # plot the trace
  if(plot){
    plot(chain[,display])
  }
  
  # Just compute some statistics (autocorrelation...)
  cat("Effective sample size:","\n")
  print(effectiveSize(chain[,display]))
  
  # return results
  invisible(chain)
  
}

###############################################################
### load chains from Bayes Traits
test_1 = tracePlots(file=paste("./BTDat/BTJaws1stChainPPCA/BTselPPC.txt.VarRates.txt"),
                    burnin=4000)
test_2 = tracePlots(file=paste("./BTDat/BTJaws2ndChainPPCA/BTselPPC.txt.VarRates.txt"),
                    burnin=4000)

my_list_of_chains = mcmc.list(list(test_1[,c("Lh", "Lh...Prior", 
                                             "No.Pram", "Alpha", "Sigma.2")],
                                   test_2[,c("Lh", "Lh...Prior", 
                                             "No.Pram", "Alpha", "Sigma.2")]))
#Then simply use any functions from "coda" that works with objects of class "mcmc.list"

### supplementary figure S2
plot(my_list_of_chains)

#Gelman and Rubin's convergence diagnostic
gelman.diag(my_list_of_chains, confidence = 0.95, transform=FALSE, 
            autoburnin=TRUE,
            multivariate=TRUE)

effectiveSize(my_list_of_chains) 

#Plotting rate onto tree#
Phy <- Phy.sel

{
  ## load tree
  tree=Phy
  
  ## EVOLUTIONARY MODELS ##
  # BAYES TRAITS- JAW ##
  bt.tree <- tree
  bt.tree <- reorder(bt.tree, "cladewise")
  bt.tree.node.no <- makeNodeLabel(bt.tree)
  rate.heatmap <- colorRampPalette(c("#440154FF","#3B528BFF",
                                     "#21908CFF","#5DC863FF","#FDE725FF"),
                                   bias=1)
  plot(bt.tree, use.edge.length = F)
  nodelabels(text=bt.tree.node.no$node.label, frame="none", cex=0.5)
}

L_trad_tree_BTraits_tout<-BTRTools::rjpp(rjlog = paste("./BTDat/BTJaws2ndChainPPCA/BTselPPC.txt.VarRates.txt", sep=""),
                                         rjtrees = paste("./BTDat/BTJaws2ndChainPPCA/BTselPPC.txt.Output.trees", sep=""),
                                         tree = bt.tree, thinning = 1, burnin=4000)

cophylotrees <- cophylo(bt.tree,L_trad_tree_BTraits_tout$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_trad_tout.a <- cophylotrees$trees[[1]]

plotTree(tree_trad_tout.a, fsize = 0.5, ftype = "i")

#### Figure 5a
plotBranchbyTrait(tree_trad_tout.a, 
                  log(L_trad_tree_BTraits_tout$data$meanRate)[-1],
                  mode="edges", cex = 0.5,
                  palette=colorRampPalette(c("#440154FF","#3B528BFF",
                                             "#21908CFF","#5DC863FF",
                                             "#FDE725FF")))

##### plot shift at the node
new.bt.tree<-plotShifts(L_trad_tree_BTraits_tout, scalebar = TRUE,
                        scalar="node", tips=TRUE, cex=0.5,
                        shp = 24, colour = "grey")


#### Histogram rates

df <- (log(L_trad_tree_BTraits_tout$data$meanRate)[-1])

mdf <- matrix(df, nrow = 176, ncol = 1)
colnames(mdf) <- c("X")
write.csv(df,file = "ER_BR_1011_1.csv")

df2 <- read.table("ER_BR_1011_1.csv", header = TRUE, sep = ",")


yticks <- c(0,20,50,60,70,80)

h <- ggplot(data = df2, aes(x = x)) +
  geom_histogram(bins = 15, color = NA, fill = viridis(15))+
  xlab("log(meanRates)")+
  theme_pubr()

### Figure 5b
h

###############################################################
###############################################################
library(Rphylopars)
recon <- anc.recon(trait_data=Jaws.land,tree=Phy)

plot(Phy,no.margin = TRUE, cex=0.5)
nodelabels(text=1:Phy$Nnode,node=1:Phy$Nnode+Ntip(Phy))

### obtain the node for the reconstruction
Ela.anc<-arrayspecs(recon,100,3)

#### Supplementary Figure S11
plotRefToTarget(Ela.anc[,,2],
                Ela.anc[,,2], 
                label = F, mag = 1,
                method = "vector",
                gridPars = gridPar(pt.bg = "firebrick3"))

###############################################################
################ ratematrix ###################################
###############################################################

##### Select components to compare rates
PaPCS <- Sel.paca$x
PaPCS <- PaPCS[match(Phy.sel$tip.label, rownames(PaPCS)),]
RMPCS <- PaPCS[,1:2]

###############################################################
###### Example with Habitat as factor

map.sym.hab <- make.simmap(Sel.trees,Hab,model="SYM")

###############################################################
rmmcmc.hab3.1 <- ratematrixMCMC(data = RMPCS, 
                                phy  = map.sym.hab,
                                gen  = 10000000, dir=tempdir())

rmmcmc.hab3.2 <- ratematrixMCMC(data = RMPCS, 
                                phy  = map.sym.hab,
                                gen  = 10000000, dir=tempdir())

post.hab3.1 <- readMCMC(rmmcmc.hab3.1, burn = 0.1, thin = 1000)
post.hab3.2 <- readMCMC(rmmcmc.hab3.2, burn = 0.1, thin = 1000)

checkConvergence(post.hab3.1, post.hab3.2)

merg.post.hab3 <- mergePosterior(post.hab3.1, post.hab3.2)

#### Supplementary figure S10
plotRatematrix(merg.post.hab3, 
               colors = c("#000004FF","#781C6DFF","#ED6925FF","#FCFFA4FF"))


plotRootValue(merg.post.hab3)
testRatematrix(chain=merg.post.hab3, par="all")
testRatematrix(chain=merg.post.hab3, par="correlation")
testRatematrix(chain=merg.post.hab3, par="rates")

###############################################################
###############################################################
#### Reply results
###############################################################
###############################################################

### remove Oxynotus and Isistius to check differences due to outliers
Red.set <- Jaws3d.lm[,,c(-74,-98,-99,-100)]

Red.gpa<-gpagen(Red.set,curves = Jcurves,surfaces = JSurf,
                 ProcD = F)

Red.cla <- class.Raw[-c(74,98,99,100), ]
Red.cla <- droplevels(Red.cla)
gp2 <- as.factor(paste(Red.cla$Order))
col2.gp <- c("darkseagreen4","black","cyan3","darkorchid",
             "gold","darkorange","lightsteelblue",
             "deepskyblue4","firebrick3","navy")
names(col2.gp) <- levels(gp2)
col2.gp <- col2.gp[match(gp2,names(col2.gp))]

#### PCA
Red.pc.all <- gm.prcomp(Red.gpa$coords)
plot(Red.pc.all, axis1 = 1, axis2 = 2, pch=21,
     cex = 1.5, bg=col2.gp)

#### PaCA
Red.ave<-two.d.array(Red.gpa$coords)
means<-(aggregate(Red.ave~Red.cla$TipLab, FUN = mean))[,-1]
rownames(means)<-levels(Red.cla$TipLab)
Red.spe<-arrayspecs(means,100,3)

set.seed(183)
Phy<-read.tree("JTree.tre")

Red.c.sp <- class.Tre[-c(50,66), ]
Red.c.sp <- droplevels(Red.c.sp)

species <- Red.c.sp$TipLab
Phy.red <- drop.tip(Phy,Phy$tip.label
                   [-match(species, Phy$tip.label)])

gpRed <- as.factor(paste(Red.c.sp$Order))
col.m.gp <- c("darkseagreen4","black","cyan3","darkorchid",
              "gold","darkorange","lightsteelblue",
              "deepskyblue4","firebrick3","navy")
names(col.m.gp) <- levels(gpRed)
col.m.gp <- col.m.gp[match(gpRed,names(col.m.gp))]

################################

Red.av.pc <- gm.prcomp(Red.spe, phy = Phy.red,
                        align.to.phy = TRUE)

summary(Red.av.pc)

plot(Red.av.pc, axis1 = 1, axis2 = 2, pch=21, cex = 1.5,
     xlab = "PaCA 1: 67.91 %", ylab = "PaCA 26.51 %",
     bg=col.m.gp, phylo = TRUE, time.plot = FALSE,
     phylo.par = list(edge.color="grey50", anc.states =FALSE,
                      node.cex=0, node.labels = FALSE,
                      tip.labels = FALSE,
                      tip.txt.cex = 0.5))

### dispRity
Sub.sel <- subset(Red.c.sp, GP1 == "Sel")
Sub.sel <- droplevels(Sub.sel)

species <- Sub.sel$TipLab
Phy.r.s <- drop.tip(Phy,Phy$tip.label
                    [-match(species, Phy$tip.label)])

Sel.r.land <- Red.spe[,,-7]
Array2D <- two.d.array(Sel.r.land[,,c(1:87)])
Red.s.land <- Array2D[Phy.r.s$tip.label,]

average.sq <- function(X) mean(pairwise.dist(X)^2)

disp.nout <- dispRity.per.group(Red.s.land,
                              list(Hexa=c(85:87),
                                   Squa=c(68:79),
                                   Squt=c(81:83),
                                   Orec=c(52:65),
                                   Lamn=c(1:11),
                                   Carc=c(12:51)),
                              metric=c(average.sq))

summary(disp.nout)

test.dispRity(disp.nout, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(disp.nout, ylab = "Ave. Squared Distances")

#### change metric to c(sum, variances) to observe the differences between metrics

#### different habitat groups
#### fishbase classification with all the shark specimens
dispHab2 <- dispRity.per.group(as.matrix(Sel.land),
                               list(BatD=c(3,19:22,31:35,45,48,55,65,71),
                                    BatP=c(37,49),
                                    BenP=c(38:40,80,83),
                                    Deme=c(13,15,23,25,26,29,42,43,46,47,50,56:59,
                                           62,67,69,70,72,74,78,81,82,84:86,89),
                                    Pelag=c(1,2,18,24,51:54,66,68,75,76,79),
                                    Reef=c(4:12,14,16,17,27,28,30,36,41,44,60,61,
                                           63,64,73,77,87,88)),
                               metric=c(sum, variances))
summary(dispHab2)

#### fishbase classification without small group number
dispHab2 <- dispRity.per.group(as.matrix(Sel.land),
                               list(BatD=c(3,19:22,31:35,45,48,55,65,71),
                                    Deme=c(13,15,23,25,26,29,42,43,46,47,50,56:59,
                                           62,67,69,70,72,74,78,81,82,84:86,89),
                                    Pelag=c(1,2,18,24,51:54,66,68,75,76,79),
                                    Reef=c(4:12,14,16,17,27,28,30,36,41,44,60,61,
                                           63,64,73,77,87,88)),
                               metric=c(sum, variances))
summary(dispHab2)

#### Martinez et al (2021) classification
dispHab3 <- dispRity.per.group(as.matrix(Sel.land),
                               list(Shal=c(2,4:18,24:30,36,38:44,46,47,50,51,53,54,
                                           56:64,66,69,70,72:79,83,86:89),
                                    Twil=c(1,19,20,22,23,32,34,35,45,48,49,52,55,65,
                                           68,71,80:82,84,85),
                                    Deep=c(3,21,31,33,37,67)),
                               metric=c(sum, variances))
summary(dispHab3)

test.dispRity(dispHab3, test = wilcox.test, 
              comparison = "pairwise",
              correction = "bonferroni")
plot(dispHab3, ylab = "Procrustes Variance")

###########################################################################
#### Mechanical Advantage and Torsional Resistance in relation to
#### principal components Supplementaru figures S12 & S13
###########################################################################  
### torsinal resistance 
lmks <- matrix(c(5,6,1,5), ncol=2, byrow=TRUE,
               dimnames = list(c("Symphysis", "DGL"),c("start", "end")))
A <- jaws3d.gpa$coords
lineardists <- interlmkdist(A, lmks) 
TR <- data.frame(lineardists)
TR2<-TR/TR[,1]
### mechanical advantage
lmks2 <- matrix(c(1,6,1,2), ncol=2, byrow=TRUE,
               dimnames = list(c("ILD1", "ILD2"),c("start", "end")))
B <- jaws3d.gpa$coords
lineardists2 <- interlmkdist(B, lmks2) 
MA <- data.frame(lineardists2)
MA2<-MA/MA[,1]
MA2<-MA[,1]/MA


df <- data.frame(Ord = class.Raw$Order,
                 ID  = class.Raw$ID,
                 Sp2 = class.Raw$Species2,
                 PC1 = Jaws.pc.all$x[,1],
                 PC2 = Jaws.pc.all$x[,2],
                 PC3 = Jaws.pc.all$x[,3],
                 LTR = log(TR2[,2]),
                 TR = TR2[,2],
                 LMA = log(MA2[,2]),
                 MA = MA2[,2])

library(car)
library(ggplot2)
library(ggpubr)
library(ggrepel)

sp <- ggscatter(df, x = "MA", y = "PC1", color = "Ord",
                add = "reg.line",  
                add.params = list(color = "blue",
                                  fill = "lightgray"), 
                conf.int = TRUE, 
                xlab = "MA", ylab = "PC1")+
  scale_color_manual(values=c("darkseagreen4","black","cyan3","darkorchid",
                              "gold","darkorange","lightsteelblue",
                              "deepskyblue4","firebrick3","navy"))+
  geom_text_repel(mapping =  aes(label = Sp2, fontface = "italic"),
                  colour = "black", size = 2)
sp
sp + stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
              show.legend = FALSE, method = "pearson")


#### finally, report evolutionary rates as fold change
fold.change(mean(ER_all$Squt),mean(ER_all$Orec))
fold.change(mean(ER_all$Hexa),mean(ER_all$Orec))
fold.change(mean(ER_all$Lamn),mean(ER_all$Orec))
fold.change(mean(ER_all$Carc),mean(ER_all$Orec))
fold.change(mean(ER_all$Squa),mean(ER_all$Orec))

fold.change(mean(ER_all$BP),mean(ER_all$CEPH))
fold.change(mean(ER_all$CR),mean(ER_all$CEPH))
fold.change(mean(ER_all$FISH),mean(ER_all$CEPH))
fold.change(mean(ER_all$GEN),mean(ER_all$CEPH))

fold.change(mean(ER_all$TP),mean(ER_all$LP))
fold.change(mean(ER_all$MP),mean(ER_all$LP))

fold.change(mean(ER_all$reef),mean(ER_all$deep))
fold.change(mean(ER_all$shelf),mean(ER_all$deep))
fold.change(mean(ER_all$pelagic),mean(ER_all$deep))

