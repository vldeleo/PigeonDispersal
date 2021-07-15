#######################################
########### Bursera ###################
#######################################

# how to get statistics and figures presented in the manuscript
# before running these, make sure the genotype files are filtered, biallelic, and deduplicated
	# see other scripts for how I did that

# set up R
library("phyclust") # for making trees
library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")
library("vegan")
library("scales")
library("vcfR") #used for Gst/Hs stats
library("diveRsity") #used for Fis stats especially
#for fst map:
library(fields)
library(raster)
library(sp)
library(SNPRelate)
library(geosphere)
library(RColorBrewer)
library(prettymapr)
library(ggplot2)
#
# a note on my color palette:
# A was "royalblue2" (#436eee) now #002766
# B is "gray43" (#6E6E6E, light) now #FFDD33 now #FFE047
# C is "deepskyblue" now #007EF5
# D is "gray23" (#3B3B3B, dark) now #F5C000
# E is "slateblue3" (#6959cd) now #004CA3
# adults are filled/patterned, young are open shapes

	#I think I want to switch A and E, but also that seems like a lot of work now

#
ocean <- shapefile('/home/vqd5075/ne_10m_ocean.shp') # a more precise outline of coastlines
colonies <- shapefile('/home/vqd5075/PigeonBanks_Jaragua_2011-2016.shp') # mapped shapefiles of colonies
colonies$color <- c("orangered", "red2", "brown")[match(colonies$year, c('2015', '2011', '2016'))]
allforest <- raster('/home/vqd5075/treecover2010_20N_080W.tif') # % tree cover, to give context on forest fragmentation
#
locations <- read.csv("/media/elgon/victoria/NY2020/Bursera/renamedfiles3/SampleCoordinatesApproximate.csv") # coordinates for my samples
locations <- locations[which(locations$Species == "Bursera simaruba"),2:8]
locations$Location <- as.character(locations$Location)
	#I'm going to reset the names to A-E so that it's easy to match up with my genotype data
locations$Location[locations$Location == "WCP Bank, 2013, small"] <- "A"
locations$Location[locations$Location == "Parque Jaragua seed traps"] <- "B"
locations$Location[locations$Location == "Maipiero Stop 1"] <- "C"
locations$Location[locations$Location == "Maipiero Stop 2"] <- "D"
locations$Location[locations$Location == "forest"] <- "E"
locations <- locations[which(locations$Location %in% c("A", "B", "C", "D", "E")),]
locations$Indv <- tolower(sapply(strsplit(as.character(locations$Sample.number), " "), `[`, 2))
locations$matchingID <- paste(locations$Location, locations$Indv, sep = "_") 
locations$Age <- substring(locations$Indv, 1, 1)
#
pop_list <- read.table("/media/elgon/victoria/NY2020/Bursera/renamedfiles3/popmap")
names(pop_list) <- c("INDIVIDUALS", "POPULATION")
pop_code <- pop_list[c(1:9,11,13:26,28:38,41:45),] #leave out Met Browneii individuals By9 and Cm14 (in this case, leave out low quality individuals)
pop_code$AGE <- locations$Age[match(pop_code$INDIVIDUALS, locations$matchingID)]
pop_code$LAT <- locations$Lat[match(pop_code$INDIVIDUALS, locations$matchingID)]
pop_code$LONG <- locations$Lon[match(pop_code$INDIVIDUALS, locations$matchingID)]
#
# genotypes
Bur <- "/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/m33c90h01BiDD.vcf"
BurV <- read.vcf(Bur, which.loci = 1:1e7)
BurV2 <- read.vcfR(Bur)

#genind
BurG1 <- loci2genind(BurV)
strata(BurG1) <- as.data.frame(pop_code)
nameStrata(BurG1) <- ~Sample/Pop/Age/lat/long
setPop(BurG1) <- ~Pop 
BurG1@other$latlong <- pop_code[,4:5]
names(BurG1@other$latlong) <- c("lat", "long")
BurG1@other$data <- pop_code[,c(2,3)]

Bur_MAF <- minorAllele(BurG1)
summary(Bur_MAF)

#####################################
### Part 1: population statistics ###
#####################################
### Population Summary Statistics: Gst, Hs, pi ###

# Figure S???: Gst violin plots... use red/gold color palette to match Figure S1


   ### Geographic Distance ###
pop_code_dist <- pop_code[,c(1,4,5)] #must be formatted name, lat, lon
names(pop_code_dist) <- c("index", "lat", "lon")
	#the following function GeoDistanceinMetresMatrix comes from https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
BurGeoDist <- GeoDistanceInMetresMatrix(pop_code_dist)/1000 # am I sure I'm using great circles???
BurGeoDist2 <- BurGeoDist[-22, -22] #double check that y13 == individual 22


   ### Genetic Distance ###
# Euclidean distance between individuals
BurEuc <- dist(BurG1, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
BurEuc <- as.matrix(BurEuc)
write.csv(BurEuc, "Bur_EucDist_6.27.21.csv")

max(BurEuc, na.rm = T) #125.0343 max distance is less distant than in Metopium
min(BurEuc[BurEuc > 0], na.rm = T) #50.61731 

#without the parent offspring pair:
BurEuc2 <- dist(BurG1[-22], method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
BurEuc2 <- as.matrix(BurEuc2)
write.table(BurEuc2, "Bur_EucDist_noy13_6.27.10.csv")


#change dataformat
BurL <- genind2loci(BurG1) #technically, this seems like it should be the same as BurV

#genetic distance as number of loci that differ between individuals
BurDifLoci <- dist.gene(BurL, method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)
	#here, I drop loci that are missing in 1 individual
BurDifLoci <- as.matrix(BurDifLoci)
#write.csv(BurDifLoci, "Bur_LociDist_6.27.21.csv")
#without the parent/offspring pair?
BurDifLoci2 <- dist.gene(BurL[-22,], method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)
BurDifLoci2 <- as.matrix(BurDifLoci2)
#write.csv(BurDifLoci, "Bur_LociDist_noy13_6.27.21.csv")


#genetic distance as allelic difference
BurDifAllele <- diss.dist(BurG1, percent = FALSE, mat = FALSE)
BurDifAllele <- as.matrix(BurDifAllele)
#write.csv(BurDifAllele, "Bur_AlleleDist_4.16.20.csv")
BurDifAllele2 <- diss.dist(BurG1[-22], percent = FALSE, mat = FALSE)
BurDifAllele2 <- as.matrix(BurDifAllele2)
#write.csv(BurDifAllele2, "Bur_AlleleDist_noy13_4.16.20.csv")



   ### Fig1a: Isolation by Distance ###
#png(filename="Bur_GeoGenNucBiDD_5.27.21.png", res = 300, height = 1200, width = 1400)
pdf("Bur_GeoGenNucBiDD_6.27.21.pdf", height = 2.2, width = 2.5) #remembering that pdf default is in inches (this is more narrow than the 2.1x2.8 I had previously
par(mai=c(.5, .5, .25, .25), cex = .6, mgp=c(1.8,.8,0))
plot(as.matrix(BurGeoDist), as.matrix(BurEuc), xlab = "Geographic Distance (km)", ylab = "Genetic Distance (Euclidean)", col = alpha("black",.2), pch = 19) # if wanting to match Metopium, set ylim = c(40,160)
#remove self comparisons
#subset points that compare between Y and M
#index old
oldindex <- which(pop_code$AGE == "m")
youngindex <- which(pop_code$AGE == "y")
agegeo <- BurGeoDist
agegen <- BurEuc
agegeo[oldindex,oldindex] <- NA
agegeo[youngindex,youngindex] <- NA
agegen[oldindex,oldindex] <- NA
agegen[youngindex,youngindex] <- NA
#points(agegeo, agegen, col = alpha("red",.2), pch = 19)

allman <- mantel(BurGeoDist, BurEuc)$statistic #change numbers
#Mantel statistic r: 0.2113
#Significance: 0.002
ageman <- mantel(agegeo, agegen, na.rm = T)$statistic #I don't like doing this, because the help guide warns me about dropping NAs causing bias
#Mantel statistic r: 0.2398
#Significance: 0.002

points(agegeo[upper.tri(agegeo)], agegen[upper.tri(agegen)], col = alpha("red",.2), pch = 19)
points(BurGeoDist[BurEuc == 0], BurEuc[BurEuc == 0], col = "white", pch = 19)
legend("bottom", legend = c(paste("All Pairwise Comparisons: Mantel r = ", formatC(round(allman, digits = 3), format = "f", digits = 3), sep = ""), paste("Comparisons Across Age-Class: Mantel r = ", formatC(round(ageman, digits = 3), format="f", digits = 3), sep = "")), text.col = c(alpha("black", .5), alpha("red", .5)), cex = .8, bty = "n")
dev.off()

#
# now to make the plots without the parent/offspring
pdf("Bur_GeoGenNucBiDD_noy13_6.27.21.pdf", height = 2.2, width = 2.5) #remembering that pdf default is in inches... need to make more narrow
par(mai=c(.5, .5, .25, .25), cex = .6, mgp=c(1.8,.8,0))
plot(as.matrix(BurGeoDist2[BurGeoDist2 > 0]), as.matrix(BurEuc2[BurGeoDist2 > 0]), xlab = "Geographic Distance (km)", ylab = "Genetic Distance (Euclidean)", col = alpha("black",.2), pch = 19) #if I want to match Metopium, I can set "ylim = c(40,160)", but Jesse said not to
#remove self comparisons
#subset points that compare between Y and M
#index old
oldindex <- which(pop_code$AGE[-22] == "m")
youngindex <- which(pop_code$AGE[-22] == "y")
agegeo <- BurGeoDist[-22, -22]
agegen <- BurEuc2
agegeo[oldindex,oldindex] <- NA
agegeo[youngindex,youngindex] <- NA
agegen[oldindex,oldindex] <- NA
agegen[youngindex,youngindex] <- NA
#points(agegeo, agegen, col = alpha("red",.2), pch = 19)

allman <- mantel(BurGeoDist[-22,-22], BurEuc2)$statistic #change numbers
#Mantel statistic r: 0.1999
#Significance: 0.001
ageman <- mantel(agegeo, agegen, na.rm = T)$statistic #I don't like doing this, because the help guide warns me about dropping NAs causing bias
#Mantel statistic r: 0.2164
#Significance: 0.006

points(agegeo[upper.tri(agegeo)], agegen[upper.tri(agegen)], col = alpha("red",.2), pch = 19)
points(BurGeoDist2[BurEuc2 == 0], BurEuc2[BurEuc2 == 0], col = "white", pch = 19)
legend("bottom", legend = c(paste("All Pairwise Comparisons: Mantel r = ", formatC(round(allman, digits = 3), format = "f", digits = 3), sep = ""), paste("Comparisons Across Age-Class: Mantel r = ", formatC(round(ageman, digits = 3), format = "f", digits = 3), sep = "")), text.col = c(alpha("black", .5), alpha("red", .5)), cex = .8, bty = "n")
dev.off()



  ### test population differences in genetic distances ###
#how to test dispersion, without a t-test
	#this uses PERMDISP2 in vegan to analyze multivariate homogeneity of group variances
BurEuc2 <- dist(BurG1[-22], method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	#because I don't want a matrix, respecify BurEuc2
BurBeta <- betadisper(BurEuc2, pop_code$POPULATION[-22]) #can specify type (spatial median or group centroid), default is spatial median
	#error logical subscript too long (I think because I forgot to drop the parent/offspring pair)
TukeyHSD # creates set of confidence intervals between mean distance to centroid)

anova(BurBeta) #maybe?
	# F value = 1.0057 p = 0.4181
BurpermBeta <- permutest(BurBeta, pairwise = T, permutations = 10000)
	# so, disappointingly, no significant differences between groups
#	Overall (F)     A-B (t)     A-C (t)     A-D (t)     A-E (t)     B-C (t)
#  1.0057262   1.8615910   0.4758812   0.2836977  -0.2331287  -1.1850483
#    B-D (t)     B-E (t)     C-D (t)     C-E (t)     D-E (t)
# -2.1885916  -1.4181520  -0.3380398  -0.5376787  -0.4600580

BurTukBeta <- TukeyHSD(BurBeta, pairwise = T)
	# also no differences here; lower and upper confidence intervals span 0 and p is high
#          diff        lwr       upr     p adj
#B-A -4.1098610 -11.140092  2.920370 0.4572506
#C-A -1.0782011  -7.559754  5.403352 0.9887776
#D-A -0.4510781  -6.579892  5.677735 0.9995272
#E-A  0.6808864  -6.349345  7.711118 0.9986062
#C-B  3.0316599  -4.172183 10.235503 0.7446204
#D-B  3.6587829  -3.229409 10.546975 0.5511428
#E-B  4.7907474  -2.910485 12.491980 0.3950115
#D-C  0.6271230  -5.700085  6.954331 0.9984742
#E-C  1.7590875  -5.444756  8.962931 0.9543663
#E-D  1.1319645  -5.756227  8.020156 0.9892815
	#there isn't a significant difference in genetic diversity (rather, genetic distances) among populations


   ### Fig1b: Neighbor-Joining tree of Genetic Distance ###
#if I want differently shaped tips for pop and fill for age, the most straightforward way might be to make a separate vector to specify that:
age <- c(rep(c("m", "y"),5))
pop <- c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E")
shape <- c(21, 1, 23, 5, 24, 2, 25, 6, 22, 0) #there's no basic filled down triangle??
shape2 <- c(16, 1, 18, 5, 17, 2, 6, 6, 15, 0)
plotshape <- cbind(age, pop, shape, shape2)
colnames(plotshape) <- c("AGE", "POPULATION", "SHAPE", "SHAPE2")
popshape <- merge(pop_code, plotshape)
popshape$SHAPE <- as.numeric(as.character(popshape$SHAPE))
popshape$SHAPE2 <- as.numeric(as.character(popshape$SHAPE2))
popshape$COLOR <- c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(popshape$POPULATION, c('A', 'B', 'C', 'D', 'E'))]


#png(filename="Bur_EucNJUnrootedBiDD_6.25.21.png", res = 300, height = 1500, width = 1500)
par(mar=c(.6, .6, .6, .6))
plotnj(nj(as.matrix(BurEuc)), show.tip.label =F) #a distance matrix
tiplabels(pch = popshape$SHAPE[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], col =  popshape$COLOR[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], bg = popshape$COLOR[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], lwd = 5, cex = 5)
#legend("bottomleft", legend = c("Juvenile", "Adult"), pch = c(1, 16), col = "#F5C000", bty = "n", cex = 1.2) # I've left off the legend here because I have a standardized legend I use for all the plots
#dev.off()


#I'd like to make friendly for b/w printing, but it is unreasonably difficult to plot patterned fill in R
	#instead, I'm going to switch to a blue/gold color palette (bc blue/black wasn't obvious enough distinction for colony/not)
	#and I've added shapes to differentiate a bit
png(filename="Bur_EucNJUnrootedBiDD_6.27.21.png", res = 300, height = 1500, width = 1500)
par(mar=c(.6, .6, .6, .6))
plotnj(nj(as.matrix(BurEuc)), show.tip.label =F) 
tiplabels(pch = popshape$SHAPE[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], col =  popshape$COLOR[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], bg = popshape$COLOR[match(as.character(rownames(BurEuc)), popshape$INDIVIDUALS)], lwd = 5, cex = 5)
#tiplabels(text = pop_code$INDIVIDUALS, frame = 'none', bg = NA, cex = 0.3, col = gray(0.2))
#legend("bottomleft", legend = c("Juvenile", "Adult"), pch = c(1, 16), col = "#F5C000", bty = "n", cex = 1.2) # I've left off the legend here because I have a standardized legend I use for all the plots
dev.off()
	

# what about the likely parent/offspring pair C_m13 and C_y13?
	#I can drop C_y13 to see how the tree changes
#min(BurEuc2[BurEuc2 > 0], na.rm = T) #82.22552 #so, more distant, but not radically so
#png(filename="Bur_EucNJUnrootedBiDD_noy13_5.27.21.png", res = 300, height = 1500, width = 1500)
pdf("Bur_EucNJUnrootedBiDD_noy13_5.27.21.pdf")
par(mar=c(.6, .6, .6, .6))
plotnj(nj(as.matrix(BurEuc2)), show.tip.label =F) 
tiplabels(pch = popshape$SHAPE[match(as.character(rownames(BurEuc2)), popshape$INDIVIDUALS)], col =  popshape$COLOR[match(as.character(rownames(BurEuc2)), popshape$INDIVIDUALS)], bg = popshape$COLOR[match(as.character(rownames(BurEuc2)), popshape$INDIVIDUALS)], lwd = 5, cex = 5)
dev.off()


#nj trees of the other measures of genetic distances (not included in manuscript, but here for comparison)
#
plotnj(nj(BurDifAllele), show.tip.label =F) 
tiplabels(pch = c(1, 16)[match(pop_code$AGE, c('y', 'm'))], col =  c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(pop_code$POPULATION, c('A', 'B', 'C', 'D', 'E'))], bg = 'white', lwd = 2, cex =2)
	#
plotnj(nj(BurDifAllele2), show.tip.label =F) 
tiplabels(pch = c(1, 16)[match(pop_code$AGE[-22], c('y', 'm'))], col =  c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(pop_code$POPULATION[-22], c('A', 'B', 'C', 'D', 'E'))], bg = 'white', lwd = 2, cex =2)
# loci dif
plotnj(nj(BurDifLoci), show.tip.label =F) # a distance matrix
tiplabels(pch = c(1, 16)[match(pop_code$AGE, c('y', 'm'))], col =  c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(pop_code$POPULATION, c('A', 'B', 'C', 'D', 'E'))], bg = 'white', lwd = 2, cex =2)
	# loci dif, without offspring
plotnj(nj(BurDifLoci2), show.tip.label =F) # a distance matrix
tiplabels(pch = c(1, 16)[match(pop_code$AGE[-22], c('y', 'm'))], col =  c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(pop_code$POPULATION[-22], c('A', 'B', 'C', 'D', 'E'))], bg = 'white', lwd = 2, cex =2)

# making a legend
pdf("BlueGoldLegend.pdf", height = 3, width = 3)
par(mar=c(.6, .6, .6, .6))
#plotnj(nj(BurEuc), show.tip.label =F) 
plot.new()
#legend("bottomleft", legend = c("Juvenile", "Adult"), pch = c(1, 16, 21, 22, 24, 23, 25), col = "gray23", bty = "n", cex = 1.2)
legend("bottomright", legend = c("Juvenile", "Adult", "2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony"), pch = c(1, 16, 21, 22, 24, 23, 25), col =  c("black", "black", "#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000"), bty = "n", cex = 1.2)
dev.off()
# this legend takes some fixing in Inkscape to make the shapes filled and to reveal the cutoff sections



   ### Figure 1c: PCoA of Genetic Distances ###
#pdf("Bur_PCoABiDD_6.27.21_2.pdf", height = 1.5, width = 1.5)
par(cex = c(.7))
BurMDS <- scaleGen(BurG1, scale=FALSE, NA.method="mean") 
pco.bur <- dudi.pco(dist(BurMDS), scannf=FALSE, nf=3)
#s.class(pco.bur$li, fac=pop(BurG1), xax=1, yax=2, col=transp(c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3"), .75), axesel=FALSE, cstar=0, cpoint=3, clabel = "", addaxes = TRUE, pch = c(1, 16)[match(pop_code$AGE, c('y', 'm'))], xlim = c(min(pco.bur$li$A1) - 8, max(pco.bur$li$A1) + 5), ylim = c(min(pco.bur$li$A2) - 5, max(pco.bur$li$A2) + 5)) 
	#okay, what I've been doing wrong is that I specified the fillable shape instead of the pre-filled shape
	#I did this because there is no basic upside down filled triangle. however, this means that the bg command doesn't work.
	#in the end, I made a new vector using basic empty and filled shapes, but had to fill the upside down triangles manually based on where the filled circles from ^ were on my plot
s.class(pco.bur$li, fac=pop(BurG1), xax=1, yax=2, col=transp(c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3"), .75), axesel=FALSE, cstar=0, cpoint=3, clabel = "", addaxes = TRUE, pch = popshape$SHAPE2[match(rownames(pco.bur$tab), popshape$INDIVIDUALS)], xlim = c(min(pco.bur$li$A1) - 8, max(pco.bur$li$A1) + 5), ylim = c(min(pco.bur$li$A2) - 5, max(pco.bur$li$A2) + 5))  # bg = popshape$COLOR[match(rownames(pco.bur$tab), popshape$INDIVIDUALS)], don't know if I need this
#s.class(pco.bur$li, fac=paste(popshape$POPULATION, popshape$AGE), xax=1, yax=2, col=transp(c("#002766","#002766", "#FFE047","#FFE047", "#007EF5","#007EF5", "#F5C000","#F5C000", "#004CA3", "#004CA3"), .75), axesel=FALSE, cstar=0, cpoint=3, clabel = "", addaxes = TRUE, pch = popshape$SHAPE2[match(rownames(pco.bur$tab), popshape$INDIVIDUALS)], bg = transp(c("#002766",NA,"#FFE047",NA,"#007EF5",NA,"#F5C000",NA,"#004CA3",NA), .75), xlim = c(min(pco.bur$li$A1) - 8, max(pco.bur$li$A1) + 5), ylim = c(min(pco.bur$li$A2) - 5, max(pco.bur$li$A2) + 5))  # bg = popshape$COLOR[match(rownames(pco.bur$tab), popshape$INDIVIDUALS)], don't know if I need this
#dev.off()
	#there's no obvious way to fix this.

#making shapes smaller so easier to see?
pdf("Bur_PCoABiDD_6.27.21_3.pdf", height = 1.5, width = 1.5)
par(cex = c(.6))
BurMDS <- scaleGen(BurG1, scale=FALSE, NA.method="mean") 
pco.bur <- dudi.pco(dist(BurMDS), scannf=FALSE, nf=3)
s.class(pco.bur$li, fac=pop(BurG1), xax=1, yax=2, col=transp(c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3"), .75), axesel=FALSE, cstar=0, cpoint=3, clabel = "", addaxes = TRUE, pch = popshape$SHAPE2[match(rownames(pco.bur$tab), popshape$INDIVIDUALS)], xlim = c(min(pco.bur$li$A1) - 8, max(pco.bur$li$A1) + 5), ylim = c(min(pco.bur$li$A2) - 5, max(pco.bur$li$A2) + 5))  
dev.off()


# could I do this with plain ole Euclidean distances? (no, but not sure why. I think I need to change the function I'm using for that, and ade4 is complicated.)
pdf("Bur_PCoAEucBiDD_6.4.21.pdf", height = 1.5, width = 1.5)
par(cex = c(.7))
pco.bur <- dudi.pco(BurEuc, scannf=FALSE, nf=3)
s.class(pco.bur$li, fac=pop(BurG1), xax=1, yax=2, col=transp(c("002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3"), .75), axesel=FALSE, cstar=0, cpoint=3, clabel = "", addaxes = TRUE, pch = c(1, 16)[match(pop_code$AGE, c('y', 'm'))], xlim = c(min(pco.bur$li$A1) - 8, max(pco.bur$li$A1) + 5), ylim = c(min(pco.bur$li$A2) - 5, max(pco.bur$li$A2) + 5)) 
#dev.off()



   ### Figure 1d: Loci Fis ###
#see code for 83% CI bootstrap
	#to get 95% confidence intervals
CIstats <- divBasic("/media/bale/victoria/NY2020/Bursera/filtertests/nucreads/missing33/01_radiator_genomic_converter_20200212@1331/m33c90h01BiDD_genepop.gen", bootstraps = 1000, fis_ci = TRUE)
GPstats <- divBasic("/media/bale/victoria/NY2020/Bursera/filtertests/nucreads/missing33/01_radiator_genomic_converter_20200212@1331/m33c90h01BiDD_genepop.gen", bootstraps = 1000, fis_ci = TRUE)
write.csv(GPstats$fis, "Bur_FisBiDD_5.4.20.csv")
ciTable <- lapply(GPstats$fis, function(x){
  return(c(x["overall", "fis"], x["overall", "BC_lower_CI"],
           x["overall", "BC_upper_CI"]))
})

	#to get 83% confidence intervals
	# first, grab divBasic83 function from divBasic83 script
GPstats83 <- divBasic83("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/01_radiator_genomic_converter_20200212@1331/m33c90h01BiDD_genepop.gen", bootstraps = 1000)
write.csv(GPstats83$fis, "Bur_FisBiDD83_6.28.21.csv")
	#followed same code as for 95% to make a whiskers plot with 83% confidence intervals
# convert this into a dataframe
ciTable <- lapply(GPstats83$fis, function(x){
  return(c(x["overall", "fis"], x["overall", "BC_lower_CI"],
           x["overall", "BC_upper_CI"]))
})
ciTable <- as.data.frame(do.call("rbind", ciTable))
dimnames(ciTable) <- list(paste("pop_", 1:nrow(ciTable), sep = ""),
                          c("Fis", "lower", "upper"))
ciTable

#if have done divBasic already (it takes forever), can load as
#GPstats83 <- read.csv("Bur_FisBiDD83_6.28.21.csv")

write.csv(ciTable, "Bur_FisCIBiDD83_6.28.21.csv")
#ciTable <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/RAnalysis/Bur_FisCIBiDD83_6.28.21.csv")
library(data.table)
#change population names?
ciTable$pop <- c("2011 Colony", "Parque Jaragua NC", "Mapioró C", "Mapioró NC", "Parque Jaragua C")
setDT(ciTable, keep.rownames = TRUE)[]

ggplot(ciTable, aes(x = pop, y = Fis)) +
  coord_cartesian(ylim = c(-.1, .3)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90)) + xlab("")

#rawpoints <- head(GPstats83$fis, -1)
rawpoints <- GPstats83$fis #this works if you're using GPstats directly from the divBasic function
		#if you've read GPstats$fis back in after saving as csv, code will be a little different
		#might want to skip down to the rawlong file if you've already saved it
raw1 <- as.data.frame(GPstats83$fis[1]) #or, raw1 <- as.data.frame(GPstats83[,2])
raw2 <- as.data.frame(GPstats83$fis[2]) #raw2 <- as.data.frame(GPstats83[,7])
raw3 <- as.data.frame(GPstats83$fis[3]) #raw3 <- as.data.frame(GPstats83[,12])
raw4 <- as.data.frame(GPstats83$fis[4]) #raw4 <- as.data.frame(GPstats83[,17])
raw5 <- as.data.frame(GPstats83$fis[5]) #raw5 <- as.data.frame(GPstats83[,22])
rawpoints <- as.data.frame(cbind(raw1[,1], raw2[,1], raw3[,1], raw4[,1], raw5[,1]))
names(rawpoints)[1:5] <- c("2011 Colony", "Parque Jaragua NC", "Mapioró C", "Mapioró NC", "Parque Jaragua C")
#perfect, now to make it long format
library(tidyr)
rawlong <- gather(rawpoints, pop, Fis, "2011 Colony":"Parque Jaragua C", factor_key=TRUE)
rawlong$color <- c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3")[match(rawlong$pop, c("2011 Colony", "Parque Jaragua NC", "Mapioró C", "Mapioró NC", "Parque Jaragua C"))]
write.csv(rawlong, "Bur_FisLociLong83_6.30.21.csv")
	#fwiw, I tried to do this within the ggplot2 call, but something was not working

group.colors <- c("2011 Colony" = "#002766", "Parque Jaragua NC" = "#FFE047", "Mapioró C" ="#007EF5", "Mapioró NC" = "#F5C000", "Parque Jaragua C" = "#004CA3")
	#make sure the populations show up on the plot in the order I want:
#rawlong$popfake <- c("A", "B", "C", "D", "E")[match(rawlong$pop, c("2011 Colony", "Parque Jaragua C", "Mapioró C", "Parque Jaragua NC", "Mapioró NC"))]
#ciTable$popfake <- c("A", "B", "C", "D", "E")[match(ciTable$pop, c("2011 Colony", "Parque Jaragua C", "Mapioró C", "Parque Jaragua NC", "Mapioró NC"))]
rawlong$pop <- factor(rawlong$pop, levels = c("2011 Colony", "Parque Jaragua C", "Mapioró C", "Parque Jaragua NC", "Mapioró NC"))
ciTable$pop <- as.factor(ciTable$pop)
ciTable$pop <- factor(ciTable$pop, levels = c("2011 Colony", "Parque Jaragua C", "Mapioró C", "Parque Jaragua NC", "Mapioró NC"))

#need to set y axis to (-.1 to .3) for both in order to compare across species
	#might not be useful

#make pdf
#rawlong <- read.csv("Bur_FisLociLong83_6.30.21.csv")
#ciTable <- read.csv("Bur_FisCIBiDD_6.28.21.csv")[,c(2:5)]
library(data.table)
#fix formatting, if reading in from saved file 
rownames(ciTable) <- ciTable[,1]
ciTable$pop <- c("2011 Colony", "Parque Jaragua NC", "Mapioro C", "Mapioro NC", "Parque Jaragua C")
setDT(ciTable, keep.rownames = TRUE)


#since 95% confidence intervals are more conservative than 5% significance,
	#I had to calculate 83% confidence intervals from my bootstrap
pdf("Bur_popFis_BiDD83_6.27.21.pdf", height = 2, width = 3, pointsize = .5)
par(cex = .6)
ggplot(ciTable, aes(x = pop, y = Fis, group = pop)) +
  geom_point() + 
  theme_bw() + xlab("") + ylab("Loci Fis") + coord_cartesian(ylim = c(-.1, .3)) +  
  #geom_jitter(data = rawlong, mapping = aes(x = pop, y = Fis, color= pop)) +
  geom_violin(data = rawlong, aes(x = pop, y = Fis, color = pop))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  scale_color_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000")) +
#theme(axis.text.x=element_text(angle = 90, vjust = -.0001), legend.position="none") 
theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="none")
dev.off()



   ### Figure S8b: Hs ###
BurV2 <- read.vcfR(Bur)
BurD <- genetic_diff(BurV2, pops = as.factor(pop_code$POPULATION), method = "jost")
 knitr::kable(head(BurD[1:13]))
knitr::kable(round(colMeans(BurD[,c(3:13)], na.rm = TRUE), digits = 3))
#|Hs_A      | 0.220|	
#|Hs_B      | 0.213|	
#|Hs_C      | 0.215|	
#|Hs_D      | 0.217|	
#|Hs_E      | 0.208|	
#|a         | 3.830|	
#|b         | 3.850|	
#|Dest_Chao | 0.005|	
#|Da        | 1.327|	
#|Dg        | 1.363|	
#|Db        | 1.024|	
#|Ht       | 0.230|
#|Gst      | 0.064|
#|Gprimest | 0.088|

BurN <- genetic_diff(BurV2, pops = as.factor(pop_code$POPULATION), method = "nei")
 knitr::kable(head(BurN[1:13]))
knitr::kable(round(colMeans(BurN[,c(3:17)], na.rm = TRUE), digits = 3))


#this does successfully reorder the populations, but I wanted a violin...
library(tidyr)
BurDlong <- gather(BurD[,c(1:7)], population, H, Hs_A, Hs_B, Hs_C, Hs_D, Hs_E)
BurDlong$population <- factor(BurDlong$population, levels = c("Hs_A", "Hs_E", "Hs_B", "Hs_C", "Hs_D"))
BurDlong$population <- c("2011 Colony", "Parque Jaragua NC", "Mapioro C", "Mapioro NC", "Parque Jaragua C")
pdf("Bur_poH_pointsBiDD_6.28.21.pdf", height = 2, width = 3, pointsize = .5)
par(cex = .6)
ggplot(BurDlong[,3:4], aes(x=population, y=H, fill=population)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#002766", "#007EF5", "#F5C000", "#004CA3", "#FFE047"))
dev.off()
	
BurNlong <- gather(BurN[,c(1:8,17)], population, H, Hs_A, Hs_B, Hs_C, Hs_D, Hs_E, Ht, Gprimest)
BurNlong$population <- factor(BurNlong$population, levels = c("Hs_A", "Hs_E", "Hs_C", "Hs_B", "Hs_D", "Ht", "Gprimest"))
pdf("Bur_LociHsBiDD_6.28.21.pdf", height = 2, width = 3, pointsize = .5)
par(cex = .6)
ggplot(BurNlong[,3:4], aes(x=population, y=H, fill=population)) + 
	geom_violin() +
	scale_fill_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000", "darkgray", "darkgray")) +
	theme_bw() + theme(legend.position="none") +
	xlab("") + ylab("Loci Hs")
  	#geom_violin(data = BootL, aes(x = pop, y = Hs, color = pop))
dev.off()
	#I ended up fixing the labels in inkscape


## confidence intervals ##
#another way to calculate Hs (with significance test) is using adegenet
Hs(BurG1, pop = pop_code$POPULATION) #or pop_code?
#       A         B         C         D         E
#0.2203450 0.2132709 0.2153276 0.2165826 0.2083183
#gives the same result as the other function, which is comforting

Hs.test(BurG1[pop="A"], BurG1[pop="B"], n.sim = 1000, alter = "two-sided")
#hypothesis: pigeon sites have more heterozygosity, so could test "greater"
	#could also test "two-sided"
#these tests take a long time, so set them up and do something else
	#Monte-Carlo test Observation: 0.0070741
	#p-value: 0.6133866
	# Std.Obs.y  Expectation     Variance
	#7.059980e-01 2.337825e-03 4.500562e-05

Hs.test(BurG1[pop="A"], BurG1[pop="C"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.005017343, p-value = 0.2367632
Hs.test(BurG1[pop="A"], BurG1[pop="D"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.003762423, p-value = 0.3986014
Hs.test(BurG1[pop="A"], BurG1[pop="E"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.0120267, p-value = 0.2127872
Hs.test(BurG1[pop="B"], BurG1[pop="C"], n.sim = 1000, alter = "two-sided")
	#test Observation: -.0.002056757, p-value = 0.5774226
Hs.test(BurG1[pop="B"], BurG1[pop="D"], n.sim = 1000, alter = "two-sided")
	#test Observation: -0.003311677, p-value = 0.5124875
Hs.test(BurG1[pop="B"], BurG1[pop="E"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.004952601, p-value = 0.4325674
Hs.test(BurG1[pop="C"], BurG1[pop="D"], n.sim = 1000, alter = "two-sided")
	#test Observation: -0.00125492, p-value = 0.959041
Hs.test(BurG1[pop="C"], BurG1[pop="E"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.007009358, p-value = 0.4275724
Hs.test(BurG1[pop="D"], BurG1[pop="E"], n.sim = 1000, alter = "two-sided")
	#test Observation: 0.008264277, p-value = 0.6213786
	
###### ************correcting for multiple testing *****################





###### Figure 2e: Hs ########
# I used Hs and Hs.test in adegenet originally to test if differences were significant.
# Hs.test uses a permutation test of individuals in populations to see if Hs is different between populations
# when I tested in 2020, I used a two-sided test (fairly conservative, since I had a hypothesis that nesting sites 
	# would be more diverse)
# H0: Hs will be the same irrespective of nesting colonies
	Hs(popx) - Hs(popy) = 0
# Ha: Hs will be greater in nesting colonies
	Hs(popx) - Hs(popy) > 0

# I want to sample individuals (within populations, with replacement) and calculate Hs
	#but this is a bootstrap, not a permutation test, so maybe a little dif than the results from Hs.test

BootHs <- vector()
#BootHs <- read.table(BootHs_2.26.21, sep = ",", col.names = T)
#BootHs <- BootHs[2:499, 2:6]
#colnames(BootHs) <- c("A", "B", "C", "D", "E")
for (i in 1:1000){ # yeah, yeah, I know for loops are bad
mySamp <- lapply(seppop(BurG1_s), function(x) x[sample(1:nrow(x$tab), nrow(x$tab), replace = T)])
	row.names(mySamp$A@tab) <- paste("A", 1:nrow(mySamp$A@tab), sep = "_") #can't repool with duplicated sample rownames
	row.names(mySamp$B@tab) <- paste("B", 1:nrow(mySamp$B@tab), sep = "_")
	row.names(mySamp$C@tab) <- paste("C", 1:nrow(mySamp$C@tab), sep = "_")
	row.names(mySamp$D@tab) <- paste("D", 1:nrow(mySamp$D@tab), sep = "_")
	row.names(mySamp$E@tab) <- paste("E", 1:nrow(mySamp$E@tab), sep = "_")
SampRep <- repool(mySamp)
BootHs <- rbind(BootHs, Hs(SampRep))
write.csv(BootHs, "BootHs_2.26.21.csv") #because writing this as a for loop takes a stupidly long time, need to save for when I lose connection to server
print(i)
}
#cp BootHs_2.26.21 BootHs_2.26.21.csv #BootHs_2.26.21 is good up until iteration 499; not sure what happened after that

BootDF <- as.data.frame(BootHs)
for (i in 1:ncol(BootDF)){
BootDF[,i] <- as.numeric(as.character(BootDF[,i]))
}
HsCI83 <- lapply(BootDF, quantile, probs = c(0.083, 0.50,
            0.917), na.rm = TRUE) # how to break this down by population?

write.csv(HsCI83, "HsBootQ_2.26.21.csv")

# visualization 
	#make violin in the style of my Fis plot
	#make long for geom_violin
library(tidyr)
BootL <- gather(BootDF, pop, Hs, A:E, factor_key = TRUE)
CIdf <- as.data.frame(HsCI83) %>% t() %>% as.data.frame() #getting fancy with dplyr?
#CIdf$pop <- rownames(CIdf)
CIdf$pop <- c("2011 Colony", "Parque Jaragua Non-Colony", "Mapioró Colony", "Mapioró Non-Colony", "Parque Jaragua Colony")
colnames(CIdf) <- c("lower", "Hs", "upper", "pop")
BootL$pop <- as.character(BootL$pop)
BootL$pop[which(BootL$pop == "A")] <- "2011 Colony" 
BootL$pop[which(BootL$pop == "B")] <- "Parque Jaragua Non-Colony" 
BootL$pop[which(BootL$pop == "C")] <- "Mapioró Colony"
BootL$pop[which(BootL$pop == "D")] <- "Mapioró Non-Colony"
BootL$pop[which(BootL$pop == "E")] <- "Parque Jaragua Colony" 
BootL$pop <- factor(BootL$pop, levels = c("2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony")) #this FINALLY allows me to match the order on the Fis plot

rownames(CIdf) <- NULL
CIdf$pop <- as.factor(CIdf$pop)
CIdf$pop <- factor(CIdf$pop, levels = c("2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony")) #this FINALLY allows me to match the order on the Fis plot
#CIdf <- CIdf[c(1,3,4,5,2),]
#levels(CIdf$pop) <- c("2011 Colony", "Parque Jaragua Non-Colony", "Mapioro Colony", "Mapioro Non-Colony", "Parque Jaragua Colony")


png(filename = "Bur_popHs_violinBiDD3.2.21.png", res = 300, height = 1000, width = 1000)
ggplot(CIdf, aes(x = pop, y = Hs, group = pop)) + #I tried leaving out ", group = pop", but made no difference
  geom_point(size = .3) +  #this plots the median
  theme_bw() + xlab("") + ylab("Hs") + #coord_cartesian(ylim = c(-.1, .3)) +
  geom_violin(data = BootL, aes(x = pop, y = Hs, color = pop))+
  #scale_fill_manual(c("royalblue2", "gray43", "deepskyblue", "gray23", "slateblue3")) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + #this gives whiskers from 8.3% to 91.7% quantile
  #geom_pointrange(aes(ymin = lower, ymax = upper), data = CIdf)
  #scale_fill_manual(values=c(), drop = FALSE,   #this was SUPPOSED to give me color coding to match my other plots but I didn't figure it out
  #          name  ="Population",breaks=c(1,2,3,4,5),
  #          labels=c("2011 Colony", "Parque Jaragua NC", "Mapioro C", "Mapioro NC", "Parque Jaragua C")) +
  theme(axis.text.x=element_text(angle = 90), legend.position="none"))
dev.off()


pdf("Bur_popHs_violinBiDD6.29.21.pdf", height = 2, width = 3, pointsize = .5)
par(cex = .6)
ggplot(CIdf, aes(x = pop, y = Hs, group = pop)) +
  geom_point() +  #should plot the median
  theme_bw() + xlab("") + ylab("Hs") + #coord_cartesian(ylim = c(,)) +  
  geom_violin(data = BootL, aes(x = pop, y = Hs, color = pop))+
  scale_fill_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000")) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + #this gives whiskers from 8.3% to 91.7% quantile
  geom_pointrange(aes(ymin = lower, ymax = upper)) + #an attempt to get around my missing center point
  scale_color_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000")) +
#theme(axis.text.x=element_text(angle = 90, vjust = -.0001), legend.position="none") 
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), # not ready for this until I get the order correct
        axis.ticks.x=element_blank(), legend.position="none")
dev.off()



### alternatively ###
### compare pigeon vs no-pigeon and young vs old ###
# to minimize the problem of having small sample sizes, I can lump together all my colony and non-colony pops and my young and old (for colony/non-colony)












   ### Figure 1f: pi ###
# I can use nuc.div to calculate "variances", but I run into the same multiple testing problem

BurBin <- vcfR2DNAbin(BurV2, unphased_as_NA = FALSE)
pegas::nuc.div(BurBin, pairwise.deletion = T) #0.2286449

# nuc.div() only works with BurBin, but I'm not sure that I can actually convert to a DNAbin object
	# DNAbin is supposed to be good because it takes into account invariant locations, but if this is coming from a vcfR, there's no point
#pairwise deletion drops missing or ambiguous nucleotides
#nuc.div has a built in fuction to calculate variance, but only works if 
	# sequences are the same length

#if I don't use pairwise deletion (can't use if calculating variance), estimated pi is higher
pegas::nuc.div(BurBin[1:18,], variance = T) #0.22961600 0.01321669
pegas::nuc.div(BurBin[19:30,], variance = T) #0.22809328 0.01393375
pegas::nuc.div(BurBin[31:48,], variance = T) #0.22335726 0.01250619
pegas::nuc.div(BurBin[49:68,], variance = T) #0.22439401 0.01246252
pegas::nuc.div(BurBin[69:80,], variance = T) #0.2233198 0.0133568
	#not sig: A/B, A/C, A/D, A/E, B/C, B/D, B/E, all of them


# first, I resampled from each population a bunch of times to get an approximation of confidence intervals
# the code for this is in BootPiCodeCommandLine.sh
	# and the output is in a handful of files called eg PopA_BootstrapPi.txt
# then I can plot: (note, I originally typed up this code in PlottingTable1_sigtests.txt

#from inside the /SingleVCFs folder
PopA <- read.table("PopA_BootstrapPi.txt", header = F)
PopB <- read.table("PopB_BootstrapPi.txt", header = F)
PopC <- read.table("PopC_BootstrapPi.txt", header = F)
PopD <- read.table("PopD_BootstrapPi.txt", header = F)
PopE <- read.table("PopE_BootstrapPi.txt", header = F)
AllPi <- cbind (PopA, PopB, PopC, PopD, PopE)
names(AllPi) <- c("A", "B", "C", "D", "E")

PiCI83 <- lapply(AllPi, quantile, probs = c(0.083, 0.50,
            0.917), na.rm = TRUE) # how to break this down by population?
AllPiL <- gather(AllPi, pop, Pi, A:E, factor_key = TRUE)
CIPi <- as.data.frame(PiCI83) %>% t() %>% as.data.frame() #getting fancy with dplyr?
CIPi$pop <- rownames(CIPi)
#CIPi$pop <- factor(CIPi$pop, levels = c("A", "C", "D", "E", "B")) 
CIPi$pop <- c("2011 Colony", "Parque Jaragua Non-Colony", "Mapioró Colony", "Mapioró Non-Colony", "Parque Jaragua Colony")
CIPi$pop <- as.factor(CIPi$pop)
CIPi$pop <- factor(CIPi$pop, levels = c("2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony")) 
colnames(CIPi) <- c("lower", "Pi", "upper", "pop") #because it's tricky to use special characters like % in a column name

#get data in form that ggplot approves of
AllPiL$pop <- as.character(AllPiL$pop)
AllPiL$pop[which(AllPiL$pop == "A")] <- "2011 Colony" 
AllPiL$pop[which(AllPiL$pop == "B")] <- "Parque Jaragua Non-Colony" 
AllPiL$pop[which(AllPiL$pop == "C")] <- "Mapioró Colony"
AllPiL$pop[which(AllPiL$pop == "D")] <- "Mapioró Non-Colony"
AllPiL$pop[which(AllPiL$pop == "E")] <- "Parque Jaragua Colony" 
AllPiL$pop <- factor(AllPiL$pop, levels = c("2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony")) 


pdf("Bur_popPi_violinBiDD6.29.21.pdf", height = 2, width = 3, pointsize = .5)
par(cex = .6)
ggplot(CIPi, aes(x = pop, y = Pi, group = pop)) +
  geom_point() +  #should plot the median
  theme_bw() + xlab("") + ylab("Pi") + scale_y_continuous(breaks = scales::pretty_breaks(n=5)) + #coord_cartesian(ylim = c(,)) +  
				#I want the same number of tick marks as I get in other plots
  geom_violin(data = AllPiL, aes(x = pop, y = Pi, color = pop))+
  scale_fill_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000")) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + #this gives whiskers from 8.3% to 91.7% quantile
  geom_pointrange(aes(ymin = lower, ymax = upper)) + #an attempt to get around my missing center point
  scale_color_manual(values = c("#002766", "#004CA3", "#007EF5", "#FFE047", "#F5C000")) +
#theme(axis.text.x=element_text(angle = 90, vjust = -.0001), legend.position="none") 
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), legend.position="none")
dev.off()







	### break to be more statistically rigorous ###
Genhet?









########################
### Figure 2: clusters
# see other script


########################
### Figure 3: Fst ###

snpgdsVCF2GDS(Bur, "BurNuc.gds", method = "biallelic.only")
BurNuc <- snpgdsOpen("BurNuc.gds")
summary(BurNuc)

#load my populations
sample.id <- read.gdsn(index.gdsn(BurNuc, "sample.id"))
BurPop <- read.table("/media/elgon/victoria/NY2020/Bursera/renamedfiles3/popmap")
BurPop <- as.data.frame(BurPop)
BurPop <- BurPop[c(1:9,11,13:26,28:38,41:45),]  #this removes all my filtered samples
names(BurPop) <- c("sample.id", "pop.id")
bpop_code <- BurPop[,2]


#matrix of pairwise population Fst
Burfstmat <- matrix(data = NA, nrow = 5, ncol = 5)
rownames(Burfstmat) <- c("A", "B", "C", "D", "E")
colnames(Burfstmat) <- c("A", "B", "C", "D", "E")

for (j in 1:nrow(Burfstmat)){
  for (i in 1:nrow(Burfstmat)){
    pops <- c("A", "B", "C", "D", "E")
    flag <- bpop_code %in% as.factor(c(pops[j], pops[i]))
    if (pops[j] == pops[i]){
      Burfstmat[i,j] <- NA
    } else {
      samp.sel <- sample.id[flag]
      pop.sel <- as.character(bpop_code[flag])
      Burfstmat[i,j] <- snpgdsFst(BurNuc, sample.id=samp.sel, population = as.factor(pop.sel), method= "W&C84", autosome.only = FALSE)$Fst
    }
  }
}
write.csv(Burfstmat, "Bur_fstBiDD6.30.21.csv")

#old
snpgdsFst(BurNuc, sample.id=pop_code$INDIVIDUALS[c(1:3,10:12,16:20,25:28,35:36)], population = pop_code$POPULATION[c(1:3,10:12,16:20,25:28,35:36)], method= "W&C84", autosome.only = FALSE)$Fst
#0.01947112
#young
snpgdsFst(BurNuc, sample.id=pop_code$INDIVIDUALS[c(4:9,13:15,21:24,29:34,37:40)], population = pop_code$POPULATION[c(4:9,13:15,21:24,29:34,37:40)], method= "W&C84", autosome.only = FALSE)$Fst
#0.005575933
#all
snpgdsFst(BurNuc, sample.id=pop_code$INDIVIDUALS, population = pop_code$POPULATION, method= "W&C84", autosome.only = FALSE)$Fst
#0.01686899

#make a quick PCA:
pca <- snpgdsPCA(BurNuc, autosome.only=FALSE) 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2)) ### checking how much each PC explains, in this case very little
clean_totalRaw <- data.frame(sample.id = pca$sample.id,
                   pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                   EV1 = pca$eigenvect[,1],    # the first eigenvector
                   EV2 = pca$eigenvect[,2],    # the second eigenvector
                   EV3 = pca$eigenvect[,3],    # the third eigenvector
                   EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                   EV5 = pca$eigenvect[,5],    # the fifth eigenvector
                   EV6 = pca$eigenvect[,6],    # the sixth eigenvector
		   EV7 = pca$eigenvect[,7],    # the seventh eigenvector
		   EV8 = pca$eigenvect[,8],    # the eighth eigenvector
stringsAsFactors = FALSE) ### file with indv/pop/eingenvector1/eigenvenctor2
write.csv(clean_totalRaw, file = "Bur_PCA_6.30.21.csv" ) ### saving results

pdf("Bur_PCABiDD_6.30.21.pdf", height = 1.5, width = 1.5)
par(cex = c(.6))
plot(clean_totalRaw$EV1, clean_totalRaw$EV2, xlab="eigenvector 1", ylab="eigenvector 2", pch = popshape$SHAPE, col = popshape$COLOR, bg = popshape$COLOR)
dev.off()
# I think PC1 might have to do with missingness more than anything else
pdf("Bur_PCABiDD_23_6.30.21.pdf", height = 2, width = 2.5)
par(cex = c(.6))
plot(clean_totalRaw$EV2, clean_totalRaw$EV3, xlab="eigenvector 2", ylab="eigenvector 3", pch = popshape$SHAPE, col = popshape$COLOR, bg = popshape$COLOR)
dev.off()


snpgdsClose(Bur)


#plotting fst
mbsites <- c("Maipiero Stop 2", "Maipiero Stop 1", "WCP Bank, 2013, small", "Parque Jaragua camping site", "Parque Jaragua seed traps", "el monte")
locations$colloc <- NA
locations$colloc[locations$Location == "E"] <- "#004CA3"
locations$colloc[locations$Location == "A"] <- "#002766"
locations$colloc[locations$Location == "B"] <- "#FFE047"
locations$colloc[locations$Location == "C"] <- "#007EF5"
locations$colloc[locations$Location == "D"] <- "#F5C000"


# Find a mean for each population
BurPopA <- c(mean(locations$Lon[which(locations$Location == "A")]), mean(locations$Lat[which(locations$Location == "A")]))
BurPopB <- c(mean(locations$Lon[which(locations$Location == "B")]), mean(locations$Lat[which(locations$Location == "B")]))
BurPopC <- c(mean(locations$Lon[which(locations$Location == "C")], na.rm = T), mean(locations$Lat[which(locations$Location == "C")], na.rm = T))
BurPopD <- c(mean(locations$Lon[which(locations$Location == "D")]), mean(locations$Lat[which(locations$Location == "D")]))
BurPopE <- c(mean(locations$Lon[which(locations$Location == "E")]), mean(locations$Lat[which(locations$Location == "E")]))

# Data frame
GeoPops <- rbind(BurPopA, BurPopB, BurPopC, BurPopD, BurPopE) %>% 
  as.data.frame()
colnames(GeoPops) <- c("lon","lat")


# Compute the connection between populations
#might need the code from NuclearAnalysisAllTests.R
interAB <- gcIntermediate(BurPopA,  BurPopB, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interAC <- gcIntermediate(BurPopA,  BurPopC, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interAD <- gcIntermediate(BurPopA,  BurPopD, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interAE <- gcIntermediate(BurPopA,  BurPopE, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interBC <- gcIntermediate(BurPopB,  BurPopC, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interBD <- gcIntermediate(BurPopB,  BurPopD, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interBE <- gcIntermediate(BurPopB,  BurPopE, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interCD <- gcIntermediate(BurPopC,  BurPopD, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interCE <- gcIntermediate(BurPopC,  BurPopE, n=50, addStartEnd=TRUE, breakAtDateLine=F)
interDE <- gcIntermediate(BurPopD,  BurPopE, n=50, addStartEnd=TRUE, breakAtDateLine=F)

poplines <- c(interAB, interAC, interAD, interAE, interBC, interBD, interBE, interCD, interCE, interDE)


#pdf version
pdf("Bur_MapFstBiDDforest_6.30.21.pdf", height = 5.2, width = 5.65) #pdf defaults to size in inches, which is ridiculous. this is pretty good at 5.2 and 5.65, but I am not sure how well it will shrink down.
par(mai=c(.25, .25, .25, .25), cex = 1.2)
plot(allforest, ylim = c(17.7, 17.95), xlim = c(-71.52, -71.28), col = brewer.pal(9, "BuGn"), alpha = .4, legend.lab = "", xlab = "", ylab = "", legend.shrink=.75, axes = F)
mtext(side = 4, line = 4, "% Tree Cover", cex = 1.2)
plot(ocean, col = "paleturquoise1", add = TRUE)

plot(colonies, border = colonies$color, lwd = 3, add = T)
lines(interAB, col="black", lwd=(Burfstmat[1,2] * 100))
lines(interAC, col="black", lwd=(Burfstmat[1,3] * 100))
lines(interAD, col="black", lwd=(Burfstmat[1,4] * 100))
lines(interAE, col="black", lwd=(Burfstmat[1,5] * 100))
lines(interBC, col="black", lwd=(Burfstmat[2,3] * 100))
lines(interBD, col="black", lwd=(Burfstmat[2,4] * 100))
lines(interBE, col="black", lwd=(Burfstmat[2,5] * 100))
lines(interCD, col="black", lwd=(Burfstmat[3,4] * 100))
lines(interCE, col="black", lwd=(Burfstmat[3,5] * 100))
lines(interDE, col="black", lwd=(Burfstmat[4,5] * 100))

points(x=GeoPops$lon, y=GeoPops$lat, col = c("#002766", "#FFE047", "#007EF5", "#F5C000", "#004CA3"), cex=1.75, pch = 19)
points(x=GeoPops$lon, y=GeoPops$lat, col = c("white"), cex=2, pch = 1)
addscalebar(plotunit="latlon", pos = "topright")
box()
dev.off()

#make a same size legend?
pdf("Bur_MapFstBiDDforest_legend.pdf", height = 5.2, width = 5.65) #pdf defaults to size in inches, which is ridiculous. this is pretty good at 5.2 and 5.65, but I am not sure how well it will shrink down.
par(mai=c(.25, .25, .25, .25), cex = 1.2)
plot(allforest, ylim = c(17.7, 17.95), xlim = c(-71.52, -71.28), col = brewer.pal(9, "BuGn"), alpha = 0, legend.lab = "", xlab = "", ylab = "", legend.shrink=.75, axes = F)
legend("topleft", legend=c("Fst of .05", "Fst of .01", "2011 Colony", "2015 Colony", "2016 Colony"), lwd = c(5, 1, 3, 3, 3), col = c("black", "black", "red2", "orangered", "brown"), bty = "n", cex = 2)
dev.off()

png(filename = "Bur_Fst_PrimaryForest.png",res = 300, height = 2000, width = 2200) 
plot(allforest, ylim = c(17.6, 18.0), xlim = c(-71.65, -71.19), col = brewer.pal(9, "BuGn"), legend.lab = "% Tree Cover", axes = F, xlab='',ylab='', legend.line = -2, legend.cex = 2, legend.width=3, legend.shrink=.75, legend.mar=10)
#primary forest has to go first, since it's not a transparent layer
#plot(ocean, col = "cadetblue3", ylim = c(17.7, 17.89), xlim = c(-71.339, -71.211), add = TRUE)
plot(ocean, col = "gray2", ylim = c(17.7, 17.89), xlim = c(-71.339, -71.211), add = TRUE)
points(rbind(BurPopA, BurPopB, MetPopB, BurPopC, BurPopD, BurPopE), xlab = "Longitude", ylab = "Latitude", main = "Sampling Locations", ylim = c(17.5, 18.1), xlim = c(-71.6, -71.15), col =  c("royalblue2", "gray43", "gray43", "deepskyblue", "gray23", "slateblue3"), lwd = 2, cex = 1.5, pch = 19)
box()
legend("bottomright", legend=c("Nests", "No Nests", "Primary Forest"), pch = 21, col = c("royalblue", "black", "darkolivegreen4"), cex = 1.5, pt.bg  = c("dodgerblue2", "black", "green4"))
#the legend doesn't actually add any information needed for this poster
boxpar <- par(fig=c(0,0.4, 0.6, 0.95), new = T)
plot(ocean, col = "cadetblue3", ylim = c(17.8, 19.1), xlim = c(-75, -68))
rect(-71.6, 17.55, -71.15, 18, border = "black", lwd = 2)
box()
dev.off()
#note: I did some finagling in Inkscape to make points larger/text more readable.

######################################################################
 ### Figure 4: age ###

BurEuc2 <- read.csv("C:/Users/torie/Documents/DRgeneflow/FebruaryDraft/Bur_EucDist_noy13_4.16.10.csv", sep = " ")
#I had to drop 22 because it's one of the parent/offspring pair, so these numbers don't align with pop_code exactly
BurEuc2[BurEuc2 == 0] <- NA
BurEuc2 <- as.matrix(BurEuc2)
#young pop A
median(BurEuc2[4:9,4:9], na.rm = T) #106.1181 #mean is 107.666
#old pop A
median(BurEuc2[1:3,1:3], na.rm = T) #103.4366 #mean is 103.6698
#all pop A
median(BurEuc2[1:9,1:9], na.rm = T) #104.5366 #mean is 105.9627
#young pop B
median(BurEuc2[13:15,13:15], na.rm = T) #104.7911 #mean is 105.4736
#old pop B
median(BurEuc2[10:12,10:12], na.rm = T) #104.0115 #mean is 103.5139
#all pop B
median(BurEuc2[10:15,10:15], na.rm = T) #104.665 #mean is 102.9861
#young pop C 
median(BurEuc2[21:23, 21:23], na.rm = T) #101.1322 #mean is 101.7356
#old pop C
median(BurEuc2[16:20,16:20], na.rm = T) #105.6156 #mean is 105.6127
#all pop C
median(BurEuc2[16:23,16:23], na.rm = T) #105.6793 #mean is 105.1034
#young pop D
median(BurEuc2[28:33, 28:33], na.rm = T) #104.6638 #mean is 105.0329
#old pop D
median(BurEuc2[24:27,24:27], na.rm = T) #104.796 #mean is 105.2257
#all pop D
median(BurEuc2[24:33,24:33], na.rm = T) #105.0282 #mean is 104.5824
#young pop E
median(BurEuc2[36:39,36:39], na.rm = T) #105.068 #mean is 105.0272
#old pop E
median(BurEuc2[34:35,34:35], na.rm = T) #123.4912 #mean is 123.4912 ***only two samples here
#all pop E
median(BurEuc2[34:39,34:39], na.rm = T) #113.3017 #mean is 110.5799

#worth noting? for pop B and D the mean distance between young and old is smaller
	#than the mean distance within age classes
	#this is even when I take out the parent-offspring pair

#I tried to use a t-test on the difference between young and old, but these values aren't independent
	#so don't take these test results as meaningful
BurEuc2[upper.tri(BurEuc2, diag=T)] <- NA
 #take out the upper triangle of the matrix
t.test(BurEuc2[4:9,4:9], BurEuc2[1:3,1:3], alternative = "two.sided")
#t = 2.764, df = 15.942, p-value = 0.01386, 95 percent confidence interval: 0.930327 7.061929
	#at population A, young Bursera are more genetically dissimilar than mature Bursera
t.test(BurEuc2[13:15,13:15], BurEuc2[10:12,10:12], alternative = "two.sided")
#t = 1.657, df = 3.8411, p-value = 0.1758, 95 percent confidence interval: -1.37828  5.29771
	#at population B, young and mature Bursera are not more dissimilar
t.test(BurEuc2[21:23, 21:23], BurEuc2[16:20,16:20], alternative = "two.sided")
#t = -1.1092, df = 8.9431, p-value = 0.2963, 95 percent confidence interval: -11.791265   4.037264
	#at population C, young and mature Bursera are not more dissimilar
t.test(BurEuc2[28:33, 28:33], BurEuc2[24:27,24:27], alternative = "two.sided")
#t = -0.18948, df = 9.713, p-value = 0.8536, 95 percent confidence interval: -2.469429  2.083778
	#at population D, young and mature Bursera are not more dissimilar
t.test(BurEuc2[36:39,36:39], BurEuc2[34:35,34:35], alternative = "two.sided")
#not enough observations to compare!


#does colony presence affect genetic distances of Bursera populations?
t.test(BurEuc2[c(1:9,16:23,34:39),c(1:9,16:23,34:39)], BurEuc2[c(10:15,24:33),c(10:15,24:33)], alternative = "two.sided")
#t = 5.6498, df = 355.59, p-value = 3.296e-08, 95 percent confidence interval: 1.810509 3.744019
	#surprisingly, yes. colonies have higher genetic dissimilarity than non-colonies.
	#is this driven by a single outlier? no. the mature individuals at the colony in PJ are both distant from each
	#other and distant from the seedlings.

#what about ALL young trees vs all mature trees? 
t.test(BurEuc2[c(4:9,13:15,21:23,28:33,36:39),c(4:9,13:15,21:23,28:33,36:39)], BurEuc2[c(1:3,10:12,16:20,24:27,34:35),c(1:3,10:12,16:20,24:27,34:35)], alternative = "two.sided")
#t = -3.1026, df = 201.86, p-value = 0.002193, 95 percent confidence interval: -2.766139 -0.616410
	#young trees may be less dissimilar.
	#HOWEVER, this is only because of popE, and if I remove those individuals the difference disappears



	

   ### age differences in dispersion ###
BurEuc2 <- dist(BurG1[-22], method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# can't just read in Bur_EucDist_noy13_6.27.21.csv because I need a ~dist~ object

BurBeta <- betadisper(BurEuc2, c(paste(pop_code$POPULATION[-22], pop_code$AGE[-22], sep = "_"))) #can specify type (spatial median or group centroid), default is spatial median

anova(BurBeta) # groups may be significant
BurpermBeta <- permutest(BurBeta, pairwise = T, permutations = 10000)
	#returns a matrix with p value on bottom diagonal, permuted p-value above
	#observed p-values with pops?
	#A_m vs A_y (0.02) D_m vs D_Y (0.04 actual, 0.03 permut)
BurTukBeta <- TukeyHSD(BurBeta, pairwise = T) # creates set of confidence intervals between mean distance to centroid)
	#the comparisons that matter:
#		dif	lower		upper	p adjusted
#A_y-A_m   9.58113655  -0.6858344 19.8481075 0.0829828
#B_y-B_m   1.13298399 -10.7222929 12.9882609 0.9999989
#C_y-C_m  -8.00045848 -18.6041405  2.6032236 0.2710374
#D_y-D_m   3.36233139  -6.0100879 12.7347507 0.9615670
#E_y-E_m   2.56520894 -10.0092111 15.1396290 0.9993209
	#so despite differences in populations, no significant differences between ages within populations
#and if I do this as ALL young vs ALL old?
BurBetaAge <- betadisper(BurEuc2, c(pop_code$AGE[-22])) #can specify type (spatial median or group centroid), default is spatial median
anova(BurBetaAge) # age class is probably not significant
BurpermBetaAge <- permutest(BurBetaAge, pairwise = T, permutations = 10000)
BurTukBetaAge <- TukeyHSD(BurBetaAge, pairwise = T)
	# sure enough, a permutation test suggests that young trees are not significantly different from mature trees

#finally, visualize pairwise genetic distance difs using bootstrapped tree
set.seed(1005)
strata(BurG1) <- data.frame(other(BurG1))
nameStrata(BurG1) <- ~Age/Population
genind2genpop(pop = ~Age/Population) %>%
aboot(BurEuc2, dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE) #it doesn't look like I can just use Euc? used provesti because it works for SNPs with missing data

#but how to get branch length?

strata(BurG1) <- data.frame(other(BurG1))
nameStrata(BurG1) <- ~Age/Population
BurAgeTree <- BurG1 %>%
genind2genpop(pop = ~data.AGE/data.POPULATION) %>%
  aboot(dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE) #it doesn't look like I can just use Euc? used provesti because it works for SNPs with missing data


#should really do more samples than that, and also use euclidean distance
BurG1 <- loci2genind(BurG) #remake to fix strata mismatch
strata(BurG1) <- as.data.frame(pop_code[,2:3])
#strata(BurG1) <- data.frame(other(BurG1))
nameStrata(BurG1) <- ~POPULATION/AGE
set.seed(1100)
BurAgeTreeP <- BurG1 %>%
  genind2genpop(pop = ~POPULATION/AGE) %>%
  aboot(dist = provesti.dist, sample = 1000, tree = "nj", cutoff = 50, quiet = TRUE) 
pdf("Bur_RootedAgeTreeP_6.30.21.pdf", height = 5, width = 6, pointsize = .5)
plot(BurAgeTreeP, show.node.label = T, direction = "downwards", tip.color = c(rep("#002766",2),rep("#FFE047",2),rep("#007EF5",2),rep("#F5C000",2),rep("#004CA3",2)))
add.scale.bar(cex = 1.5, length=c(0,0.03), lwd = 3.5) #will have to reposition in Inkscape
dev.off() #[match(pop_code$POPULATION, c('A', 'B', 'C', 'D', 'E'))]


set.seed(1120)
BurAgeTree <- BurG1 %>%
genind2genpop(pop = ~POPULATION/AGE) %>%
  aboot(dist = rogers.dist, sample = 1000, tree = "nj", cutoff = 50, quiet = TRUE) 
pdf("Bur_RootedAgeTreeR_6.30.21.pdf", height = 5, width = 6, pointsize = .5)
plot(BurAgeTree, show.node.label = T, direction = "downwards", tip.color = c(rep("#002766",2),rep("#FFE047",2),rep("#007EF5",2),rep("#F5C000",2),rep("#004CA3",2)))
#add.scale.bar()
add.scale.bar(cex = 1.5, length=c(0,0.03), lwd = 3.5) #will have to reposition in Inkscape
dev.off() #[match(pop_code$POPULATION, c('A', 'B', 'C', 'D', 'E'))]


   ###


   ### young vs old nucleotide diversity 
#(in command line)
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/1young.keep --site-pi --out PopAyoung_Pi
awk '{ total += $3 } END { print total/NR }' PopAyoung_Pi.sites.pi #0.232089
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/1old.keep --site-pi --out PopAold_Pi
awk '{ total += $3 } END { print total/NR }' PopAold_Pi.sites.pi #0.241163
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/2young.keep --site-pi --out PopByoung_Pi
awk '{ total += $3 } END { print total/NR }' PopByoung_Pi.sites.pi #0.23437
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/2old.keep --site-pi --out PopBold_Pi
awk '{ total += $3 } END { print total/NR }' PopBold_Pi.sites.pi #0.23285
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/3young.keep --site-pi --out PopCyoung_Pi
awk '{ total += $3 } END { print total/NR }' PopCyoung_Pi.sites.pi #0.234266
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/3old.keep --site-pi --out PopCold_Pi
awk '{ total += $3 } END { print total/NR }' PopCold_Pi.sites.pi #0.229081
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/4young.keep --site-pi --out PopDyoung_Pi
awk '{ total += $3 } END { print total/NR }' PopDyoung_Pi.sites.pi #0.230358
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/4old.keep --site-pi --out PopDold_Pi
awk '{ total += $3 } END { print total/NR }' PopDold_Pi.sites.pi #0.227984
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/5young.keep --site-pi --out PopEyoung_Pi
awk '{ total += $3 } END { print total/NR }' PopEyoung_Pi.sites.pi #0.235365
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/5old.keep --site-pi --out PopEold_Pi
awk '{ total += $3 } END { print total/NR }' PopEold_Pi.sites.pi #0.186056

#all together
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/Allyoung.keep --site-pi --out Allyoung_Pi
awk '{ total += $3 } END { print total/NR }' Allyoung_Pi.sites.pi #0.233632
vcftools --vcf m33c90h01BiDD.vcf --keep /media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/Allold.keep --site-pi --out Allold_Pi
awk '{ total += $3 } END { print total/NR }' Allold_Pi.sites.pi #0.232912
#nucleotide diversity is pretty similar in young and old

#getting significance via variance: 
pegas::nuc.div(BurBin[1:6,], variance = T) #0.23690694 0.01871747
pegas::nuc.div(BurBin[7:18,], variance = T) #0.22814094 0.01393957
 #old then young of PopB
pegas::nuc.div(BurBin[19:24,], variance = T) #0.23474459 0.01837743
pegas::nuc.div(BurBin[25:30,], variance = T) #0.22958987 0.01757939
 #old then young of PopC
pegas::nuc.div(BurBin[31:40,], variance = T) #0.22333198 0.01392389
pegas::nuc.div(BurBin[41:48,], variance = T) #0.22975062 0.01571809
 #old then youn of PopD
pegas::nuc.div(BurBin[49:56,], variance = T) #0.22382319 0.01491772
pegas::nuc.div(BurBin[57:68,], variance = T) #0.22622669 0.01370669
 #old then young of PopE
pegas::nuc.div(BurBin[69:72,], variance = T) #0.19699104 0.01653733
pegas::nuc.div(BurBin[73:80,], variance = T) #0.23048958 0.01581934


 #violin plots (figure 4b)
Ayoung <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopAyoung_Pi.sites.pi", sep = "")
Byoung <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopByoung_Pi.sites.pi", sep = "")
Cyoung <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopCyoung_Pi.sites.pi", sep = "")
Dyoung <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopDyoung_Pi.sites.pi", sep = "")
Eyoung <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopEyoung_Pi.sites.pi", sep = "")
Aold <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopAold_Pi.sites.pi", sep = "")
Bold <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopBold_Pi.sites.pi", sep = "")
Cold <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopCold_Pi.sites.pi", sep = "")
Dold <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopDold_Pi.sites.pi", sep = "")
Eold <- read.csv("/media/elgon/victoria/NY2020/Bursera/filtertests/nucreads/missing33/PopEold_Pi.sites.pi", sep = "")

AgePi <- cbind(Ayoung[,1:3], Aold[,3], Byoung[,3], Bold[,3], 
	Cyoung[,3], Cold[,3], Dyoung[,3], Dold[,3],
	Eyoung[,3], Eold[,3])
colnames(AgePi) <- c("Chrom", "Pos", "A_young", "A_old", "B_young", "B_old",
	"C_young", "C_old", "D_young", "D_old", "E_young", "E_old")
summary(AgePi) #less similar across populations and age classes than metopium


library(tidyr)
library(ggplot2)
OldPi <- cbind(Aold[,1:3], Bold[,3], Cold[,3], Dold[,3], Eold[,3])
colnames(OldPi)[3:7] <- c("2011 Colony", "Parque Jaragua Non-Colony", "Mapioró Colony", "Mapioró Non-Colony", "Parque Jaragua Colony")
YoungPi <- cbind(Ayoung[,1:3], Byoung[,3], Cyoung[,3], Dyoung[,3], Eyoung[,3])
colnames(YoungPi)[3:7] <- c("2011 Colony", "Parque Jaragua Non-Colony", "Mapioró Colony", "Mapioró Non-Colony", "Parque Jaragua Colony")

OldPi <- as.data.frame(pivot_longer(OldPi, names_to = "Pop", values_to = "Pi", -c("CHROM", "POS")))
OldPi$Age <- "Old"
YoungPi <- as.data.frame(pivot_longer(YoungPi, names_to = "Pop", values_to = "Pi", -c("CHROM", "POS")))
YoungPi$Age <- "Young"
LongAgePi <- rbind(OldPi, YoungPi)
LongAgePi$Pop <- as.factor(LongAgePi$Pop)
LongAgePi$Pop <- factor(LongAgePi$Pop, levels = c("2011 Colony", "Parque Jaragua Colony", "Mapioró Colony", "Parque Jaragua Non-Colony", "Mapioró Non-Colony"))

LongAgeMean <- LongAgePi %>% group_by(Pop, Age) %>% mutate(AgeMean = mean(Pi))
LongAgeMedian <- LongAgePi %>% group_by(Pop, Age) %>% mutate(AgeMedian = median(Pi))

#I could change the levels encoded by the age*population classes so those show up underneath the correct branch on the tree
	#instead I might just plot them in the same order as my other plots


pdf(file = "Bur_AgePiViolin_6.30.21.pdf", width = 6, heigh = 4)
#violin plots (will adjust colors in inkscape, because it's hard to get both age AND population)
LongAgePi %>%
ggplot(aes(fill=Pop, y=Pi, x=Pop)) + 
    geom_violin(aes(x = Pop, y = Pi, color = Age), position="dodge") +
    scale_color_manual(values = c("#002766","#004CA3","#007EF5", "#FFE047", "#F5C000")) +    
    scale_fill_manual(values = c("#002766", "#004CA3","#007EF5", "#FFE047","#F5C000")) +    
    #stat_summary(fun.y = mean, geom = "errorbar", aes(group = Pop), width = 1, linetype = "solid") +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
	 geom = "crossbar", width = 1, linetype = "solid", aes(x = Pop, y = Pi, group=Age), position=position_dodge(.9)) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 1, 
	 linetype = "dotted", aes(x = Pop, y = Pi, group=Age), position=position_dodge(.9)) +
    theme_bw()  +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="none")
dev.off()



###########################################################################
  ### Figure 5: chloroplast ###
### figure 5a: map of Fst
### figure 5b: chloroplast map
### figure 5c: genetic by geographic distance

# None of these really matter, since I didn't have enough SNPs to consider this worth it.