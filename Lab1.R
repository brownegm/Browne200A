#####Marvin Browne#####
#####Lab Exercise 1: Introduction to estimating diversification####

####Lineage through time plots####

##Load phytools & dependencies and setwd()
library(phytools)
setwd( "C:/Users/Marvin/Documents/UCLA/Fall 2017/200A")

##Load a tree from the file "etheostoma_percina_chrono.tre".
darter.tree<-read.tree("etheostoma_percina_chrono.tre")#Tree is from Near et al. (2011). They sampled 201 of 216 species from the group. We assume 100% sampling. 

##Plot tree
plotTree(darter.tree,ftype="i", fsize=0.4, type="fan", lwd=1)

##use the function ltt in the phytools package to generate a LTT plot
obj<-ltt(darter.tree,log.lineages = FALSE)

##Why are we getting NA for gamma? Our tree is not fully bifurcating
##Fixing that....
darter.tree

## A fully bifurcating tree with N tips will have N-1 internal nodes. 
is.binary(darter.tree)#false= not bifurcating 
##***darter.tree has 201 tips and 198 internal nodes so that needs to be fixed****

##bifurcate tree
darter.tree<-multi2di(darter.tree)
darter.tree

## test for bifurcation
is.binary(darter.tree)#TRUE>>> We have randomly resolved the interanl nodes that were multifurcating

##Create LTT plot from new bifurcating tree
obj<-ltt(darter.tree,plot = FALSE)
obj# gamma statistic is no longer NA for the tree!!! 

##Generate plot of above 
plot(obj, log.lineages=FALSE, main="LTT Plot for Darters")
#plot the tree and the LTT plot together
plot(obj,show.tree=TRUE,log.lineages=FALSE,main="LTT Plot for Darters")

##replotting the above graph
plot(obj,log.lineages=FALSE,log="y",main="LTT plot for darters",
     ylim=c(2,Ntip(darter.tree)))
## we can overlay the pure-birth prediction:
h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
lines(x,2*exp(b*x),col="red",lty="dashed",lwd=2)

##We can compare the observe LTT, to simulated LTTs assuming a pure-birth process of the same duration& resulting the same total number 
## of species. We can do this using the phytools function pbtree:
trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",
              quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT of darters compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(darter.tree))),lty="dashed",lwd=2,col="red")
## now let's overlay our original tree
ltt(darter.tree,add=TRUE,lwd=2)

## Alternatively, the function ltt95 gives a (1-a)% CI for the LTT based on a set of trees. We use it here using a simulated distribution of trees, 
## but this might also be useful if we have a set of trees from the posteriod distribution of a Bayesian analysis
ltt95(trees,log=TRUE)
title(main="LTT of darters compared to simulated LTTs")
ltt(darter.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2)

####Fitting pure birth and birth-death models to trees####

##birthdeath() reports compound parameters, d/b(extinction fraction) and b-d(net diversification).
library(phytools)
fitbd<-birthdeath(darter.tree)
fitbd

##To get the b and d rates separeately we will use the phytools function, bd()
bd(fitbd)
## One problem with this function is that the birth and death rates are conditioned on full sampling of the tree.
##THis is not the case the darter tree but we will not pursue more appopriate rate estimation here.
####The y-statistic####

## The gamma statistic of P&H 2010 was designed to have a standard normal distribution for trees generated under a pure-birth speciation model. 
#Let's test this with our 100 simulated pure-birth phylogenies
g<-sapply(trees,function(x) ltt(x,plot=FALSE)$gamma)
hist(g,main=expression(paste("Distribution of ",gamma," from simulation")))
mean(g)
var(g)

##We can test hypotheses about y. This is done automatically with ltt:
obj<-ltt(darter.tree,plot=FALSE)
print(obj)

##Let's simulate a tree under a different model of lineage accumulation- the coalescent- and see what the result is.
##Theoretically, this should result in a significantly positive gamma.
coal.tree<-rcoal(n=100)
plotTree(coal.tree,ftype="off")
obj<-ltt(coal.tree,log.lineages=FALSE,log="y")

darter.gamma<-obj$gamma #ltt returns a gamma value as one of its elements

## Let's compare the mean gamma value of 200 pb trees to gamma from the coalescent tree. 
trees<-pbtree(n=100,nsim=200,scale=max(nodeHeights(coal.tree)))
ltt95(trees,log=TRUE)
title(main="Simulated coalescent trees compared to pure-birth LTTs")
ltt(coal.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2,lty="dashed")

####Incomplete Sampling####
###The gamma value we recieved before assumed that the entire clade was sampled. So to test the 
###effect of incomplete sampling we can simulate a tree with the expect amount of branches and trim the 
### the tree to what we have and see if we get similar results. This test is called the MCCR test.

##performs the MCCR test;we use a pure birth estimate of the tree based upon their total age and richness. 
library(geiger)
age <- 25.91862
richness <- 216
darterbirth =  (log(richness) - log(2))/age
darterbirth

##This simulates gamma values when trees are undersampled.
##we will grow trees with n=34 and prune them down to 13 taxa

missing<- 15 

num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations)#g1_null will hold the simulated gamma values
for(i in 1:num_simulations){
  sim_tree<-sim.bdtree(darterbirth, d=0, stop = "taxa", n=richness)
  prune<-drop.random(sim_tree, missing)#prune down the tree to the number of taxa in the phylogeny
  g1_null[i]<-gammaStat(prune)
}#for 200 simulations simulate a tree using values of darterbirth, then take out/prune tree by the amount of taxing in the missing vector(15) to get to the number of taxa in the phylogeny
#then take the gamma stat of the pruned tree and store it in the vector named g1_null at the end of each simulation.

##Create a histogram of the null distribution(expected values derived from for loop)
hist(g1_null)

#arrow indicates where the observed gamma falls in the null dist just created
arrows(darter.gamma, 40, darter.gamma, 0, col="red", lwd=2)

## Which of the null values are smaller(more negative) than the data?
smallerNull<-g1_null<=darter.gamma

##How many TRUEs are there?
count<-sum(smallerNull)

##finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval# p-value=1 so there is exactly 0 difference between the observed and the simulations


####Exercises####
##1. Calculate the gamma statistic for this phylogeny of homalopsid snakes from Alfaro et al., 2008

#load packages and set working directory
library(phytools)
setwd( "C:/Users/Marvin/Documents/UCLA/Fall 2017/200A")

#Load a tree from the file "homalops.phy".
snake.tree<-read.tree("homalops.phy")#Tree from Alfaro et al., 2008

#perform ltt on snake tree to get gamma statistic
obj<-ltt(snake.tree,log.lineages = FALSE)
obj#Gamma statistic for this tree is -3.2411 with a p-value of 0.0012

## 2. Given this gamma value, what would you conclude about the tempo of speciation in this clade?
#I would conclude that the majority of the speciation has occured further back in the tree. The lineages
# through time plot illustrates that the lineages increased rapidly at first and became more steady closer to the 
#present. 

## 3. Given that the crown age of the snake radiation is 22MY and the total richness of the clade is 34 species, determine 
##whether the observed gamma could be due to the amount of incomplete sampling in the empirical tree. On the basis of the MCCR
## test what can you conclude about the tempo of homalopsid snake diversification?

library(geiger)
age <- 22
richness <- 34
snakebirth =  (log(richness) - log(2))/age
snakebirth

#This simulates gamma values when trees are undersampled.

missing<- 13 

num_simulations<-21 #number of simulations
g1_null<-numeric(num_simulations)#g1_null will hold the simulated gamma values
for(i in 1:num_simulations){
  sim_tree<-sim.bdtree(snakebirth, d=0, stop = "taxa", n=richness)
  prune<-drop.random(sim_tree, missing)#prune down the tree to the number of taxa in the phylogeny
  g1_null[i]<-gammaStat(prune)
}

snake.gamma <- obj$gamma
smallerNull<-g1_null<=snake.gamma

#How many TRUEs are there?
count<-sum(smallerNull)

#finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval#P-value is 0.004975124

#The tempo is significant different between the simulated values and the observed values. This may indicate
#that the diversification of lineages may actually start off slowly and increase significantly towards the ends of the tree

## 4. What is the birth rate and death rate of the homalopsid tree. 
fitbd<-birthdeath(snake.tree)
fitbd

bd(fitbd)#birth rate=0.1005303  death rate=0.0000000 
## 5. Find a time-calibrated phylogeny for a group that interests you(ideally with more than 30 tips and fewer than 200). Do the following:
##(1) Describe the clade(including a description of the number of tips in the tree and the total number of species in the clades) and provide
## a reference or citation to the source. (2) Fit a birthdeath model to this tree and report b and d.(3) Perform and MCCR test and describe whether 
## the gamma value is extreme or not given the level of sampliing in the tree. 

####Data is from:Fine PVA, Zapata F , and Daly DC. (2014) Investigating processes of Neotropical rain forest tree diversification by examining the####
####evolution and historical biogeography of the Protieae (Burseraceae). Evolution 68: 1988-2004. doi:10.1111/evo.12414.####
tropical.tree<-read.nexus("AllNormal_MCCT.tre")

outgroups <- c("Bo_sacra", "Bu_simaruba", "Ca_pilosum", "Co_edulis", "Pi_mexicana", "Rh_trilobata", "Sa_griffithii")

# Modify trees to work only with ingroup

tropical.tree <- drop.tip( tropical.tree, outgroups )

plotTree(tropical.tree,ftype="i", fsize=0.4, type="fan", lwd=1)

##(1) The clade has 111 tips with 110 internal nodes. Seven outgroups were removed by the authors as outgroups in the data.I removed the same groups from the data set here. 
##If the those species were not removed there would be 118 species total in the tree. Furthermore, the authors further trim the tree to 102(73%) sampled for their analyses but I did not. 
##The tree describes the evolutionary trajectory of Neo tropical trees within a tribe of the family Burseraceae. 
##
##(2) The birth rate in this group is 0.1228332 and the death rate is 0.00000000  
birthdeath<-birthdeath(tropical.tree)
birthdeath

bd(birthdeath)
## (3)
age <- 55
richness <- 111
tropicalbirth =  (log(richness) - log(2))/age
tropicalbirth

#This simulates gamma values when trees are undersampled.

missing<- 7 

num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations)#g1_null will hold the simulated gamma values
for(i in 1:num_simulations){
  sim_tree<-sim.bdtree(tropicalbirth, d=0, stop = "taxa", n=richness)
  prune<-drop.random(sim_tree, missing)#prune down the tree to the number of taxa in the phylogeny
  g1_null[i]<-gammaStat(prune)
}

tropical.gamma <- obj$gamma
smallerNull<-g1_null<=tropical.gamma

#How many TRUEs are there?
count<-sum(smallerNull)

#finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval#p-value=0.004975124

##This result does indicate that the clade was incompletely sampled. However, the authors of the paper did a good job of indicating that their results were only from a portion of the possible species within the clade

