#####################
# Phylogeny
#####################

# R markdownとapeの相性が悪い？ので、R scriptで分析をすすめる。

# install & load packages
rm(list=ls())  # reset workspace
library(tidyverse)  # dataset modification
library(ape)        # Phylogeny
library(phytools)   # Phylogeny
setwd("/media/sf_Dropbox/Meta_F1HybridVariation/Data")

##############################################################
# Tree for dataset "Recip lnRR"
##############################################################

##### extracting names of species #####
species <- read.table("../Analysis/sexual.ES.lnRR.recip.csv", head = TRUE, sep = "\t") %>%
  select(spL.name) %>% # selecting column with species names
  unique(.) # extract unique species names
species.recip <- as.vector(species$spL.name) 
write.table(species, "../Data/species.recip.lnRR.names.txt", 
            col.names = F, row.names = F, quote = F)

############################################################################
##### Linux command #####
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'
sed -i "s/_/ /g" species.recip.lnRR.names.txt


##### Downloading taxonomical tree from BLAST taxonomy #####
@ https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi?session=1Za3H4MydSSQb2OL02M43p_XUlMnWWaeFmTB5ADsumEgV42K-Exx4OhN8EWTxYM-uj2pvh_GW67DmKHDrBMeSk_xOvrQ4N7krqsLzm-9nH5EEVp30LOzNg1WbXITHmvG2_4QL8&en=588766#588766

##### Linux command #####
# remove all "'"" & replace spaces " " with underscores "_"  
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i -e "s/\"//g" -e "s/'//g" -e "s/\ /_/g" PhylTree.spL.recip.lnRR.rawBLAST.phy  
   #  " を置換する前には、\  
   #  ' では、s~g全体を””で囲む  
   #  スペース　では、\  
# change species name
sed -i -e "s/Coturnix_japonica/Coturnix_coturnix_japonica/g" PhylTree.spL.recip.lnRR.rawBLAST.phy
sed -i -e "s/Dryophytes/Hyla/g" PhylTree.spL.recip.lnRR.rawBLAST.phy
############################################################################

##### Load original tree on R #####
tree.spL.recip.lnRR <- read.tree("../Data/PhylTree.spL.recip.lnRR.rawBLAST.phy")
par(oma=c(0.1,0.1,0.1,0.1), ps=14)
is.binary.tree(tree.spL.recip.lnRR) # 二値的な系統樹になっているか確認
tree.spL.recip.lnRR[["edge.length"]] <- NULL #remove edge lengths!
tree.spL.recip.lnRR <- collapse.singles(tree.spL.recip.lnRR) #colapse singletons
plot(tree.spL.recip.lnRR)

##### Add missing species #####
# add "Chorthippus_parallelus_eythropus" as a sister species of "Chorthippus_parallelus_parallelus"
tree.spL.recip.lnRR.add <- bind.tip(tree.spL.recip.lnRR, "Chorthippus_parallelus_eythropus", 
                         where = which(tree.spL.recip.lnRR$tip=="Chorthippus_parallelus_parallelus"))
plot(tree.spL.recip.lnRR.add)


##### Randomly solve Polytomy #####
set.seed(2014) #making it replicable
# set.seed()は乱数種を指定する関数で、常に同じ乱数を発生させられる
tree.spL.recip.lnRR.random <- multi2di(tree.spL.recip.lnRR.add,random=TRUE)
is.binary.tree(tree.spL.recip.lnRR.random)
plot(tree.spL.recip.lnRR.random)
write.tree(tree.spL.recip.lnRR.random, file="Phyltree.spL.recip.lnRR.random.tre")
png('Phyltree.spL.recip.lnRR.random.png')        #  デバイスドライバの用意．最初の引数に作製するファイル名を与える
plot(tree.spL.recip.lnRR.random)   #  グラフを描く
dev.off()              #  デバイスを閉じる


##############################################################
# Tree for dataset "Recip SMD"
##############################################################

##### extracting names of species #####
species <- read.table("../Analysis/sexual.ES.recip.csv", head = TRUE, sep = "\t") %>%
  select(spL.name) %>% # selecting column with species names
  unique(.) # extract unique species names
species.recip.SMD <- as.vector(species$spL.name) 
write.table(species.recip.SMD, "../Data/species.recip.SMD.names.txt", 
            col.names = F, row.names = F, quote = F)

############################################################################
##### Linux command #####
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i "s/_/ /g" species.recip.SMD.names.txt


##### Downloading taxonomical tree from BLAST taxonomy #####
@ https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi?session=1Za3H4MydSSQb2OL02M43p_XUlMnWWaeFmTB5ADsumEgV42K-Exx4OhN8EWTxYM-uj2pvh_GW67DmKHDrBMeSk_xOvrQ4N7krqsLzm-9nH5EEVp30LOzNg1WbXITHmvG2_4QL8&en=588766#588766

##### Linux command #####
# remove all "'"" & replace spaces " " with underscores "_"  
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i -e "s/\"//g" -e "s/'//g" -e "s/\ /_/g" PhylTree.spL.recip.SMD.rawBLAST.phy  
   #  " を置換する前には、\  
   #  ' では、s~g全体を””で囲む  
   #  スペース　では、\  
# change species name
sed -i -e "s/Coturnix_japonica/Coturnix_coturnix_japonica/g" PhylTree.spL.recip.SMD.rawBLAST.phy
sed -i -e "s/Dryophytes/Hyla/g" PhylTree.spL.recip.lnRR.rawBLAST.phy
############################################################################

##### Load original tree on R #####
tree.spL.recip.SMD <- read.tree("../Data/Phyltree.spL.recip.SMD.rawBLAST.phy")
par(oma=c(0.1,0.1,0.1,0.1), ps=14)
is.binary.tree(tree.spL.recip.SMD) # 二値的な系統樹になっているか確認
tree.spL.recip.SMD[["edge.length"]] <- NULL #remove edge lengths!
tree.spL.recip.SMD <- collapse.singles(tree.spL.recip.SMD) #colapse singletons
plot(tree.spL.recip.SMD)

##### Add missing species #####
# add "Chorthippus_parallelus_eythropus" as a sister species of "Chorthippus_parallelus_parallelus"
tree.spL.recip.SMD.add <- bind.tip(tree.spL.recip.SMD, "Chorthippus_parallelus_eythropus", 
                               where = which(tree.spL.recip.SMD$tip=="Chorthippus_parallelus_parallelus"))
# add "Gryllus_armatus" as a sister species of "Gryllus_rubens"
tree.spL.recip.SMD.add <- bind.tip(tree.spL.recip.SMD.add, "Gryllus_armatus", where = which(tree.spL.recip.SMD.add$tip=="Gryllus_rubens"))
plot(tree.spL.recip.SMD.add)


##### Randomly solve Polytomy #####
set.seed(2014) #making it replicable
# set.seed()は乱数種を指定する関数で、常に同じ乱数を発生させられる
tree.spL.recip.SMD.random <- multi2di(tree.spL.recip.SMD.add,random=TRUE)
is.binary.tree(tree.spL.recip.SMD.random)
plot(tree.spL.recip.SMD.random)
write.tree(tree.spL.recip.SMD.random, file="Phyltree.spL.recip.SMD.random.tre")
png('Phyltree.spL.recip.SMD.random.png')        #  デバイスドライバの用意．最初の引数に作製するファイル名を与える
plot(tree.spL.recip.SMD.random)   #  グラフを描く
dev.off()              #  デバイスを閉じる



##############################################################
# Tree for dataset "Epistasis"
##############################################################

# extracting names of species
species <- read.table("../Analysis/sexual.ES.epistasis.csv", head = TRUE, sep = ",") %>%
  select(spL.name) %>% # selecting column with species names
  unique(.) # extract unique species names
species.epistasis <- as.vector(species$spL.name) 
write.table(species.epistasis, "../Data/species.epistasis.names.txt", 
            col.names = F, row.names = F, quote = F)

#--Linux command--
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i "s/_/ /g" species.epistasis.names.txt

# remove all "'"" & replace spaces " " with underscores "_"  
# ~Linux command~  
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i -e "s/\"//g" -e "s/'//g" -e "s/\ /_/g" PhylTree.spL.epistasis.rawBLAST.phy  
# change species name
sed -i -e "s/Coturnix_japonica/Coturnix_coturnix_japonica/g" PhylTree.spL.epistasis.rawBLAST.phy
sed -i -e "s/Dryophytes/Hyla/g" PhylTree.spL.epistasis.rawBLAST.phy

#--------Load original tree----------

tree.spL.epistasis <- read.tree("../Data/PhylTree.spL.epistasis.rawBLAST.phy")
par(oma=c(0.1,0.1,0.1,0.1), ps=14)
is.binary.tree(tree.spL.epistasis) # 二値的な系統樹になっているか確認
tree.spL.epistasis[["edge.length"]] <- NULL #remove edge lengths!
tree.spL.epistasis <- collapse.singles(tree.spL.epistasis) #colapse singletons
plot(tree.spL.epistasis)

#--------Add missing species-------------

# add "Gryllus_armatus" as a sister species of "Gryllus_rubens"
tree.spL.epistasis.add <- bind.tip(tree.spL.epistasis, "Gryllus_armatus", where = which(tree.spL.epistasis$tip=="Gryllus_rubens"))
plot(tree.spL.epistasis.add)

# add "Chorthippus_parallelus_eythropus" as a sister species of "Chorthippus_parallelus_parallelus"
tree.spL.epistasis.add2 <- bind.tip(tree.spL.epistasis.add, "Chorthippus_parallelus_eythropus", 
                                where = which(tree.spL.epistasis.add$tip=="Chorthippus_parallelus_parallelus"))
plot(tree.spL.epistasis.add2)

#--------Randomly solve Polytomy-------------

set.seed(2014) #making it replicable
# set.seed()は乱数種を指定する関数で、常に同じ乱数を発生させられる
tree.spL.epistasis.random <- multi2di(tree.spL.epistasis.add2,random=TRUE)
is.binary.tree(tree.spL.epistasis.random)
plot(tree.spL.epistasis.random)
write.tree(tree.spL.epistasis.random, file="PhylTree.spL.epistasis.random.tre")
png('PhylTree.spL.epistasis.random.png')        #  デバイスドライバの用意．最初の引数に作製するファイル名を与える
plot(tree.spL.epistasis.random)   #  グラフを描く
dev.off()              #  デバイスを閉じる




##############################################################
# Tree for dataset "Paternal/Maternal"
##############################################################

##### extracting names of species #####
species <- read.table("../Analysis/sexual.ES.matpat.2.csv", head = TRUE, sep = "\t") %>%
  select(spL.name) %>% # selecting column with species names
  unique(.) # extract unique species names
species.matpat <- as.vector(species$spL.name) 
write.table(species.matpat, "../Data/species.matpat.names.txt", 
            col.names = F, row.names = F, quote = F)


##### Linux command #####
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -e "s/_/ /g" species.matpat.names.txt > species.matpat.names2.txt


##### Downloading taxonomical tree from BLAST taxonomy #####
https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi?session=1Za3H4MydSSQb2OL02M43p_XUlMnWWaeFmTB5ADsumEgV42K-Exx4OhN8EWTxYM-uj2pvh_GW67DmKHDrBMeSk_xOvrQ4N7krqsLzm-9nH5EEVp30LOzNg1WbXITHmvG2_4QL8&en=588766#588766

##### Linux command #####
# remove all "'"" & replace spaces " " with underscores "_"  
# ~Linux command~  
cd '/media/sf_Dropbox/Meta_F1HybridVariation/Data'  
sed -i -e "s/\"//g" -e "s/'//g" -e "s/\ /_/g" PhylTree.spL.matpat.rawBLAST.phy  
# change species name
sed -i -e "s/Coturnix_japonica/Coturnix_coturnix_japonica/g" PhylTree.spL.matpat.rawBLAST.phy


##### Load original tree on R #####
tree.spL.matpat <- read.tree("../Data/PhylTree.spL.matpat.rawBLAST.phy")
par(oma=c(0.1,0.1,0.1,0.1), ps=14)
is.binary.tree(tree.spL.matpat) # 二値的な系統樹になっているか確認
tree.spL.matpat[["edge.length"]] <- NULL #remove edge lengths!
tree.spL.matpat <- collapse.singles(tree.spL.matpat) #colapse singletons
plot(tree.spL.matpat)

##### Add missing species #####
# add "Gryllus_armatus" as a sister species of "Gryllus_rubens"
# tree.add <- bind.tip(tree, "Gryllus_armatus", where = which(tree$tip=="Gryllus_rubens"))
# plot(tree.add)
# add "Chorthippus_parallelus_eythropus" as a sister species of "Chorthippus_parallelus_parallelus"
tree.spL.matpat.add <- bind.tip(tree.spL.matpat, "Chorthippus_parallelus_eythropus", 
                                where = which(tree.spL.matpat$tip=="Chorthippus_parallelus_parallelus"))
plot(tree.spL.matpat.add)


##### Randomly solve Polytomy #####
set.seed(2014) #making it replicable
# set.seed()は乱数種を指定する関数で、常に同じ乱数を発生させられる
tree.spL.matpat.random <- multi2di(tree.spL.matpat.add,random=TRUE)
is.binary.tree(tree.spL.matpat.random)
plot(tree.spL.matpat.random)
write.tree(tree.spL.matpat.random, file="PhylTree.spL.matpat.random.tre")
png('PhylTree.spL.matpat.random.png')        #  デバイスドライバの用意．最初の引数に作製するファイル名を与える
plot(tree.spL.matpat.random)   #  グラフを描く
dev.off()              #  デバイスを閉じる



