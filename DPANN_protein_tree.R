library(readxl)
library(xlsx)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(treeio)
library(cowplot)
library(ape)
library(phylotools)

#DPANN nosZ, cytochrome C, and pmo protein tree

setwd("/Users/irene/Desktop/DPANN/") 

DPANN_protein_MAD <- read.tree("new_DPANN_protein_tree/all_pmo_Cox_nosZ_new_MAFFT_trimal_MAD.treefile")
DPANN_protein_MAD %>% ggtree()

DPANN_protein_tibble <- as_tibble(DPANN_protein_MAD)

#trying to incorporate bootstraps
boot <- tibble(
  # hard-code node ID: internal nodes start after tip nodes,
  # and phy$node.label is in the same order as internal nodes
  node=1:Nnode(DPANN_protein_MAD) + Ntip(DPANN_protein_MAD), 
  # Don't print BS < 50%
  # (or to show all BS, just do `bootstrap = phy$node.label`)
  bootstrap = DPANN_protein_MAD$node.label)

boot_split<-separate(data = boot, col = bootstrap, into = c("SHalrt", "UFboot"), sep = "\\/")

DPANN_merge_boot<-merge(DPANN_protein_tibble,boot,by="node", all.x=TRUE)
DPANN_merge_boot_split<-merge(DPANN_protein_tibble,boot_split,by="node", all.x=TRUE)


setwd("/Users/irene/Desktop/DPANN/Excel_files/")
#write.xlsx(DPANN_merge_boot, file = "protein_tree_annotation.xlsx", sheetName="DPANN_boot", append=TRUE)

DPANN_annote <- read_excel("protein_tree_annotation.xlsx", 
                           sheet="annotation")

protein_merge<-left_join(DPANN_protein_tibble,DPANN_annote,by="label")
protein_merge_boot<-left_join(DPANN_merge_boot_split,DPANN_annote,by="label")
#write.xlsx(protein_merge, file = "protein_tree_annotation.xlsx", sheetName="tree_merge", append=TRUE)

DPANN_clean <- read_excel("protein_tree_annotation.xlsx", sheet="tree_merge",col_types = c("numeric","numeric","numeric",
                                                                                           "text","numeric","numeric","text"))
test_merge<- left_join(protein_merge,boot,by="node")
DPANN_treedata_test<-as.treedata(test_merge) #test this out 
DPANN_treedata<-as.treedata(protein_merge)
DPANN_treedata2<-as.treedata(protein_merge_boot)
#DPANN_treedata2<-as.treedata(DPANN_clean)


DPANN_treedata %>% ggtree()
DPANN_treedata2 %>% ggtree() ##finally got this to work!!
#DPANN_treedata_test %>% ggtree()

#tree with bootstraps and SHalrt
DPANN_tree_boots<-DPANN_treedata2 %>% ggtree() + 
  geom_tippoint(aes(fill=annotation),shape=21,size=1,stroke=NA)+ 
  scale_fill_brewer(palette='Set2',direction=-1)+
  theme(legend.position = "right")+
  geom_point2(aes(subset = !isTip,color  = as.numeric(UFboot)),size=2,shape=18)+
  scale_color_continuous(name='UFboot', limits=c(0, 100),
                         oob=scales::squish,
                         low='white', high='black') +
  geom_text(aes(label=SHalrt), hjust=-.25, size = 2)+
  ggtitle('DPANN protein tree')

DPANN_protein_tree<-DPANN_treedata2 %>% ggtree() + 
  geom_tippoint(aes(fill=annotation),shape=21,size=1,stroke=NA)+ 
  scale_fill_brewer(palette='Set2',direction=-1)+
  theme(legend.position = "right")+
  geom_text(aes(label=node), hjust=-.25, size = 2)+
  ggtitle('DPANN protein tree')

rerooted_treedata<-root(DPANN_treedata2, node=1871)

rerooted_tree<-rerooted_treedata2 %>% ggtree() + 
  geom_tippoint(aes(fill=annotation),shape=21,size=2,stroke=NA)+ 
  scale_fill_brewer(palette='Set2',direction=-1)+
  theme(legend.position = "right")+
  ggtitle('DPANN protein tree')

rerooted_tree

rerooted_treedata2<-root(DPANN_treedata2, node=3475)
rerooted_tree_2<-rerooted_treedata2 %>% ggtree() + 
  geom_tippoint(aes(fill=annotation),shape=21,size=3,stroke=NA)+ 
  scale_fill_brewer(palette='Set2',direction=-1)+
  geom_point2(aes(subset = !isTip,color  = as.numeric(UFboot)),size=2,shape=18)+
  scale_color_continuous(name='UFboot', limits=c(0, 100),
                         oob=scales::squish,
                         low='white', high='black') +
  geom_text(aes(label=SHalrt), hjust=-.25, size = 2.5)+
  theme(legend.position = "right")+
  ggtitle('DPANN protein tree')

rerooted_tree_nodes<-rerooted_treedata2 %>% ggtree() + 
  geom_tippoint(aes(fill=annotation),shape=21,size=2,stroke=NA)+ 
  scale_fill_brewer(palette='Set2',direction=-1)+
  geom_text(aes(label=node), hjust=-.25, size = 2)+
  theme(legend.position = "right")+
  ggtitle('DPANN protein tree')

rerooted_tree_2

viewClade(rerooted_tree, node=2888)

scaled_tree<-scaleClade(rerooted_tree_2, 3475, 0.1)%>% ##cytochrome c oxidase
  scaleClade(3495, 0.1)%>% ##DPANN
  scaleClade(3581, .01)%>% ##Sec nitrous oxide
  scaleClade(3553, .01)%>% ##Sec nitrous oxide
  scaleClade(1854, .01)%>% ##Sec nitrous oxide
  scaleClade(3434, .01)%>% ##Sec nitrous oxide
  scaleClade(3389, .01)%>% ##Sec nitrous oxide
  scaleClade(3575, .1)%>% ##Sec nitrous oxide
  scaleClade(2875, .01)%>% ##Sec nitrous oxide
  scaleClade(1930, .01)%>% ##Sec nitrous oxide
  scaleClade(1892, .01)%>% ##Sec nitrous oxide
  scaleClade(3360, .01)%>% ##Sec nitrous oxide
  scaleClade(3417, .01)%>% ##Sec nitrous oxide
  scaleClade(3457, .01)%>% ##Sec nitrous oxide
  scaleClade(3570, .01)%>% ##Sec nitrous oxide
  scaleClade(3130, .01)%>% ##Sec nitrous oxide
  scaleClade(2893, .01)%>% ##Sec nitrous oxide
  scaleClade(2881, .01)%>% ##Sec nitrous oxide
  scaleClade(3311, .01)%>% ##Sec nitrous oxide
  scaleClade(1907, .01)%>% ##Sec nitrous oxide
  scaleClade(1925, .01)%>% ##Sec nitrous oxide
  scaleClade(3549, .01)%>% ##Sec nitrous oxide
  scaleClade(1922, .01)%>% ##Sec nitrous oxide
  scaleClade(3303, .01)%>% ##Sec nitrous oxide
  scaleClade(3300, .01)%>% ##Sec nitrous oxide
  scaleClade(1858, .01)%>% ##Sec nitrous oxide
  scaleClade(3297, .01)%>% ##Sec nitrous oxide
  scaleClade(2754, .01)%>% ##Sec nitrous oxide
  scaleClade(2825, .01)%>% ##Sec nitrous oxide
  scaleClade(2821, .01)%>% ##Sec nitrous oxide
  scaleClade(1954, .005)%>% ##TAT nitrous oxide
  collapse(3581, 'min', fill='orange',color='black')%>% ##good
  collapse(1854, 'min', fill='orange',color='black')%>% ##good
  collapse(3553, 'min', fill='orange',color='black')%>% ##good
  collapse(3434, 'min', fill='orange',color='black')%>% ##good
  collapse(3389, 'min', fill='orange',color='black')%>% ##good
  collapse(3575, 'min', fill='orange',color='black')%>%
  collapse(2821, 'max', fill='orange',color='black')%>%
  collapse(2874, 'min', fill='orange',color='black')%>%
  collapse(1930, 'max', fill='orange',color='black')%>%
  collapse(1892, 'min', fill='orange',color='black')%>%
  collapse(1907, 'min', fill='orange',color='black')%>%
  collapse(1925, 'min', fill='orange',color='black')%>%
  collapse(1922, 'min', fill='orange',color='black')%>%
  collapse(2754, 'min', fill='orange',color='black')%>%
  collapse(2825, 'min', fill='orange',color='black')%>%
  collapse(1858, 'max', fill='pink',color='black')%>%##
  collapse(3570, 'max', fill='pink',color='black')%>%
  collapse(3549, 'max', fill='pink',color='black')%>%
  collapse(3580, 'max', fill='pink',color='black')%>%
  collapse(3495, 'min', fill='green',color='black')%>%
  collapse(1954, 'min', fill='blue',color='black')%>%
  collapse(3475, 'min', fill='yellow',color='black')

scaled_tree
scaled_tree+ylim(200,250)+xlim(0,45)

##collapse(1871, 'max', fill='grey',color='black')%>%
##  collapse(1950, 'max', fill='pink',color='black')%>%

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="DPANN_protein_tree.svg", scaled_tree)




