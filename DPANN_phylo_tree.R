library(readxl)
library(xlsx)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(treeio)
library(cowplot)
library(ape)
library(phylotools)
library(RColorBrewer)

setwd("/Users/irene/Desktop/DPANN/") 
DPANN_MAD<-read.tree("DPANN_derep_GToTree_out/DPANN_derep_GToTree_MAD.tre")
DPANN_MAD$node.label = as.numeric(DPANN_MAD$node.label)
DPANN_MAD$node.label = round(DPANN_MAD$node.label * 100)
DPANN_MAD$node.label
DPANN_MAD$node.label = as.character(DPANN_MAD$node.label)
DPANN_MAD%>% ggtree()

DPANN_tibble <- as_tibble(DPANN_MAD)

DPANN_annote <- read_excel("Excel_files/DPANN_tree_annotation.xlsx", 
                           sheet="annotations")

DPANN_merge<-left_join(DPANN_tibble,DPANN_annote,by="label")
#write.xlsx(DPANN_merge, file = "DPANN_merge_tibble.xlsx", sheetName="DPANN")

DPANN_treedata<-as.treedata(DPANN_merge)

DPANN_treedata %>% ggtree()

#view color palette
display.brewer.pal(8,name='Dark2')
brewer.pal(n=8,"Dark2")

DPANN_tree_taxa<-DPANN_treedata %>% ggtree() + 
  geom_tippoint(aes(fill=Taxa),shape=21,size=2,stroke=NA)+ 
  scale_fill_brewer(palette='Dark2',direction=-1)+
  theme(legend.position = "right")+
  geom_text(aes(label=label), hjust=-.25, size = 2)+
  ggtitle('DPANN tree')

DPANN_tree<-DPANN_treedata %>% ggtree() + 
  geom_tippoint(shape=21, aes(color=source), size=3,stroke=1) + 
  scale_color_manual(values=c("public"='NA',"ODZ"="black",TARA="blue")) +
  ggnewscale::new_scale_fill() +
  geom_tippoint(aes(fill=Taxa),shape=21,size=3, stroke=NA)+ 
  scale_fill_brewer(palette='Dark2',direction=-1)+
  theme(legend.position = "right")+
  geom_text(aes(label=label), hjust=-.25, size = 2)+
  ggtitle('Phylogenetic tree of DPANN archaea')+ylim(0,800)


DPANN_tree_source<-DPANN_treedata %>% ggtree() + 
  geom_tippoint(aes(fill=source),shape=21,size=2,stroke=NA)+ 
  theme(legend.position = "right")+
  geom_text(aes(label=node), hjust=-.25, size = 2)+
  ggtitle('DPANN tree')


scaled_tree<-scaleClade(DPANN_tree, 868, .2)%>%scaleClade(916, .3)%>%
  scaleClade(591, .3)%>%scaleClade(831, .3)%>%
  scaleClade(604, .3)%>%scaleClade(664, .3)%>%
  scaleClade(491, .1)%>%scaleClade(680, .3)%>%
  scaleClade(553, .3)%>%scaleClade(770, .3)%>%
  scaleClade(698, .3)%>%scaleClade(783, .3)%>%
  scaleClade(825, .3)%>%scaleClade(755, .3)%>%
  scaleClade(725, .3)%>%scaleClade(791, .3)%>%
  scaleClade(713, .3)%>%scaleClade(731, 1.3)%>%
  scaleClade(831, .3)%>%scaleClade(605, .3)%>%
  scaleClade(803, 1.5)%>%scaleClade(807, 2)%>%
  scaleClade(767, 3)%>%scaleClade(587, 2)%>%
  scaleClade(586, 2.2)%>%scaleClade(734, 2.2)%>%
  scaleClade(794, 3.5)%>%scaleClade(488, 2.8)%>%
  scaleClade(818, 1.8)%>%
  collapse(868, 'min', fill='darkgrey',color='black') %>% 
  collapse(916, 'min', fill='darkgrey',color='black') %>%
  collapse(798, 'min', fill='darkgrey',color='black') %>%
  collapse(604, 'min', fill="#66A61E",color='black')%>%
  collapse(591, 'min', fill='#66A61E',color='black')%>%
  collapse(664, 'min', fill='#66A61E',color='black')%>%
  collapse(491, 'min', fill='#66A61E',color='black')%>%
  collapse(680, 'min', fill='#66A61E',color='black')%>%
  collapse(553, 'min', fill='#66A61E',color='black')%>%
  collapse(580, 'min', fill='#66A61E',color='black')%>%
  collapse(590, 'min', fill='#66A61E',color='black')%>%
  collapse(589, 'min', fill='#66A61E',color='black')%>%
  collapse(576, 'min', fill='#66A61E',color='black')%>%
  collapse(484, 'min', fill='#66A61E',color='black')%>%
  collapse(490, 'min', fill='#66A61E',color='black')%>%
  collapse(693, 'min', fill='#66A61E',color='black')%>%
  collapse(663, 'min', fill='#66A61E',color='black')%>%
  collapse(588, 'min', fill='#66A61E',color='black')%>%
  collapse(770, 'min', fill='#1B9E77',color='black')%>%
  collapse(783, 'min', fill='#1B9E77',color='black')%>%
  collapse(755, 'min', fill='#1B9E77',color='black')%>%
  collapse(725, 'min', fill='#1B9E77',color='black')%>%
  collapse(762, 'min', fill='#1B9E77',color='black')%>%
  collapse(732, 'min', fill='#1B9E77',color='black')%>%
  collapse(698, 'min', fill='#E6AB02',color='black')%>%
  collapse(831, 'min', fill='#7570B3',color='black')%>%
  collapse(825, 'min', fill='#A6761D',color='black')%>%
  collapse(815, 'min', fill='#A6761D',color='black')%>%
  collapse(821, 'min', fill='#A6761D',color='black')%>%
  collapse(820, 'min', fill='#A6761D',color='black')%>%
  collapse(791, 'min', fill='#E6AB02',color='black')%>%
  collapse(721, 'min', fill='#E6AB02',color='black')%>%
  collapse(789, 'min', fill='#E6AB02',color='black')%>%
  collapse(719, 'min', fill='#E6AB02',color='black')%>%
  collapse(804, 'min', fill='#D95F02',color='black')%>%
  collapse(808, 'min', fill='#D95F02',color='black')%>%
  collapse(713, 'min', fill='#E6AB02',color='black')

scaled_tree<-scaled_tree+ylim(50,330)+xlim(0.25,2)
scaled_tree

setwd("/Users/irene/Desktop/DPANN/figures")
ggsave(filename="DPANN_phylogenetic_tree_labels.svg", scaled_tree, height=11.4, width =  8.97)
