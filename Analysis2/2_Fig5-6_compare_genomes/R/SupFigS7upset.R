library(tidyverse)
library(clipr)
library(patchwork)
library(cowplot)
library(viridis)
library(ComplexUpset)
library(tidyverse)
library(phytools)

dat<-read.table('OrthogroupCounts_Interpro.tsv', sep = '\t', header = T, quote = '', comment.char = '', fill = T)
tree<-read.tree('SpeciesTreeAlignment.fa.timetree.nwk')

taxa<-colnames(dat)[2:13]
taxa2<-tree$tip.label
dat_upset<-dat %>% select(Aedes_aegypti_L5:Sabethes_cyaneus, Orthogroup)
dat_upset2<-dat_upset
dat_upset[taxa] <- dat_upset[taxa] > 0
dat_upset2[taxa2] <- dat_upset2[taxa2] > 0


upset.plot<-upset(dat_upset2[taxa2], taxa2, name = 'Taxon',
              intersections = list(
                'Culex_quinquefasciatus',
                'Culex_pipiens_pallens',
                'Sabethes_cyaneus',
                'Armigeres_subalbatus',
                'Aedes_mascarensis',
                'Aedes_aegypti_L5',
                'Aedes_aegypti_formosus',
                'Anopheles_ziemanni',
                'Anopheles_gambiae',
                'Anopheles_darlingi',
                'Anopheles_cruzii',
                'Phlebotomus_papatasi',
                taxa2[grepl('Aedes', taxa2)],
                taxa2[grepl('Anopheles', taxa2)],
                taxa2[grepl('Culex', taxa2)],
                c(taxa2[grepl('Aedes', taxa2)], taxa2[grepl('Culex', taxa2)], 'Armigeres_subalbatus', 'Sabethes_cyaneus'),
                c(taxa2[grepl('Aedes', taxa2)], 'Armigeres_subalbatus'),
                c('Aedes_aegypti_L5', 'Aedes_aegypti_formosus'),
                c('Aedes_aegypti_formosus', 'Aedes_mascarensis'),
                c('Aedes_aegypti_L5', 'Aedes_mascarensis'),
                taxa2[-12],
                taxa2
              ), sort_intersections = 'ascending', sort_sets = F,
              set_sizes = (upset_set_size() +
                             geom_text(aes(label = ..count..), hjust = 1.1, stat = 'count') +
                             theme(axis.text.x = element_text(angle = 45)))
              )
ggsave2(plot = upset.plot, filename = 'upset_plot3.pdf', width = 8, height = 6, units = 'in')


save.image('upset.RData')
