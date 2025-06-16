library(tidyverse)
library(ape)
library(simplifyEnrichment)
library(viridis)
library(cowplot)
library(magick)
library(patchwork)


gamma.changes<-read.table('Gamma_change.tab', header = T, sep = '\t')
ogs<-read.table('OrthogroupCounts_Interpro.tsv', header = T, sep = '\t', quote = '')

gamma.changes<-gamma.changes %>%
  select(FamilyID:Phlebotomus_papatasi.12.)
names(gamma.changes)<-gsub('\\.[0-9]+\\.', '', names(gamma.changes))

####expanded orthogroups for taxa of interest
L5.expanded<-gamma.changes %>%
  select(FamilyID, Aedes_aegypti_L5) %>%
  filter(Aedes_aegypti_L5 > 0) %>%
  left_join(ogs %>% select(Orthogroup, panther_id, panther_fam, go), by = c('FamilyID' = 'Orthogroup'))

aaf.expanded<-gamma.changes %>%
  select(FamilyID, Aedes_aegypti_formosus) %>%
  filter(Aedes_aegypti_formosus > 0) %>%
  left_join(ogs %>% select(Orthogroup, panther_id, panther_fam, go), by = c('FamilyID' = 'Orthogroup'))

am.expanded<-gamma.changes %>%
  select(FamilyID, Aedes_mascarensis) %>%
  filter(Aedes_mascarensis > 0) %>%
  left_join(ogs %>% select(Orthogroup, panther_id, panther_fam, go), by = c('FamilyID' = 'Orthogroup'))


ag.expanded<-gamma.changes %>%
  select(FamilyID, Anopheles_gambiae) %>%
  filter(Anopheles_gambiae > 0) %>%
  left_join(ogs %>% select(Orthogroup, panther_id, panther_fam, go), by = c('FamilyID' = 'Orthogroup')) %>% nrow


L5.go<-L5.expanded$go
L5.go<-L5.go[which(L5.go != '-')]
L5.go<-unlist(lapply(strsplit(L5.go, split = '\\|'), `[`))
L5.go<-L5.go[-grep('InterPro', L5.go)]
L5.go<-gsub('\\([A-z]+\\)', '', L5.go)

aaf.go<-aaf.expanded$go
aaf.go<-aaf.go[which(aaf.go != '-')]
aaf.go<-unlist(lapply(strsplit(aaf.go, split = '\\|'), `[`))
aaf.go<-aaf.go[-grep('InterPro', aaf.go)]
aaf.go<-gsub('\\([A-z]+\\)', '', aaf.go)

am.go<-am.expanded$go
am.go<-am.go[which(am.go != '-')]
am.go<-unlist(lapply(strsplit(am.go, split = '\\|'), `[`))
am.go<-am.go[-grep('InterPro', am.go)]
am.go<-gsub('\\([A-z]+\\)', '', am.go)


aaf.go %in% go.simp$id %>% sum
am.go %in% go.simp$id %>% sum
L5.go %in% go.simp$id %>% sum

sum(!(go.simp$id %in% c(aaf.go[aaf.go %in% go.simp$id], L5.go[L5.go %in% go.simp$id]))) ##unique to Am 
sum(!(go.simp$id %in% c(aaf.go[aaf.go %in% go.simp$id], am.go[am.go %in% go.simp$id]))) ##unique to L5 
sum(!(go.simp$id %in% c(am.go[am.go %in% go.simp$id], L5.go[L5.go %in% go.simp$id]))) ##unique to aaf



#pdf('kmeans_simplify_go_heat.pdf', width = 6.89, height = 6.89) #doesn't work well
#just export from GUI
go.heat<-simplifyGOFromMultipleLists(lt = list(L5 = L5.go, Aaf= aaf.go, Am = am.go), ont = 'BP', 
                            db = 'org.Ag.eg.db', method = 'kmeans', 
                            padj_cutoff = 0.01, measure = 'Sim_XGraSM_2013',
                            show_barplot = F,
                            max_words = 7,
                            min_term = 7,
                            heatmap_param = list(use_raster = F)
                           )
#dev.off()

go.mat<-GO_similarity(go_id = unique(c(L5.go, aaf.go, am.go)), ont = 'BP', db = 'org.Ag.eg.db',
              measure = 'Sim_XGraSM_2013')

go.simp<-simplifyGO(mat = go.mat, method = 'kmeans', plot = F)
tile.dat<-data.frame(go.mat) %>%
  rownames_to_column('y') %>%
  pivot_longer(cols = -y, names_to = 'x', values_to = 'val') %>% 
  left_join(go.simp, by = c('y' = 'id')) %>% 
  rename(y.cluster = cluster) %>%
  mutate(x = gsub('\\.', ':', x)) %>%
  left_join(go.simp, by = c('x' = 'id')) %>%
  rename(x.cluster = cluster) %>%
  mutate(y = factor(y, levels = rev(unique(y[order(y.cluster)])), ordered = T),
         x = factor(x, levels = unique(x[order(x.cluster)]), ordered = T))

y.bounds<-tile.dat %>% 
  group_by(y.cluster) %>%
  mutate(row = row_number()) %>% 
  slice_max(row, n = 1)

x.bounds<-tile.dat %>% 
  group_by(x.cluster) %>%
  mutate(row = row_number()) %>% 
  slice_max(row, n = 1)

word.grob<-anno_word_cloud_from_GO(align_to = tile.dat$y.cluster,
                        go_id = tile.dat$y,
                        stat = 'count')

#creates heatmap in cluster order with cluster bounds  
 go.heat<-ggplot() +
  geom_tile(data = tile.dat, aes(x = x, y = y, fill = val)) +
  geom_hline(data = y.bounds, aes(yintercept = y)) +
  geom_vline(data = x.bounds, aes(xintercept = x)) +
  scale_fill_gradient(low = 'white', high = '#08306b') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

setOldClass('AnnotationFunction', )
go.heat + annotation_custom(word.grob)
 
 summarizeGO(go_id = )


  #theme_void()
#methods
# dynamicTreeCut 69
# louvain 7
# kmeans 18
# apcluster 43
# walktrap 16
# MCL 9
# mclust 49
# pam 62
# binary_cut 26
# hdbscan 36