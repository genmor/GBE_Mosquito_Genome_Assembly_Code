library(tidyverse)
library(viridis)
library(cowplot)
library(phytools)
library(patchwork)
library(ggridges)


og.counts<-read.table('OrthogroupCounts_Interpro.tsv', header = T, quote = '',
                      comment.char = '', sep = '\t')
tree.list<-read.nexus('Gamma_asr.tre')
tree<-read.tree('SpeciesTreeAlignment.fa.timetree.nwk')
sig.og.dat<-read.table('Gamma_family_results.txt', sep = '\t', header =  F)
rate.cats<-read.table('Gamma_family_likelihoods.txt', sep = '\t', header= F)
branch.prob<-read.table('Gamma_branch_probabilities.tab', sep = '\t', header = T)
change.dat<-read.table('Gamma_change.tab', sep = '\t', header = T)
gamma.counts<-read.table('Gamma_count.tab', sep = '\t', header = T)


#create key for cafe5 node labels
cafe5.key<-rbind(data.frame(cafe.lab = as.numeric(gsub('(.+)(<)([0-9]+)(.+)', '\\3', tree.list[[1]]$tip.label)),
                 node.lab = gsub('(.+)(<)([0-9]+)(.+)', '\\1', tree.list[[1]]$tip.label)),
      data.frame(cafe.lab = gsub('(<)([0-9]+)(.+)', '\\2', tree.list[[1]]$node.label),
                 node.lab = (Ntip(tree.list[[1]])+1):((Ntip(tree.list[[1]])-1)+Ntip(tree.list[[1]]))
      )) %>%
  arrange(as.numeric(cafe.lab)) %>%
  mutate(col_labs = case_when(grepl('_', node.lab) ~ node.lab,
                              .default = paste0('node_',node.lab)))

#change or add column names
names(sig.og.dat)<-c('Orthogroup', 'pval', 'sig_0.01')
names(rate.cats)<-c('Orthogroup', 'gamma_cat', 'cat_lik', 'fam_lik', 'posterior', 'sig')
names(change.dat)<-c('Orthogroup', cafe5.key$col_labs)
names(gamma.counts)<-c('Orthogroup', cafe5.key$col_labs)
names(branch.prob)<-c('Orthogroup', cafe5.key$col_labs)

sp.change.long<-change.dat %>%
  select(Orthogroup, !contains('node')) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'taxon', values_to = 'change')


gamma.counts<-gamma.counts %>%
  select(Orthogroup, tree$tip.label)

names(branch.prob)<-c('Orthogroup', cafe5.key$col_labs)

sp.sig.change<-sapply(tree.list, function(x) {
  sp<-gsub('(.+)(<.+)', '\\1', x$tip.label[grep('\\*', x$tip.label)])
  x<-tree$tip.label %in% sp
  names(x)<-tree$tip.label
  return(x)
})

sp.sig.change.long<-data.frame(sp.sig.change) %>%
  rownames_to_column('taxon') %>%
  pivot_longer(cols = -taxon, names_to = 'Orthogroup', values_to = 'sig_change')

branch.prob.long<-branch.prob %>%
  dplyr::select(!contains('node')) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'taxon', values_to = 'pval') %>%
  mutate(fdr = p.adjust(p = pval, method = 'fdr'))

gamma.counts.long<-gamma.counts %>%
  select(Orthogroup, !contains('node')) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'taxon', values_to = 'ortho_count') %>%
  left_join(og.counts %>% select(Orthogroup, sig.description, GOannotation))


dat.overall<-gamma.counts.long%>%
  left_join(sp.change.long) %>%
  left_join(branch.prob.long) %>%
  left_join(sp.sig.change.long) %>%
  mutate(rapid = case_when(fdr <= 0.01 ~ T,
                           is.na(fdr) ~ F,
                           .default = F),
         taxon = factor(taxon, levels = tree$tip.label)
  )


cafe5.ridges<-ggplot(dat.overall, aes(y = taxon, x = change, fill = rapid)) +
  geom_density_ridges(alpha = 0.5, scale = 0.95, rel_min_height = 0.01) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), 'Copy number change',
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_viridis(discrete = T, begin = 0, end = 1, 'Rate category',
                     labels = c('Not rapid', 'Rapid')) +
  theme_minimal_grid() +
  theme(axis.title.y = element_blank(),
        aspect.ratio = 1,
        legend.position = 'bottom')

ggsave2(filename = 'cafe5_ridges.pdf', plot = cafe5.ridges, width = 8, height = 6, units = 'in')


ex.fams<-dat.overall %>%
  filter(sig.description != '', sig.description != '-', 
         sig.description != 'UNCHARACTERIZED', sig.description != 'FAMILY NOT NAMED', 
         GOannotation != '', GOannotation != '-', rapid == T, sig_change == T) %>%
  group_by(sig.description, GOannotation) %>%
  summarise(count = n()) %>%
  ungroup %>%
  slice_max(count, n = 5) %>%
  pull(sig.description)

tmp<-dat.overall %>% 
  filter(fdr < 0.01, sig.description %in% ex.fams) %>% 
  group_by(Orthogroup, sig.description) %>% 
  summarise(count = n()) %>%
  ungroup() %>% 
  group_by(sig.description) %>% mutate(shape = row_number()) %>%
  ungroup() %>%
  select(Orthogroup, shape)


(cafe5.jitter<-ggplot() +
  geom_jitter(data = dat.overall %>% filter(fdr > 0.01), aes(x = change, y = taxon),
              position = position_jitter(width = 0, height = 0.35), size = 2, shape = 21, fill = 'white') +
  geom_jitter(data = dat.overall %>% filter(fdr < 0.01, !(sig.description %in% ex.fams)), aes(x = change, y = taxon),
              position = position_jitter(width = 0, height = 0.35), size = 2, shape = 21, fill = 'grey78') +
  geom_jitter(data = left_join(dat.overall %>% filter(fdr < 0.01, sig.description %in% ex.fams), tmp), aes(x = change, y = taxon, fill = sig.description, shape = as.factor(shape)),
              position = position_jitter(width = 0, height = 0.35), size = 3) +
  scale_fill_viridis(discrete = T, option = 'turbo') +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), 'Copy number change') +
  guides(fill = guide_legend(nrow = 2, title.position = 'top', override.aes = list(shape = 21))) +
  theme_minimal_grid() +
  theme(axis.title.y = element_blank(),
        aspect.ratio = 0.75,
        legend.position = 'bottom',
        axis.text = element_text(color = 'black', size = 9),
        legend.text = element_text(size = 9)))

ggsave2(filename = 'cafe5_jitter.pdf', plot = cafe5.jitter, height = 6,
        width = 8, units = 'in')



save.image('cafe5_ridge_jitter.plots.RData')

