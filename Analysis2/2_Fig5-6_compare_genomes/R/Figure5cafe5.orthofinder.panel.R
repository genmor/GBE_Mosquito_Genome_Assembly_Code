library(tidyverse)
library(patchwork)
library(cowplot)
library(phytools)
library(viridis)
library(ggtree)
library(ggridges)
library(ggnewscale)
library(readr)


tree.list<-read.nexus('Gamma_asr.tre')
tree<-read.tree('SpeciesTreeAlignment.fa.timetree.nwk')
# ex.con.dat<-read.table('CONTRACTION_EXPANSION.txt', header = T)
gamma.clade.result<-read.table('Gamma_clade_results.txt', sep = '\t')
names(gamma.clade.result)<-c('node', 'increase', 'decrease')
gamma.branch.probs<-read.table('Gamma_branch_probabilities.tab', header = T)
names(gamma.branch.probs)[grep('X\\.[0-9]+\\.', names(gamma.branch.probs))]<-gsub('(X\\.)([0-9]+)(\\.)','node_\\2', names(gamma.branch.probs)[grep('X\\.[0-9]+\\.', names(gamma.branch.probs))])
names(gamma.branch.probs)[grep('\\.[0-9]+\\.', names(gamma.branch.probs))]<-gsub('(\\.)([0-9]+)(\\.)','', names(gamma.branch.probs)[grep('\\.[0-9]+\\.', names(gamma.branch.probs))])
gamma.fam.lik<-read.table('Gamma_family_likelihoods.txt')
names(gamma.fam.lik)<-c('Orthogroup', 'gamma_cat', 'cat_lik', 'fam_lik', 'post_prob', 'sig')
ogs<-read.table('OrthogroupCounts_Interpro.tsv', sep = '\t', comment.char = '', quote = '', header= T)
og.summary<-read.delim('clipboard', header = T, sep = '\t') ###copy the first table from statistics_perspecies.tsv from orthofinder

og.summary2<-og.summary %>%
  pivot_longer(cols = -X, names_to = 'sp', values_to = 'count') %>%
  filter(X == 'Number of genes' | X == 'Number of unassigned genes') %>%
  mutate(X = case_when(X == 'Number of genes' ~ 'num.genes',
                       .default = 'unassigned')) %>%
  pivot_wider(names_from = 'X', values_from = 'count')


#count single copy HOGs
single.copy<-ogs %>%
  select(-c(Total, sig.accession, sig.description, GOannotation)) %>%
  filter(if_all(where(is.numeric), ~  .x == 1)) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'sp', values_to = 'count') %>%
  group_by(sp) %>%
  summarise(single.copy = sum(count))

og.total<-ogs %>%
  select(-c(Total, sig.accession, sig.description, GOannotation)) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'sp', values_to = 'count') %>%
  group_by(sp) %>%
  summarise(og.total = sum(count))

uniq.para<-ogs %>%
  filter(if_any(contains('_'), ~ .x == Total)) %>% 
  select(-c(sig.accession, sig.description, GOannotation)) %>%
  pivot_longer(cols = -Orthogroup, names_to = 'sp', values_to = 'count') %>%
  group_by(sp) %>%
  summarise(uniq.para = sum(count))

og.summary.all<-og.summary2 %>%
  left_join(og.total) %>%
  left_join(single.copy) %>%
  # left_join(multi.copy) %>%
  left_join(uniq.para) %>%
  mutate(multi.copy = og.total - (single.copy + unassigned + uniq.para)) %>% 
  arrange(match(sp, tree$tip.label)) %>%
  mutate(tip.order = row_number())



lab.key<-rbind(data.frame(cafe.lab = as.numeric(gsub('(.+)(<)([0-9]+)(.+)', '\\3', tree.list[[2]]$tip.label)),
                          node.lab = gsub('(.+)(<)([0-9]+)(.+)', '\\1', tree.list[[2]]$tip.label)),
               data.frame(cafe.lab = as.numeric(gsub('(<)([0-9]+)(.+)', '\\2', tree.list[[2]]$node.label)),
                          node.lab = (Ntip(tree.list[[2]]) + 1):(Ntip(tree.list[[2]]) -1 + Ntip(tree.list[[2]])))
) %>%
  arrange(as.numeric(cafe.lab)) %>%
  mutate(col_labs = case_when(grepl('_', node.lab) ~ node.lab,
                              .default = paste0('node_', node.lab))) %>%
  left_join(data.frame(tree.node = 1:(Ntip(tree) + Nnode(tree)),
                       tree.node.lab = c(tree$tip.label, rep('', Ntip(tree)- 1))
  ) %>%
    mutate(tree.node.lab = case_when(tree.node.lab == '' ~ paste0('node_', tree.node),
                                     .default = tree.node.lab)),
  by = c('col_labs' = 'tree.node.lab')
  )
names(gamma.branch.probs)<-c('Orthogroup', lab.key$col_labs)


all.res<-left_join(gamma.clade.result %>%
            separate(node, sep = '<', into = c('taxon', 'cafe_node')) %>%
            mutate(cafe_node = as.numeric(gsub('>', '', cafe_node))) %>%
            left_join(lab.key, by = c('cafe_node' = 'cafe.lab')) %>%
            rename(label = taxon) %>%
            rowwise %>%
            mutate(ratio = case_when(increase > decrease ~ increase/decrease,
                                     decrease > increase ~ -1 * decrease/increase),
                   ratio = case_when(is.infinite(ratio) ~ NA, .default = ratio)) %>%
            ungroup,
          gamma.branch.probs %>%
            mutate(node_13 = NA) %>%
            pivot_longer(cols = -Orthogroup, names_to = 'nodes', values_to = 'pval') %>%
            mutate(fdr = p.adjust(p = pval, method = 'fdr'),
                   rapid = fdr < 0.01) %>%
            group_by(nodes) %>%
            summarise(rapid = sum(rapid)),
          by = c('col_labs' = 'nodes')) %>%
  select(label, tree.node, increase, decrease, ratio, rapid) %>%
  mutate(label2 = gsub('_', ' ', label)) %>%
  rename(node = tree.node)

all.res2<-all.res %>%
  mutate(ratio = case_when(ratio > 20 ~ NA, .default = ratio))






(cafe5.plot<-ggtree(tree, right = T) %<+% all.res2 +
  # geom_tiplab(aes(label = label2), offset = 55, size = 3) +
  geom_tiplab(aes(label = rapid), offset = 0.15, size = 3, color = '#21908CFF') +
  geom_tiplab(aes(label = paste0('+',increase)), size  = 3, nudge_x = 0.05, nudge_y = 0.15) +
  geom_tiplab(aes(label = paste0('-',decrease)), size  = 3, nudge_x = 0.05, nudge_y = -0.15) +
  geom_nodelab(aes(label = paste0('+',increase)), size  = 3, nudge_x = 0.05, nudge_y = 0.15) +
  geom_nodelab(aes(label = paste0('-',decrease)), size  = 3, nudge_x = 0.05, nudge_y = -0.15) +
  geom_nodelab(aes(label = rapid), size = 3, nudge_x = 0.125, color = '#21908CFF') +
  geom_tippoint(aes(fill = ratio, size = abs(ratio)), shape = 21) +
  geom_nodepoint(aes(fill = ratio, size = abs(ratio)), shape = 21) +
  scale_fill_viridis(option = 'rocket', 'Ratio', breaks = scales::pretty_breaks()) +
  xlim(0, 1.25) +
  geom_treescale(x  = 0.1, y = 7.5) +
  guides(size = 'none',
         fill = guide_colorbar(direction = 'horizontal', position = 'inside')) +
  theme(aspect.ratio = 1.25,
        legend.position.inside = c(0.15, 0.3))
)

(og.summary.plot<-og.summary.all %>% 
  rename(`Single copy` = single.copy, `Unique paralogs` = uniq.para, `Unassigned genes`= unassigned, `Multi-copy` = multi.copy) %>%
  pivot_longer(cols = c(-sp, -og.total, -num.genes, -tip.order), names_to = 'cat', values_to = 'vals') %>%
  mutate(cat = factor(cat, ordered = T, levels = rev(c('Single copy', 'Multi-copy', 'Unique paralogs', 'Unassigned genes'))),
         sp = factor(sp, levels = ladderize(tree, right = F)$tip.label)) %>% 
  ggplot(aes(x = tip.order, y = vals/1000, fill = cat)) + 
  geom_bar(position = 'stack', stat = 'identity', color = 'black') +
  geom_text(data = og.summary.all %>% mutate(vals = 0, cat = 'Single copy', sp2 = gsub('_', ' ', sp)),
            aes(x = tip.order, y = vals, label = sp2), color = 'white', hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = scales::pretty_breaks(n = 13), 'Number of orthologs (x 1000)') +
  scale_x_continuous(labels = og.summary.all$sp, breaks = 1:12, expand = expansion(mult = c(0, 0.0))) +
  scale_fill_viridis(option = 'viridis', discrete = T, direction = -1, 'Orthogroup category') +
  coord_flip() +
  guides(fill = guide_legend(position = 'inside')) +
  theme_minimal_vgrid() +
  theme(legend.position.inside = c(0.6, 0.9),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1
        )
)

cafe5.plot<-cafe5.plot + scale_y_continuous(expand = c(0, 0.1, 0, 0.1)) + theme(aspect.ratio = 1)

og.summary.plot.legend<-orthogroup_summary %>% 
  rename(`Single copy` = single.copy, `Unique paralogs` = uniq.para, `Unassigned genes`= unassigned, `Multi-copy` = multi.copy) %>%
  pivot_longer(cols = c(-sp, -total, -n.seq, -tip.order), names_to = 'cat', values_to = 'vals') %>%
  mutate(cat = factor(cat, ordered = T, levels = rev(c('Single copy', 'Multi-copy', 'Unique paralogs', 'Unassigned genes'))),
         sp = factor(sp, levels = ladderize(tree, right = F)$tip.label)) %>% 
  ggplot(aes(x = tip.order, y = vals/1000, fill = cat)) + 
  geom_bar(position = 'stack', stat = 'identity', color = 'black') +
  geom_text(data = orthogroup_summary %>% mutate(vals = 0, cat = 'Single copy', sp2 = gsub('_', ' ', sp)),
            aes(x = tip.order, y = vals, label = sp2), color = 'white', hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = scales::pretty_breaks(n = 12), 'Number of orthologs (x 1000)') +
  scale_x_continuous(labels = orthogroup_summary$sp, breaks = 1:12, expand = expansion(mult = c(0, 0.0))) +
  scale_fill_viridis(option = 'viridis', discrete = T, direction = -1) +
  coord_flip() +
  theme_minimal_vgrid() +
  theme(legend.position = 'left',
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1)


(cafe5.plot.legend<-cafe5.plot + theme(legend.position = 'left'))


(panel1<-cafe5.plot + plot_spacer() + og.summary.plot + plot_layout(ncol = 3, widths = c(1, -0.4, 1)))
ggsave2(filename = 'cafe5.orthofinder.panel.pdf', plot = panel1, 
        width = 175, height = 175, units = 'mm')

# ggsave2(filename = 'cafe5.legend.pdf', plot = cafe5.plot.legend, 
#         width = 175, height = 175, units = 'mm')
# ggsave2(filename = 'og.summary.plot.legend.pdf', plot = og.summary.plot.legend, 
#         width = 175, height = 175, units = 'mm')
save.image('cafe5.orthofinder.panel.RData')
