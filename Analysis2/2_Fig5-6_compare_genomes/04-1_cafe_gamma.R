library(tidyverse)

dat<-read.table('Orthogroups.GeneCount.tsv', header = T, sep = '\t')
interpro<-read.table('culicidae.tsv', header = F, sep = '\t', quote = '', comment.char = '')
names(interpro)<-c('accession', 'md5', 'seq.length', 'analysis', 'sig.accession', 'sig.description',
                   'start', 'end', 'score', 'status', 'date', 'interpro.anno.accession',
                   'interpro.anno.description', 'GOannotation' ,'pathways')

interpro2<-interpro %>%
  select(accession, seq.length, sig.accession:score, GOannotation) %>%
  mutate(og = gsub('(OG[0-9]+)(_.+)', '\\1', accession),
         taxon = gsub('(OG[0-9]+_)([A-z]+_[a-z]+_?L?5?)(_g[0-9]+\\..+)', '\\2', accession, perl = T),
         gene = gsub('(OG[0-9]+_)([A-z]+_[a-z]+_?L?5?)(_g[0-9]+\\..+)', '\\3', accession, perl = T),
         gene = gsub('_', '', gene)) %>%
  select(-accession) %>%
  group_by(og) %>%
  slice_min(order_by = score, n = 1, with_ties = F) %>%
  relocate(c(og, taxon, gene), .before = seq.length)

dat2<-dat %>% 
  left_join(interpro2 %>% select(og, sig.accession, sig.description, GOannotation), by = c('Orthogroup' = 'og'))


og_variances<-dat2 %>%
  select(-Total) %>%
  rowwise %>%
  mutate(var = var(c_across(Aedes_aegypti_L5:Sabethes_cyaneus))) %>%
  ungroup %>%
  select(Orthogroup, var)

og_diffs<-dat2 %>%
  select(-Total) %>%
  rowwise %>%
  mutate(diffs = max(c_across(where(is.numeric))) - min(c_across(where(is.numeric)))) %>%
  ungroup %>%
  select(Orthogroup, diffs)


og_var_keep<-og_variances %>% 
  filter(var < quantile(var, 0.95)) %>%
  pull(Orthogroup)

og_diff_keep<-og_diffs %>%
  filter(diffs < quantile(diffs, 0.95)) %>%
  pull(Orthogroup)


#output combined OrthogroupCounts and in Interpro
write.table(dat2, file = 'OrthogroupCounts_Interpro.tsv', sep = '\t', row.names = F, col.names = T,
            quote = F)


#output filtered list for Cafe5 analysis
# filter out 95th percentile of highest variance

dat2_var_filt<-dat2 %>%
  filter(Orthogroup %in% og_var_keep) %>%
  select(Orthogroup:Sabethes_cyaneus, sig.description) %>%
  relocate(sig.description, .before = Orthogroup) %>%
  rename(Description = sig.description)

write.table(dat2_var_filt, file = 'OrthogroupCounts_Interpro_var_filt.tsv', sep = '\t', row.names = F, quote = F)
