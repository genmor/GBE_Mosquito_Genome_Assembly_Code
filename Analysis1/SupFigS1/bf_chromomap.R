library(chromoMap)
library(viridis)
library(cowplot)
library(tidyverse)
library(ggrepel)
#copy output from megablast
dat<-read.table('2024_05_29_Bf05_dc_canu_ph_insp_yahs_jb_no_debris_tgs.blast.out', header = F, sep = '\t')
names(dat)<-c('bac', 'bf_seq','percent_match','overlap_length', 'num_mismatch', 'gapopen', 'tran_start', 'tran_end', 'bf_start', 'bf_end','e_val', 'bitscore')
tr_top<-dat %>%
  group_by(bac) %>%
  top_n(n = 1, bitscore) 

L5_map<-read.table('bac_pos.txt', header = T, sep = '\t') #bac_pos.txt

tr_top2<-tr_top %>%
  mutate(bac = gsub('(*)(\\.[0-9])', '\\1', bac)) %>%
  left_join(L5_map) %>%
  ungroup() %>%
  select(region, bf_seq, bf_start, bf_end, arm)

#create a tab delimited file consisting of sequence name, start, and end from genome fasta. copy and read in
chr.dat<-read.table('bf_scaf.bed', header = F, sep = '\t')
# shift start and end by +1
chr.dat$V2<-chr.dat$V2 + 1
chr.dat$V3<-chr.dat$V3 + 1

chr.dat2<-chr.dat %>%
  filter(V1 %in% tr_top2$bf_seq)


chromoMap(list(chr.dat2), list(tr_top2),
          data_based_color_map = T,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          title = 'Bf',
          data_colors = list(viridis(n = 6, begin = 0.25, end = 1, option = 'turbo'))
          )

dat2<-dat %>%
  mutate(bac = gsub('(*)(\\.[0-9])','\\1', bac))


tr_top3<-tr_top2 %>%
  mutate(bf_seq = case_when(bf_seq == 'scaffold_1' ~ 'chr_2',
                              bf_seq == 'scaffold_2' ~ 'chr_3',
                              bf_seq == 'scaffold_3' ~ 'chr_1',
                              .default = bf_seq))

chr.dat3<-chr.dat2 %>%
  mutate(V1 = case_when(V1 == 'scaffold_1' ~ 'chr_2',
                            V1 == 'scaffold_2' ~ 'chr_3',
                            V1 == 'scaffold_3' ~ 'chr_1',
                            .default = V1))

chromoMap(list(chr.dat3), list(tr_top3), 
          data_based_color_map = T,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Bf',
          chr_color = list('grey77'),
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)

bad_map<-tr_top3 %>%
  mutate(chr_match = case_when(gsub('chr_', '', bf_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F),
         arm_map = case_when(bf_seq == 'chr_1' ~ ifelse(bf_end/1000000 > 152, 'q', 'p'),
                             bf_seq == 'chr_2' ~ ifelse(bf_end/1000000 > 230, 'q', 'p'),
                             bf_seq == 'chr_3' ~ ifelse(bf_end/1000000 > 198, 'q', 'p'),
                             grepl('scaffold_', bf_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F)) %>%
  filter(arm_match == F | chr_match == F)




good_map<-tr_top3 %>%
  mutate(chr_match = case_when(gsub('chr_', '', bf_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F),
         arm_map = case_when(bf_seq == 'chr_1' ~ ifelse(bf_end/1000000 > 152, 'q', 'p'),
                             bf_seq == 'chr_2' ~ ifelse(bf_end/1000000 > 230, 'q', 'p'),
                             bf_seq == 'chr_3' ~ ifelse(bf_end/1000000 > 198, 'q', 'p'),
                             grepl('scaffold_', bf_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F)) %>%
  filter(arm_match == T & chr_match == T)


chromoMap(list(chr.dat3 %>% filter(!grepl('scaffold_', V1)),
               chr.dat3 %>% filter(!grepl('scaffold_', V1))),
          list(bad_map %>% select(region:arm), good_map %>% select(region:arm)), 
          data_based_color_map = T,
          ploidy = 2,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Aaf',
          chr_color = list('grey77'),
          export.options = T,
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)


good_bad<-tr_top3 %>%
  mutate(chr_match = case_when(gsub('chr_', '', bf_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F),
         arm_map = case_when(bf_seq == 'chr_1' ~ ifelse(bf_end/1000000 > 152, 'q', 'p'),
                             bf_seq == 'chr_2' ~ ifelse(bf_end/1000000 > 230, 'q', 'p'),
                             bf_seq == 'chr_3' ~ ifelse(bf_end/1000000 > 198, 'q', 'p'),
                             grepl('scaffold_', bf_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', bf_seq) ~ NA,
                               .default = F))
  

dat3<-dat2 %>%
  mutate(bf_seq = case_when(bf_seq == 'scaffold_1' ~ 'chr_2',
                            bf_seq == 'scaffold_2' ~ 'chr_3',
                            bf_seq == 'scaffold_3' ~ 'chr_1',
                            .default = bf_seq)) %>%
  left_join(L5_map) %>%
  select(region, bf_start, bf_end, percent_match, overlap_length, bitscore)

(marker_scatter<-left_join(good_bad, dat3) %>%
    filter(!grepl('scaffold', bf_seq)) %>%
    mutate(good_bad = case_when(chr_match == T & arm_match == T ~ 'good',
                                .default = 'bad')) %>%
    ggplot(aes(x = bitscore, y = percent_match, fill = overlap_length, shape = good_bad)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = region), size = 2) +
    scale_shape_manual(values = c(24, 21), 'BAC mapping position', labels = c('Incorrect', 'Correct')) +
    scale_fill_viridis(option = 'cividis', 'Overlap length (bp)') +
    labs(x = 'Bitscore', y = 'Percent match') +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text = element_text(color = 'black'))
)
ggsave2(plot = marker_scatter, filename = 'bf_marker_scatter.pdf',
        width = 6.5, height = 6.5, units = 'in')

#create bed file with tags for plotsr
bf05_bacmap<-good_map %>%
  select(bf_seq, bf_start, bf_end, region) %>%
  mutate(bf_seq = case_when(bf_seq == 'chr_1' ~ 'scaffold_3',
                            bf_seq == 'chr_2' ~ 'scaffold_1',
                            bf_seq == 'chr_3' ~ 'scaffold_2'),
         start = case_when(bf_start < bf_end ~ bf_start,
                           bf_start > bf_end ~ bf_end),
         end = case_when(bf_end > bf_start ~ bf_end,
                         bf_end < bf_start ~ bf_start),
         genome_id = 'Aaf',
         tags = paste('mt:_','mc:black',paste('tt', region, sep = ':'), 'tp:0.2', sep = ';')
         ) %>%
  rename(chr = bf_seq) %>%
  select(chr, start, end, genome_id, tags)

write.table(bf05_bacmap, file = 'bf05_bacmap.txt', quote = F, sep = '\t', row.names = F)
save.image('bf_chromomap.RData')
