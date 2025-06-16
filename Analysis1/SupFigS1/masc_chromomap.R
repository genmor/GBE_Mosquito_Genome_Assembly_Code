library(tidyverse)
library(cowplot)
library(chromoMap)
library(ggrepel)
library(viridis)
chr.dat<-read.table('masc_scaf.bed', header = F)
names(chr.dat)<-c('masc.seq', 'start', 'end')
chr.dat<-chr.dat %>%
  mutate(start = start + 1,
         end = end + 1)

L5_pos<-read.table('bac_pos.txt', header = T)

dat<-read.table('2024_06_04_MascCH02-dc_hifiasm_ph_yahs_jb_no_debris_tgs.blast.out', header = F)
names(dat)<-c('bac', 'masc_seq','percent_match','overlap_length', 'num_mismatch', 'gapopen', 'tran_start', 'tran_end', 'masc_start', 'masc_end','e_val', 'bitscore')
tr_top<-dat %>%
  group_by(bac) %>%
  top_n(n = 1, bitscore)  %>%
  mutate(bac = gsub('(*)(\\.1)', '\\1', bac))

dat2<-dat %>% 
  mutate(bac = gsub('(*)(\\.[0-9])', '\\1', bac)) %>%
  left_join(L5_pos)
  

tr_top2<-tr_top %>% 
  mutate(bac = gsub('(*)(\\.[0-9])', '\\1', bac))



tr_top %>%
  left_join(L5_pos)


tr_top<-tr_top %>%
  left_join(L5_pos) %>%
  ungroup() %>%
  select(region, masc_seq, masc_start, masc_end, arm)

chr.dat2<-chr.dat %>%
  filter(masc.seq %in% tr_top2$masc_seq)

chr.dat3<-chr.dat2 %>%
  mutate(masc.seq = case_when(masc.seq == 'scaffold_1' ~ 'chr_3',
                              masc.seq == 'scaffold_2' ~ 'chr_1',
                              masc.seq == 'scaffold_3' ~ 'chr_2',
                              .default = masc.seq))


tr_top2<-tr_top %>%
  mutate(masc_seq = case_when(masc_seq == 'scaffold_1' ~ 'chr_3',
                              masc_seq == 'scaffold_2' ~ 'chr_1',
                              masc_seq == 'scaffold_3' ~ 'chr_2',
                              .default = masc_seq))


chromoMap(list(chr.dat3), list(tr_top2), 
          data_based_color_map = T,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Masc',
          chr_color = list('grey77'),
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)


chromoMap(list(chr.dat3), list(tr_top3), 
          data_based_color_map = T,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Masc',
          chr_color = list('grey77'),
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)

tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               .default = F)) %>%
  filter(chr_match == F) %>%
  left_join(L5_pos) %>%
  left_join(dat2 %>% select(bac, masc_start, masc_end, percent_match, overlap_length))


bad_map<-tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               .default = F),
         arm_map = case_when(masc_seq == 'chr_1' ~ ifelse(masc_end/1000000 > 154, 'q', 'p'),
                             masc_seq == 'chr_2' ~ ifelse(masc_end/1000000 > 230, 'q', 'p'),
                             masc_seq == 'chr_3' ~ ifelse(masc_end/1000000 > 200, 'q', 'p'),
                             grepl('scaffold_', masc_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F)
  ) %>%
  filter(arm_match == F | chr_match == F)


good_map<-tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               .default = F),
         arm_map = case_when(masc_seq == 'chr_1' ~ ifelse(masc_end/1000000 > 154, 'q', 'p'),
                             masc_seq == 'chr_2' ~ ifelse(masc_end/1000000 > 230, 'q', 'p'),
                             masc_seq == 'chr_3' ~ ifelse(masc_end/1000000 > 200, 'q', 'p'),
                             grepl('scaffold_', masc_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F)
  ) %>%
  filter(arm_match == T & chr_match == T)

chromoMap(list(chr.dat3 %>% filter(!grepl('scaffold_', masc.seq)),
               chr.dat3 %>% filter(!grepl('scaffold_', masc.seq))),
          list(bad_map %>% select(region:arm), good_map %>% select(region:arm)), 
          data_based_color_map = T,
          ploidy = 2,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Masc',
          chr_color = list('grey77'),
          export.options = T,
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)


chromoMap(list(chr.dat3), list(good_map %>% select(region:arm)), 
          data_based_color_map = T,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Masc',
          chr_color = list('grey77'),
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)

good_bad<-tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F),
         arm_map = case_when(masc_seq == 'chr_1' ~ ifelse(masc_end/1000000 > 154, 'q', 'p'),
                             masc_seq == 'chr_2' ~ ifelse(masc_end/1000000 > 230, 'q', 'p'),
                             masc_seq == 'chr_3' ~ ifelse(masc_end/1000000 > 200, 'q', 'p'),
                             grepl('scaffold_', masc_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F))

dat3<-dat2 %>%
  mutate(masc_seq = case_when(masc_seq == 'scaffold_1' ~ 'chr_3',
                              masc_seq == 'scaffold_2' ~ 'chr_1',
                              masc_seq == 'scaffold_3' ~ 'chr_2',
                              .default = masc_seq)) %>%
  left_join(L5_pos) %>%
  select(region, masc_start, masc_end, percent_match, overlap_length, bitscore)

(marker_scatter<-left_join(good_bad, dat3) %>%
  filter(!grepl('scaffold', masc_seq)) %>%
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
        axis.text = element_text(color = 'black')))

ggsave2(plot = marker_scatter, filename = 'masc_marker_scatter.pdf',
        width = 6.5, height = 6.5, units = 'in')



bad_map2<-tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               .default = F),
         arm_map = case_when(masc_seq == 'chr_1' ~ ifelse(masc_end/1000000 > 154, 'q', 'p'),
                             masc_seq == 'chr_2' ~ ifelse(masc_end/1000000 > 230, 'q', 'p'),
                             masc_seq == 'chr_3' ~ ifelse(masc_end/1000000 > 200, 'q', 'p'),
                             grepl('scaffold_', masc_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F)
  ) %>%
  filter(arm_match == T & chr_match == T) %>%
  group_by(masc_seq, arm) %>%
  arrange(masc_start, .by_group = T) %>%
  mutate(band = as.numeric(gsub('([0-9][a-z])([0-9]+)', '\\2', region))) %>%
  mutate(band_seq = case_when(arm_map == 'p' ~ band < lag(band),
                              arm_map == 'q' ~ band > lag(band)),
         band_seq = case_when(is.na(band_seq) & arm_map == 'p' & band > lead(band) ~ T,
                              is.na(band_seq) & arm_map == 'q' & band < lead(band) ~ T,
                              .default = band_seq)
  ) %>% 
  filter(band_seq == F) %>%
  select(-band_seq) %>%
  bind_rows(bad_map)


good_map2<-tr_top2 %>%
  mutate(chr_match = case_when(gsub('chr_', '', masc_seq) == gsub('([0-9])([a-z])', '\\1', arm) ~ T,
                               .default = F),
         arm_map = case_when(masc_seq == 'chr_1' ~ ifelse(masc_end/1000000 > 154, 'q', 'p'),
                             masc_seq == 'chr_2' ~ ifelse(masc_end/1000000 > 230, 'q', 'p'),
                             masc_seq == 'chr_3' ~ ifelse(masc_end/1000000 > 200, 'q', 'p'),
                             grepl('scaffold_', masc_seq) ~ NA),
         arm_match = case_when(arm_map == gsub('([0-9])([a-z])', '\\2', arm) ~ T,
                               grepl('scaffold_', masc_seq) ~ NA,
                               .default = F)
  ) %>%
  filter(arm_match == T & chr_match == T) %>%
  group_by(masc_seq, arm) %>%
  arrange(masc_start, .by_group = T) %>%
  mutate(band = as.numeric(gsub('([0-9][a-z])([0-9]+)', '\\2', region))) %>%
  mutate(band_seq = case_when(arm_map == 'p' ~ band < lag(band),
                              arm_map == 'q' ~ band > lag(band)),
         band_seq = case_when(is.na(band_seq) & arm_map == 'p' & band > lead(band) ~ T,
                              is.na(band_seq) & arm_map == 'q' & band < lead(band) ~ T,
                              .default = band_seq)) %>%
  filter(band_seq == T)
         
chromoMap(list(chr.dat3 %>% filter(!grepl('scaffold_', masc.seq)),
               chr.dat3 %>% filter(!grepl('scaffold_', masc.seq))),
          list(bad_map2 %>% select(region:arm), good_map2 %>% select(region:arm)), 
          data_based_color_map = T,
          ploidy = 2,
          data_type = 'categorical',
          labels = T,
          label_angle = -45,
          chr_length = 10,
          label_font = 12,
          title = 'Masc',
          chr_color = list('grey77'),
          export.options = T,
          data_colors = list(viridis(n = 6, begin = 0.15, end = 0.95, option = 'turbo'))
)

good_map2 %>%
  mutate(good_bad = 'good') %>%
  bind_rows(bad_map2 %>% mutate(good_bad = 'bad')) %>% 
  left_join(dat3) %>%
  ggplot(aes(x = bitscore, y = percent_match, size = overlap_length, shape = good_bad, fill = good_bad)) +
  geom_point() +
  geom_text_repel(aes(label = region), size = 2) +
  scale_shape_manual(values = c(24, 21)) +
  scale_fill_viridis(discrete = T, option = 'cividis') +
  theme_bw() +
  theme(aspect.ratio = 1)


(marker_scatter2<-good_map2 %>%
    mutate(good_bad = 'good') %>%
    bind_rows(bad_map2 %>% mutate(good_bad = 'bad')) %>% 
    left_join(dat3) %>%
    ggplot(aes(x = bitscore, y = percent_match, size = overlap_length, shape = good_bad, fill = good_bad)) +
    geom_point() +
    geom_text_repel(aes(label = region), size = 2) +
    scale_shape_manual(values = c(24, 21)) +
    scale_fill_viridis(discrete = T, option = 'cividis') +
    theme_bw() +
    theme(aspect.ratio = 1))

ggsave2(plot = marker_scatter, filename = 'masc_marker_scatter3.pdf',
        width = 6.5, height = 6.5, units = 'in')


save.image('masc_chromomap.RData')
