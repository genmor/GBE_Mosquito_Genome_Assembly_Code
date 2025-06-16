library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(cowplot)



L5.div<-read.delim('clipboard', header = T, sep = ' ') #copy 2nd half of divsum
bf.div<-read.delim('clipboard', header = T, sep = ' ') #copy 2nd half of divsum
masc.div<-read.delim('clipboard', header = T, sep = ' ') #2nd half of divsum
genome.content<-read.delim('genome.content.csv', header = T, sep = '\t') #aggregated tabulatation with specimen ids and categories
masc.size<-1288978614
bf.size<-1239189751
L5.size<-1278732104

L5.div.long<-L5.div %>% 
  select(-X) %>%
  pivot_longer(cols = -Div, names_to = 'var', values_to = 'val') %>%
  mutate(norm = val/L5.size * 100,
         class = case_when(grepl('DNA', var) ~ 'DNA',
                           grepl('LINE', var) ~ 'LINEs',
                           grepl('Helitron', var) ~ 'Helitrons',
                           grepl('MITEs', var) ~ 'MITEs',
                           grepl('LTR', var) ~ 'LTR',
                           grepl('SINE', var) ~ 'SINEs',
                           grepl('RNA', var) ~ 'RNA',
                           var %in% c('UD.UD', 'Unknown', 'Unspecified') ~ 'Unknown',
                           grepl('ARTEFACT', var) ~ 'Artefact',
                           var %in% c('Satellite', 'Simple_repeat') ~ 'Simple repeat',
                           grepl('Penelope', var) ~ 'PLE',
                           grepl('PLE', var) ~ 'PLE'))


L5.land<-L5.div.long %>%
  filter(class != 'Artefact') %>%
  group_by(Div, class) %>%
  summarise(norm = sum(norm)) %>%
  mutate(class = factor(class, levels = c('LTR', 'LINEs',
                                          'SINEs', 'PLE', 'DNA', 'MITEs',
                                          'Helitrons', 'RNA', 'Simple repeat', 'Unknown', 'Artefact')))

bf.div.long<-bf.div %>% 
  select(-X) %>%
  pivot_longer(cols = -Div, names_to = 'var', values_to = 'val') %>%
  mutate(norm = val/L5.size * 100,
         class = case_when(grepl('DNA', var) ~ 'DNA',
                           grepl('LINE', var) ~ 'LINEs',
                           grepl('Helitron', var) ~ 'Helitrons',
                           grepl('MITEs', var) ~ 'MITEs',
                           grepl('LTR', var) ~ 'LTR',
                           grepl('SINE', var) ~ 'SINEs',
                           grepl('RNA', var) ~ 'RNA',
                           var %in% c('UD.UD', 'Unknown', 'Unspecified') ~ 'Unknown',
                           grepl('ARTEFACT', var) ~ 'Artefact',
                           var %in% c('Satellite', 'Simple_repeat') ~ 'Simple repeat',
                           grepl('Penelope', var) ~ 'PLE',
                           grepl('PLE', var) ~ 'PLE'))


bf.land<-bf.div.long %>% 
  filter(class != 'Artefact') %>%
  group_by(Div, class) %>%
  summarise(norm = sum(norm)) %>%
  mutate(class = factor(class, levels = c('LTR', 'LINEs',
                                          'SINEs', 'PLE', 'DNA', 'MITEs',
                                          'Helitrons', 'RNA', 'Simple repeat', 'Unknown', 'Artefact')))

masc.div.long<-masc.div %>% 
  select(-X) %>%
  pivot_longer(cols = -Div, names_to = 'var', values_to = 'val') %>%
  mutate(norm = val/L5.size * 100,
         class = case_when(grepl('DNA', var) ~ 'DNA',
                           grepl('LINE', var) ~ 'LINEs',
                           grepl('Helitron', var) ~ 'Helitrons',
                           grepl('MITEs', var) ~ 'MITEs',
                           grepl('LTR', var) ~ 'LTR',
                           grepl('SINE', var) ~ 'SINEs',
                           grepl('RNA', var) ~ 'RNA',
                           var %in% c('UD.UD', 'Unknown', 'Unspecified') ~ 'Unknown',
                           grepl('ARTEFACT', var) ~ 'Artefact',
                           var %in% c('Satellite', 'Simple_repeat') ~ 'Simple repeat',
                           grepl('Penelope', var) ~ 'PLE',
                           grepl('PLE', var) ~ 'PLE'))


masc.land<-masc.div.long %>%
  filter(class != 'Artefact') %>%
  group_by(Div, class) %>%
  summarise(norm = sum(norm)) %>%
  mutate(class = factor(class, levels = c('LTR', 'LINEs',
                                          'SINEs', 'PLE', 'DNA', 'MITEs',
                                          'Helitrons', 'RNA', 'Simple repeat', 'Unknown', 'Artefact')))
  
  

(L5.land.plot<-ggplot(L5.land, aes(x = Div, y = norm, fill = class)) +
  geom_bar(position = 'stack', stat = 'identity', color = 'black', lwd = 0.3) +
  scale_fill_viridis(discrete = T, '') +
  labs(x = 'Kimura substitution level (CpG adjusted)', y = 'Percent of genome') +
  ggtitle('A. AaegL5') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        #axis.title.x = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 10),
        axis.text = element_text(color = 'black')
        )
)


(bf.land.plot<-ggplot(bf.land, aes(x = Div, y = norm, fill = class)) +
    geom_bar(position = 'stack', stat = 'identity', color = 'black', lwd = 0.3) +
    scale_fill_viridis(discrete = T, '') +
    labs(x = 'Kimura substitution level (CpG adjusted)', y = 'Percent of genome') +
    ggtitle('B. Aedes aegpyti formosus (Aaf)') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          #axis.title.x = element_blank(),
          legend.position = 'none',
          legend.title = element_text(size = 10),
          axis.text = element_text(color = 'black')
    )
)

(masc.land.plot<-ggplot(masc.land, aes(x = Div, y = norm, fill = class)) +
    geom_bar(position = 'stack', stat = 'identity', color = 'black', lwd = 0.3) +
    scale_fill_viridis(discrete = T, '') +
    labs(x = 'Kimura substitution level (CpG adjusted)', y = 'Percent of genome') +
    ggtitle('C. Aedes mascarensis (Am)') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::pretty_breaks(n = 8)) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          legend.position = 'none',
          legend.title = element_text(size = 10),
          axis.text = element_text(color = 'black')
    )
)



genome.content2<-genome.content %>%
  select(-percent_masked) %>%
  group_by(specimen, class) %>%
  summarise(count = sum(count),
            bp = sum(bpMasked)) %>%
  ungroup %>%
  mutate(percent = case_when(specimen == 'AaegL5' ~ bp/L5.size,
                             specimen == 'Bf05' ~ bp/bf.size,
                             specimen == 'MascCH02' ~ bp/masc.size),
         percent = percent * 100) %>%
  # filter(class != 'ARTEFACT') %>%
  mutate(class = factor(class, levels = c('LTR', 'LINEs',
                                 'SINEs', 'PLE', 'DNA', 'MITEs',
                                 'Helitrons', 'RNA', 'Simple repeat', 'Unknown', 'Artefact')
                        )
         )


unmasked<-genome.content2 %>%
  group_by(specimen) %>%
  summarise(bp = sum(bp)) %>%
  ungroup %>%
  mutate(size = case_when(specimen == 'AaegL5' ~ L5.size,
                          specimen == 'Bf05' ~ bf.size,
                          specimen == 'MascCH02' ~ masc.size),
         unmasked_bp = size - bp,
         percent = unmasked_bp/size * 100,
         class = 'Unmasked DNA') %>%
  select(specimen, unmasked_bp, percent, class) %>%
  rename(bp = unmasked_bp)


genome.content3<-genome.content2 %>%
  bind_rows(unmasked)

labs<-c('LTR', 'LINEs', 'SINEs', 'PLE', 'DNA','MITEs','Helitrons',
        'RNA', 'Simple repeat', 'Unknown', 'Unmasked DNA')

L5.percent<-genome.content3 %>%
  filter(specimen == 'AaegL5') %>%
  select(-count) %>%
  mutate(ymax = cumsum(percent),
         ymin = c(0, head(ymax, n = -1)),
         class = factor(class, levels = c('LTR', 'LINEs', 'SINEs', 'PLE', 'DNA', 'MITEs', 'Helitrons',
                                          'RNA', 'Simple repeat', 'Unknown'))
         )

L5.perc.labs<-L5.percent %>%
  pull(percent, name = class)
names(L5.perc.labs)[length(L5.perc.labs)]<-'Unmasked DNA'
  
(L5.donut<-ggplot(L5.percent, aes(x = 9, y = percent, fill = class)) +
  geom_col(color = 'black', lwd = 0.3) +
  coord_polar(theta = 'y') +
  xlim(c(1, 10)) +
  scale_fill_viridis(option = 'viridis', discrete = T, na.value = 'grey77',
                     labels = paste(labs, round(L5.perc.labs[labs],1), '%'), 'Genomic content\nof AaegL5') +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.box.margin = margin(-60, -60, -60, -60),
        legend.box.spacing = unit(1.75, 'cm'))
)

bf.percent<-genome.content3 %>%
  filter(specimen == 'Bf05') %>%
  select(-count) %>%
  mutate(ymax = cumsum(percent),
         ymin = c(0, head(ymax, n = -1)),
         class = factor(class, levels = c('LTR', 'LINEs', 'SINEs', 'PLE', 'DNA', 'MITEs', 'Helitrons',
                                          'RNA', 'Simple repeat', 'Unknown'))
  )

bf.perc.labs<-bf.percent %>%
  pull(percent, name = class)
names(bf.perc.labs)[length(bf.perc.labs)]<-'Unmasked DNA'

(bf.donut<-ggplot(bf.percent, aes(x = 9, y = percent, fill = class)) +
    geom_col(color = 'black', lwd= 0.3) +
    coord_polar(theta = 'y') +
    xlim(c(1, 10)) +
    scale_fill_viridis(option = 'viridis', discrete = T, na.value = 'grey77',
                       labels = paste(labs, round(bf.perc.labs[labs],1), '%'), 'Genomic content\nof Aaf Bf05') +
    theme_void() +
    theme(aspect.ratio = 1,
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.box.margin = margin(-60, -60, -60, -60),
          legend.box.spacing = unit(1.75, 'cm'))
)

masc.percent<-genome.content3 %>%
  filter(specimen == 'MascCH02') %>%
  select(-count) %>%
  mutate(ymax = cumsum(percent),
         ymin = c(0, head(ymax, n = -1)),
         class = factor(class, levels = c('LTR', 'LINEs', 'SINEs', 'PLE', 'DNA', 'MITEs', 'Helitrons',
                                          'RNA', 'Simple repeat', 'Unknown'))
  )

masc.perc.labs<-masc.percent %>%
  pull(percent, name = class)
names(masc.perc.labs)[length(masc.perc.labs)]<-'Unmasked DNA'

(masc.donut<-ggplot(masc.percent, aes(x = 9, y = percent, fill = class)) +
    geom_col(color = 'black', lwd = 0.3) +
    coord_polar(theta = 'y') +
    xlim(c(1, 10)) +
    scale_fill_viridis(option = 'viridis', discrete = T, na.value = 'grey77',
                       labels = paste(labs, round(masc.perc.labs[labs],1), '%'), 'Genomic content\nof Am MascCH02') +
    theme_void() +
    theme(aspect.ratio = 1,
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.box.margin = margin(-60, -60, -60, -60),
          legend.box.spacing = unit(1.75, 'cm'))
)

(L5.all<-L5.land.plot + inset_element(L5.donut, left = 0.3, bottom = 0.2, right = 1, top = 1, align_to = 'plot'))
(bf.all<-bf.land.plot + inset_element(bf.donut, left = 0.3, bottom = 0.2, right = 1, top = 1, align_to = 'plot'))
(masc.all<-masc.land.plot + inset_element(masc.donut, left = 0.3, bottom = 0.2, right = 1, top = 1, align_to = 'plot'))

(panel<-L5.all + bf.all + masc.all + plot_layout(ncol = 1) & theme(legend.text = element_text(size = 9, color = 'black'), legend.key.size = unit(6, 'pt')))
# ggsave2(file = 'rep.panel2.pdf', plot = panel, width = 6, height = 8, unit = 'in')

percs<-list(masc.perc.labs, L5.perc.labs, bf.perc.labs)
spp<-c('AaegL5', 'Bf05', 'MascCH02')
names(percs)<-spp
rep.barplot<-lapply(percs, function(x) {
  z<-data.frame(percent = x)
  z$class<-rownames(z)
  rownames(z)<-NULL
  z
}) %>%
  bind_rows(., .id = 'specimen') %>%
  mutate(class = factor(class, levels = c('LTR', 'LINEs', 'SINEs', 'PLE', 'DNA', 'MITEs', 'Helitrons',
                                          'RNA', 'Simple repeat', 'Unknown'))
  ) %>%
  ggplot(aes(x = specimen, y = percent, fill = class)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
  scale_fill_viridis(option = 'viridis', discrete = T, na.value = 'grey77',
                     labels = function(zz) {zz[is.na(zz)]<-'Non-repetitive DNA'; zz},
                     'Repeat class') +
  geom_text(aes(specimen, group = class, label = sprintf('%2.1f', percent)), 
            position = position_dodge(width = 0.9), angle = 35, size = 3) +
  ggtitle('D. Comparison of genomic content') +
  labs(y = 'Percent of genome', x = 'Assembly') +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = 'right',
        legend.title = element_text(size = 10),
        axis.text = element_text(color = 'black'),
        )
  

(panel2<-L5.land.plot + bf.land.plot + masc.land.plot  + rep.barplot + 
    plot_layout(nrow = 2, guide = 'collect') + theme(axis.title = element_text(size = 10),
                                                     axis.text = element_text(size = 9),
                                                     title = element_text(size = 10)))
# ggsave2(file = 'rep.panel3.pdf', plot = panel2, width = 6, height = 8, unit = 'in')

# save.image('repeat_content.RData')
