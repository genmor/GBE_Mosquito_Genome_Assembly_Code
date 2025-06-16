library(GENESPACE)
library(ggplot2)
### set up paths
setwd("/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae9/genespace")
wd<-getwd()
mcscan<-"/work/soghigian_lab/apps/MCScanX/MCScanX-master/"
genomeRepo<-"/work/soghigian_lab/gen.morinaga/compare_genomes5/comparative_analysis/compare_genomes3/culicidae9/genespace/genome_repo/"
cores<-25
### create bed files from gff and fasta

p1<-parse_annotations(rawGenomeRepo = genomeRepo,
	genomeDirs = c('Aedes_aegypti_L5',
	'Aedes_aegypti_formosus',
	'Aedes_mascarensis',
	'Armigeres_subalbatus',
	'Anopheles_cruzii',
'Anopheles_darlingi',
'Anopheles_gambiae',
'Anopheles_ziemanni',
'Culex_pipiens_pallens',
'Culex_quinquefasciatus',
'Phlebotomus_papatasi',
'Sabethes_cyaneus'
 ),
	 genomeIDs = c('Aedes_aegypti_L5',
	'Aedes_aegypti_formosus',
	'Aedes_mascarensis',
	'Armigeres_subalbatus',
	'Anopheles_cruzii',
'Anopheles_darlingi',
'Anopheles_gambiae',
'Anopheles_ziemanni',
'Culex_pipiens_pallens',
'Culex_quinquefasciatus',
'Phlebotomus_papatasi',
'Sabethes_cyaneus'
 ),
	 presets = 'none',
	 headerEntryIndex = 1,
	 genespaceWd = wd,
	 gffIdColumn = 'ID'
 )
 
 
### genomes annotated using braker need to be treated differently...

### initialize paths for genespace and run

gpar<-init_genespace(wd = wd, path2mcscanx = mcscan, nCores = cores)
out<-run_genespace(gpar)


save.image('genespace_run.RData')
### create custom riparian plot from output

genetheme11<-theme(panel.background = element_rect(fill = 'white')) #sets the background to be white
cust.plot2<-plot_riparian(gsParam = out,
 	chrFill = 'grey88',
 	refGenome = 'Aedes_aegypti_L5',
 	genomeIDs = c('Aedes_aegypti_L5',
 'Aedes_aegypti_formosus',
 'Aedes_mascarensis',
 'Armigeres_subalbatus',
 'Sabethes_cyaneus',
 'Culex_pipiens_pallens',
 'Culex_quinquefasciatus',
 'Anopheles_cruzii',
 'Anopheles_darlingi',
 'Anopheles_gambiae',
 'Anopheles_ziemanni',
 'Phlebotomus_papatasi'
 ),
 	useRegions = F,
 	useOrder = F,
 	minChrLen2plot = 2e6,
 	addThemes = genetheme11,
 	forceRecalcBlocks = F
 )

save.image('genespace_run.RData')
