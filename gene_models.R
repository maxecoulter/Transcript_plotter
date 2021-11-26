library(dplyr)

source('genemodel_adjusted.r')

input <- 'genemodel_input.csv'
transcripts <- read.csv(input)

domain_input <- 'genemodel_input_domain.csv'
domains <- read.csv(domain_input)

snpeff_input <- read.csv("snpeff_out.csv")

#Working on BaRT2v18chr3HG123260


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123260')

gene.transcript.model.plot(model=gene, gene_start=35765878, gene_bpstop=35767144, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123260')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35765878, gene_bpstop=35767144, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123260.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123260.1')

genemodel.plot(model=transcript, start=35765878, bpstop=35767144, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123260.1')
genemodel.plot_domain(model=t_domains, start=35765878, bpstop=35767144, orientation='forward')

snpeff <- snpeff_input %>% filter(Transcript == 'BaRT2v18chr3HG123260.1')
plot_snpeff <- function(snpeff=snpeff){
for (i in 1:nrow(snpeff)){
  start <- as.numeric(snpeff[i,"position_start"])
  end <- as.numeric(snpeff[i,"position_end"])
  label <- snpeff[i,"label"]
  depth <- snpeff[i,"depth"]
  colour <- snpeff[i,"colour"]
  mutation.plot(start, end, text=label, col='black', drop=depth, haplotypes=c(colour))
  
}
}
plot_snpeff(snpeff=snpeff)

mutation.plot(35766007, 35766007, text='', col='black', drop=-0.3527546437430802, haplotypes=c('orange'))

mutation.plot(35766138, 35766138, text='', col='black', drop=-0.12936699416704553, haplotypes=c('blue'))

mutation.plot(35766222, 35766222, text='', col='black', drop=-0.17521245583903167, haplotypes=c('blue'))

mutation.plot(35766597, 35766597, text='', col='black', drop=-0.17846992983838725, haplotypes=c('blue'))

mutation.plot(35766622, 35766622, text='', col='black', drop=-0.15273524595601415, haplotypes=c('blue'))

mutation.plot(35766680, 35766680, text='', col='black', drop=-0.37436741703148535, haplotypes=c('orange'))

mutation.plot(35766696, 35766696, text='', col='black', drop=-0.1933690687020369, haplotypes=c('blue'))

mutation.plot(35766755, 35766755, text='', col='black', drop=-0.33597200944004746, haplotypes=c('orange'))

#Working on BaRT2v18chr3HG123040


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123040')

gene.transcript.model.plot(model=gene, gene_start=33873268, gene_bpstop=33876381, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123040')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33873268, gene_bpstop=33876381, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123040.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123040.1')

genemodel.plot(model=transcript, start=33873268, bpstop=33876381, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123040.1')
genemodel.plot_domain(model=t_domains, start=33873268, bpstop=33876381, orientation='forward')

mutation.plot(33873335, 33873335, text='', col='black', drop=-0.1383768535170063, haplotypes=c('blue'))

mutation.plot(33873337, 33873337, text='', col='black', drop=-0.10380734527912723, haplotypes=c('blue'))

mutation.plot(33873343, 33873343, text='', col='black', drop=-0.10236515822925263, haplotypes=c('blue'))

mutation.plot(33873492, 33873492, text='', col='black', drop=-0.3754205153670463, haplotypes=c('orange'))

mutation.plot(33873694, 33873694, text='', col='black', drop=-0.12408252190258606, haplotypes=c('blue'))

mutation.plot(33873715, 33873715, text='', col='black', drop=-0.3885768352676785, haplotypes=c('orange'))

mutation.plot(33873802, 33873802, text='', col='black', drop=-0.35148306339481256, haplotypes=c('orange'))

mutation.plot(33873988, 33873988, text='', col='black', drop=-0.19767895871682045, haplotypes=c('blue'))

mutation.plot(33876153, 33876153, text='', col='black', drop=-0.344965806675852, haplotypes=c('orange'))

mutation.plot(33876216, 33876217, text='', col='black', drop=-0.14204932261921163, haplotypes=c('blue'))

mutation.plot(33876241, 33876241, text='', col='black', drop=-0.169397291435684, haplotypes=c('blue'))

mutation.plot(33876282, 33876282, text='', col='black', drop=-0.12718008716159412, haplotypes=c('blue'))

mutation.plot(33876372, 33876373, text='', col='black', drop=-0.10380131050935766, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123040.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123040.2')

genemodel.plot(model=transcript, start=33873273, bpstop=33876380, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123040.2')
genemodel.plot_domain(model=t_domains, start=33873273, bpstop=33876380, orientation='forward')

mutation.plot(33873335, 33873335, text='', col='black', drop=-0.16471847624672842, haplotypes=c('blue'))

mutation.plot(33873337, 33873337, text='', col='black', drop=-0.10083314739229905, haplotypes=c('blue'))

mutation.plot(33873343, 33873343, text='', col='black', drop=-0.14282207325266827, haplotypes=c('blue'))

mutation.plot(33873492, 33873492, text='', col='black', drop=-0.3528928104823901, haplotypes=c('orange'))

mutation.plot(33873694, 33873694, text='', col='black', drop=-0.18004679439905752, haplotypes=c('blue'))

mutation.plot(33873715, 33873715, text='', col='black', drop=-0.3420562573280504, haplotypes=c('orange'))

mutation.plot(33873802, 33873802, text='', col='black', drop=-0.3409843409392511, haplotypes=c('orange'))

mutation.plot(33873988, 33873988, text='', col='black', drop=-0.13453145629742735, haplotypes=c('blue'))

mutation.plot(33875039, 33875039, text='', col='black', drop=-0.16152831426560638, haplotypes=c('blue'))

mutation.plot(33876153, 33876153, text='', col='black', drop=-0.18270072844805332, haplotypes=c('blue'))

mutation.plot(33876216, 33876217, text='', col='black', drop=-0.18070721957536195, haplotypes=c('blue'))

mutation.plot(33876241, 33876241, text='', col='black', drop=-0.19066980309402937, haplotypes=c('blue'))

mutation.plot(33876282, 33876282, text='', col='black', drop=-0.14328017559397874, haplotypes=c('blue'))

mutation.plot(33876372, 33876373, text='', col='black', drop=-0.12453603601030094, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122910


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122910')

gene.transcript.model.plot(model=gene, gene_start=33174328, gene_bpstop=33177081, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122910')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33174328, gene_bpstop=33177081, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG122910.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122910.1')

genemodel.plot(model=transcript, start=33174328, bpstop=33177081, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122910.1')
genemodel.plot_domain(model=t_domains, start=33174328, bpstop=33177081, orientation='forward')

mutation.plot(33174355, 33174355, text='', col='black', drop=-0.1368624588862372, haplotypes=c('blue'))

mutation.plot(33174373, 33174373, text='', col='black', drop=-0.11792768747519018, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123210


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123210')

gene.transcript.model.plot(model=gene, gene_start=35324956, gene_bpstop=35329217, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123210')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35324956, gene_bpstop=35329217, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123210.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123210.1')

genemodel.plot(model=transcript, start=35324956, bpstop=35329217, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123210.1

#Working on BaRT2v18chr3HG123200


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123200')

gene.transcript.model.plot(model=gene, gene_start=35324086, gene_bpstop=35324557, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123200')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35324086, gene_bpstop=35324557, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123200.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123200.1')

genemodel.plot(model=transcript, start=35324086, bpstop=35324557, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123200.1

#Working on BaRT2v18chr3HG123130


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123130')

gene.transcript.model.plot(model=gene, gene_start=35119001, gene_bpstop=35120274, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123130')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35119001, gene_bpstop=35120274, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123130.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123130.1')

genemodel.plot(model=transcript, start=35119001, bpstop=35120274, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123130.1')
genemodel.plot_domain(model=t_domains, start=35119001, bpstop=35120274, orientation='forward')

mutation.plot(35119035, 35119035, text='', col='black', drop=-0.1209929926411925, haplotypes=c('blue'))

mutation.plot(35120211, 35120211, text='', col='black', drop=-0.19229718743622487, haplotypes=c('blue'))

mutation.plot(35120238, 35120239, text='', col='black', drop=-0.12388123609281225, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123130.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123130.2')

genemodel.plot(model=transcript, start=35119003, bpstop=35120226, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123130.2')
genemodel.plot_domain(model=t_domains, start=35119003, bpstop=35120226, orientation='forward')

mutation.plot(35119035, 35119035, text='', col='black', drop=-0.1160652342624112, haplotypes=c('blue'))

mutation.plot(35120211, 35120211, text='', col='black', drop=-0.12810227440735963, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123350


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123350')

gene.transcript.model.plot(model=gene, gene_start=36110458, gene_bpstop=36112988, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123350')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36110458, gene_bpstop=36112988, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123350.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123350.1')

genemodel.plot(model=transcript, start=36110458, bpstop=36112988, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123350.1')
genemodel.plot_domain(model=t_domains, start=36110458, bpstop=36112988, orientation='forward')

mutation.plot(36112349, 36112349, text='', col='black', drop=-0.1461108095665773, haplotypes=c('blue'))

mutation.plot(36112858, 36112858, text='', col='black', drop=-0.12481416724724535, haplotypes=c('blue'))

mutation.plot(36112925, 36112925, text='', col='black', drop=-0.14487562420746286, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123190


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123190')

gene.transcript.model.plot(model=gene, gene_start=35319604, gene_bpstop=35320180, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123190')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35319604, gene_bpstop=35320180, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123190.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123190.1')

genemodel.plot(model=transcript, start=35319604, bpstop=35320175, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123190.1')
genemodel.plot_domain(model=t_domains, start=35319604, bpstop=35320175, orientation='forward')

mutation.plot(35319772, 35319772, text='', col='black', drop=-0.14941444428993533, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123190.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123190.2')

genemodel.plot(model=transcript, start=35319607, bpstop=35320180, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123190.2')
genemodel.plot_domain(model=t_domains, start=35319607, bpstop=35320180, orientation='forward')

mutation.plot(35319772, 35319772, text='', col='black', drop=-0.10906857476359955, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123190.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123190.3')

genemodel.plot(model=transcript, start=35319607, bpstop=35320175, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123190.3')
genemodel.plot_domain(model=t_domains, start=35319607, bpstop=35320175, orientation='forward')

#Working on BaRT2v18chr3HG123410


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123410')

gene.transcript.model.plot(model=gene, gene_start=36310092, gene_bpstop=36311949, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123410')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36310092, gene_bpstop=36311949, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123410.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123410.1')

genemodel.plot(model=transcript, start=36310092, bpstop=36311941, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123410.1')
genemodel.plot_domain(model=t_domains, start=36310092, bpstop=36311941, orientation='forward')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.39389288095821695, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.10753374204553875, haplotypes=c('blue'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.12331906178620525, haplotypes=c('blue'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.1432287544258864, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.18273328689001664, haplotypes=c('blue'))

mutation.plot(36311012, 36311012, text='', col='black', drop=-0.18173318076349793, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123410.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123410.2')

genemodel.plot(model=transcript, start=36310096, bpstop=36311945, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123410.2')
genemodel.plot_domain(model=t_domains, start=36310096, bpstop=36311945, orientation='forward')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.392268496676563, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.15468973950184628, haplotypes=c('blue'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.13492133769057144, haplotypes=c('blue'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.14705814041051227, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.15980964300077044, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123410.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123410.3')

genemodel.plot(model=transcript, start=36310096, bpstop=36311945, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123410.3')
genemodel.plot_domain(model=t_domains, start=36310096, bpstop=36311945, orientation='forward')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.3492780560556006, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.11621145870262886, haplotypes=c('blue'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.19585544999623283, haplotypes=c('blue'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.12110597884002611, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.1314262002600664, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123410.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123410.4')

genemodel.plot(model=transcript, start=36310096, bpstop=36311949, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123410.4')
genemodel.plot_domain(model=t_domains, start=36310096, bpstop=36311949, orientation='forward')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.3596594158306172, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.1604778600464295, haplotypes=c('blue'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.11640016040237791, haplotypes=c('blue'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.11325498806349417, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.19812879069661743, haplotypes=c('blue'))

mutation.plot(36311012, 36311012, text='', col='black', drop=-0.14319930372930714, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123410.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123410.5')

genemodel.plot(model=transcript, start=36310096, bpstop=36311941, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123410.5')
genemodel.plot_domain(model=t_domains, start=36310096, bpstop=36311941, orientation='forward')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.3677639429819012, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.15187021793523073, haplotypes=c('blue'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.12177450074281204, haplotypes=c('blue'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.1889752699766528, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.12254510380344759, haplotypes=c('blue'))

mutation.plot(36311012, 36311012, text='', col='black', drop=-0.1808613284095827, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123460


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123460')

gene.transcript.model.plot(model=gene, gene_start=36511900, gene_bpstop=36512548, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123460')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36511900, gene_bpstop=36512548, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123460.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123460.1')

genemodel.plot(model=transcript, start=36511900, bpstop=36512548, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123460.1')
genemodel.plot_domain(model=t_domains, start=36511900, bpstop=36512548, orientation='forward')

#Working on BaRT2v18chr3HG123080


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123080')

gene.transcript.model.plot(model=gene, gene_start=34347021, gene_bpstop=34357701, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123080')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=34347021, gene_bpstop=34357701, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123080.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123080.1')

genemodel.plot(model=transcript, start=34347021, bpstop=34357683, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123080.1')
genemodel.plot_domain(model=t_domains, start=34347021, bpstop=34357683, orientation='forward')

mutation.plot(34347332, 34347332, text='', col='black', drop=-0.14116002687775792, haplotypes=c('blue'))

mutation.plot(34356216, 34356216, text='', col='black', drop=-0.18332530261856592, haplotypes=c('blue'))

mutation.plot(34356256, 34356256, text='', col='black', drop=-0.3259346652641342, haplotypes=c('orange'))

mutation.plot(34356257, 34356257, text='', col='black', drop=-0.3064274304567681, haplotypes=c('orange'))

mutation.plot(34356694, 34356694, text='', col='black', drop=-0.3556298404217278, haplotypes=c('orange'))

mutation.plot(34356712, 34356712, text='', col='black', drop=-0.33659976911767253, haplotypes=c('orange'))

mutation.plot(34356830, 34356830, text='', col='black', drop=-0.1902997423526742, haplotypes=c('blue'))

mutation.plot(34357135, 34357135, text='', col='black', drop=-0.3631410410289815, haplotypes=c('orange'))

mutation.plot(34357465, 34357465, text='', col='black', drop=-0.1963979189085287, haplotypes=c('blue'))

mutation.plot(34357527, 34357527, text='', col='black', drop=-0.12029778946346642, haplotypes=c('blue'))

mutation.plot(34357544, 34357544, text='', col='black', drop=-0.1315260951672301, haplotypes=c('blue'))

mutation.plot(34357622, 34357622, text='', col='black', drop=-0.17465475561605726, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123080.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123080.2')

genemodel.plot(model=transcript, start=34347021, bpstop=34357701, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123080.2')
genemodel.plot_domain(model=t_domains, start=34347021, bpstop=34357701, orientation='forward')

mutation.plot(34347332, 34347332, text='', col='black', drop=-0.1354997551008787, haplotypes=c('blue'))

mutation.plot(34347403, 34347403, text='', col='black', drop=-0.17621888165805583, haplotypes=c('blue'))

mutation.plot(34347404, 34347404, text='', col='black', drop=-0.19222587922842513, haplotypes=c('blue'))

mutation.plot(34347457, 34347457, text='', col='black', drop=-0.13796572002249038, haplotypes=c('blue'))

mutation.plot(34347624, 34347624, text='', col='black', drop=-0.12770276242436865, haplotypes=c('blue'))

mutation.plot(34347967, 34347967, text='', col='black', drop=-0.3249592052495758, haplotypes=c('orange'))

mutation.plot(34348002, 34348002, text='', col='black', drop=-0.17030412295913305, haplotypes=c('blue'))

mutation.plot(34356216, 34356216, text='', col='black', drop=-0.1569655362831333, haplotypes=c('blue'))

mutation.plot(34356256, 34356256, text='', col='black', drop=-0.30443090695971264, haplotypes=c('orange'))

mutation.plot(34356257, 34356257, text='', col='black', drop=-0.3993175528148876, haplotypes=c('orange'))

mutation.plot(34356694, 34356694, text='', col='black', drop=-0.3854813072828007, haplotypes=c('orange'))

mutation.plot(34356712, 34356712, text='', col='black', drop=-0.3869125099511583, haplotypes=c('orange'))

mutation.plot(34356830, 34356830, text='', col='black', drop=-0.138150181943027, haplotypes=c('blue'))

mutation.plot(34357135, 34357135, text='', col='black', drop=-0.39896285262509823, haplotypes=c('orange'))

mutation.plot(34357465, 34357465, text='', col='black', drop=-0.16556763506181807, haplotypes=c('blue'))

mutation.plot(34357527, 34357527, text='', col='black', drop=-0.11568350832456553, haplotypes=c('blue'))

mutation.plot(34357544, 34357544, text='', col='black', drop=-0.12660135001703085, haplotypes=c('blue'))

mutation.plot(34357622, 34357622, text='', col='black', drop=-0.13212421995442702, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123080.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123080.3')

genemodel.plot(model=transcript, start=34347021, bpstop=34357683, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123080.3')
genemodel.plot_domain(model=t_domains, start=34347021, bpstop=34357683, orientation='forward')

mutation.plot(34347332, 34347332, text='', col='black', drop=-0.19273643179115657, haplotypes=c('blue'))

mutation.plot(34347403, 34347403, text='', col='black', drop=-0.19871006989788076, haplotypes=c('blue'))

mutation.plot(34347404, 34347404, text='', col='black', drop=-0.1795747365504116, haplotypes=c('blue'))

mutation.plot(34347457, 34347457, text='', col='black', drop=-0.11863515497616486, haplotypes=c('blue'))

mutation.plot(34347624, 34347624, text='', col='black', drop=-0.16295121931264148, haplotypes=c('blue'))

mutation.plot(34347967, 34347967, text='', col='black', drop=-0.39111808386634705, haplotypes=c('orange'))

mutation.plot(34348002, 34348002, text='', col='black', drop=-0.15072309432086, haplotypes=c('blue'))

mutation.plot(34356216, 34356216, text='', col='black', drop=-0.12810841873478643, haplotypes=c('blue'))

mutation.plot(34356256, 34356256, text='', col='black', drop=-0.3903784132510583, haplotypes=c('orange'))

mutation.plot(34356257, 34356257, text='', col='black', drop=-0.35854780825876753, haplotypes=c('orange'))

mutation.plot(34356694, 34356694, text='', col='black', drop=-0.3212396731349585, haplotypes=c('orange'))

mutation.plot(34356712, 34356712, text='', col='black', drop=-0.3539000396784329, haplotypes=c('orange'))

mutation.plot(34356830, 34356830, text='', col='black', drop=-0.14469017412865595, haplotypes=c('blue'))

mutation.plot(34357135, 34357135, text='', col='black', drop=-0.18516224704688433, haplotypes=c('blue'))

mutation.plot(34357465, 34357465, text='', col='black', drop=-0.10166036830691051, haplotypes=c('blue'))

mutation.plot(34357527, 34357527, text='', col='black', drop=-0.12912620934210692, haplotypes=c('blue'))

mutation.plot(34357544, 34357544, text='', col='black', drop=-0.10488868952861588, haplotypes=c('blue'))

mutation.plot(34357622, 34357622, text='', col='black', drop=-0.17460084185305225, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123080.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123080.4')

genemodel.plot(model=transcript, start=34347069, bpstop=34357699, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123080.4')
genemodel.plot_domain(model=t_domains, start=34347069, bpstop=34357699, orientation='forward')

mutation.plot(34347332, 34347332, text='', col='black', drop=-0.11212205496258071, haplotypes=c('blue'))

mutation.plot(34347624, 34347624, text='', col='black', drop=-0.15717270804181088, haplotypes=c('blue'))

mutation.plot(34347967, 34347967, text='', col='black', drop=-0.3809859278409279, haplotypes=c('orange'))

mutation.plot(34348002, 34348002, text='', col='black', drop=-0.10633776625747224, haplotypes=c('blue'))

mutation.plot(34356216, 34356216, text='', col='black', drop=-0.10836906420326006, haplotypes=c('blue'))

mutation.plot(34356256, 34356256, text='', col='black', drop=-0.37612389627410336, haplotypes=c('orange'))

mutation.plot(34356257, 34356257, text='', col='black', drop=-0.35450106654070107, haplotypes=c('orange'))

mutation.plot(34356694, 34356694, text='', col='black', drop=-0.30459662009269134, haplotypes=c('orange'))

mutation.plot(34356712, 34356712, text='', col='black', drop=-0.3320853941596997, haplotypes=c('orange'))

mutation.plot(34356830, 34356830, text='', col='black', drop=-0.19261119984316152, haplotypes=c('blue'))

mutation.plot(34357135, 34357135, text='', col='black', drop=-0.36596365637252504, haplotypes=c('orange'))

mutation.plot(34357465, 34357465, text='', col='black', drop=-0.11842697503438726, haplotypes=c('blue'))

mutation.plot(34357527, 34357527, text='', col='black', drop=-0.14465288408538557, haplotypes=c('blue'))

mutation.plot(34357544, 34357544, text='', col='black', drop=-0.1887818024611459, haplotypes=c('blue'))

mutation.plot(34357622, 34357622, text='', col='black', drop=-0.1514340298138326, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123500


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123500')

gene.transcript.model.plot(model=gene, gene_start=36880422, gene_bpstop=36890887, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123500')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36880422, gene_bpstop=36890887, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123500.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123500.1')

genemodel.plot(model=transcript, start=36880422, bpstop=36887250, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123500.1')
genemodel.plot_domain(model=t_domains, start=36880422, bpstop=36887250, orientation='forward')

#Working on BaRT2v18chr3HG123500.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123500.2')

genemodel.plot(model=transcript, start=36880915, bpstop=36885874, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123500.2')
genemodel.plot_domain(model=t_domains, start=36880915, bpstop=36885874, orientation='forward')

#Working on BaRT2v18chr3HG123500.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123500.3')

genemodel.plot(model=transcript, start=36885068, bpstop=36890887, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123500.3')
genemodel.plot_domain(model=t_domains, start=36885068, bpstop=36890887, orientation='forward')

mutation.plot(36890032, 36890032, text='', col='black', drop=-0.17434637259834784, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123500.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123500.4')

genemodel.plot(model=transcript, start=36885068, bpstop=36890887, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123500.4')
genemodel.plot_domain(model=t_domains, start=36885068, bpstop=36890887, orientation='forward')

mutation.plot(36890032, 36890032, text='', col='black', drop=-0.1403130000448546, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123310


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123310')

gene.transcript.model.plot(model=gene, gene_start=36056239, gene_bpstop=36056821, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123310')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36056239, gene_bpstop=36056821, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123310.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123310.1')

genemodel.plot(model=transcript, start=36056239, bpstop=36056821, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123310.1')
genemodel.plot_domain(model=t_domains, start=36056239, bpstop=36056821, orientation='forward')

#Working on BaRT2v18chr3HG123230


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123230')

gene.transcript.model.plot(model=gene, gene_start=35414742, gene_bpstop=35415435, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123230')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35414742, gene_bpstop=35415435, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123230.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123230.1')

genemodel.plot(model=transcript, start=35414742, bpstop=35415435, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123230.1')
genemodel.plot_domain(model=t_domains, start=35414742, bpstop=35415435, orientation='forward')

#Working on BaRT2v18chr3HG123000


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123000')

gene.transcript.model.plot(model=gene, gene_start=33574179, gene_bpstop=33580981, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123000')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33574179, gene_bpstop=33580981, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123000.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.1')

genemodel.plot(model=transcript, start=33574179, bpstop=33580939, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.1')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580939, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.14498270915327988, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.10863328014625956, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.19966441245490824, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.11262308465199704, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.13261684764407394, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.1619680983144552, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.1062718832930631, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.12464076183574939, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.10


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.10')

genemodel.plot(model=transcript, start=33574179, bpstop=33580950, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.10')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580950, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.16382403427835168, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.14124993269849706, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.1430876097079456, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.31385898327241635, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.1660631827059283, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.12432409175578507, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.3885606295931986, haplotypes=c('orange'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.15153209266493173, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.11


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.11')

genemodel.plot(model=transcript, start=33574179, bpstop=33580974, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.11')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580974, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.1504340855204609, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1977108343951073, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.11827225207181757, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.38502394091326425, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.13855410744512914, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.10386267545917366, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.3315268812464967, haplotypes=c('orange'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.1537544518097815, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.12


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.12')

genemodel.plot(model=transcript, start=33574179, bpstop=33580970, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.12')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580970, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.12809758977061309, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1333673632159337, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.1945893811160902, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.3237984668946034, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.14069470413903673, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.1678513478527937, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.36295379559368524, haplotypes=c('orange'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.19879497373845018, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.13


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.13')

genemodel.plot(model=transcript, start=33574179, bpstop=33580935, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.13')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580935, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.11666545108231816, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1281825702145572, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.10595032920700986, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.3888293098894098, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.12260762774497547, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.14535315528587686, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.33224181987573465, haplotypes=c('orange'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.1250957993106027, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.14


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.14')

genemodel.plot(model=transcript, start=33574183, bpstop=33580950, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.14')
genemodel.plot_domain(model=t_domains, start=33574183, bpstop=33580950, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.18062257691259312, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.12284342893304767, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.11475974415444803, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.31941542461442973, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.1566145204783813, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.15176928388527336, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.12251458352980851, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.16379142139174924, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.15


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.15')

genemodel.plot(model=transcript, start=33574183, bpstop=33580954, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.15')
genemodel.plot_domain(model=t_domains, start=33574183, bpstop=33580954, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.12186277706155398, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.15293713077219817, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.13423081169034287, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.16875566106962953, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.1411090177605903, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.12677130791895225, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.17588585135432908, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.18487831619748316, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.16


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.16')

genemodel.plot(model=transcript, start=33574187, bpstop=33580904, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.16')
genemodel.plot_domain(model=t_domains, start=33574187, bpstop=33580904, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.1817974849571829, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.18321016487745817, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.1710094837845402, haplotypes=c('blue'))

mutation.plot(33576083, 33576083, text='', col='black', drop=-0.11997638460665222, haplotypes=c('blue'))

mutation.plot(33576106, 33576106, text='', col='black', drop=-0.17828128620103584, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.16159633443231025, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.1964520517525857, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.11368930173380438, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.1573438270296545, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.2')

genemodel.plot(model=transcript, start=33574179, bpstop=33580981, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.2')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580981, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.16270757133084363, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.170720289342144, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.14307397954332973, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.1154604076332969, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.18395947086919806, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.19008866867266444, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.10904033956268597, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.13337186149787958, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.3')

genemodel.plot(model=transcript, start=33574179, bpstop=33580918, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.3')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580918, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.10895005121157587, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1830460681331068, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.16493062328078933, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.11665751804811086, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.1184645831944956, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.1546512027043966, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.1212885012375088, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.16096164645541466, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.4')

genemodel.plot(model=transcript, start=33574179, bpstop=33580918, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.4')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580918, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.19233322238448572, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.10116180815628413, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.134861297463557, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.14213429731860389, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.11342345284809026, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.11147072583052821, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.16683740203885944, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.16062808165219794, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.5')

genemodel.plot(model=transcript, start=33574179, bpstop=33580954, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.5')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580954, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.14575287158517017, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.13754475869148555, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.17158503879579826, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.1759154668998636, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.10925946918855745, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.15449698562473207, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.17811009684350926, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.17411074734712673, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.6')

genemodel.plot(model=transcript, start=33574179, bpstop=33580954, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.6')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580954, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.11656334484506711, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1879363052487748, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.1078779366221761, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.1301134633903049, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.16096134018430802, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.11525755969474837, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.13412392702734327, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.14035709234445493, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.7')

genemodel.plot(model=transcript, start=33574179, bpstop=33580981, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.7')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580981, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.1385051872536039, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1748442461920937, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.10055002730301163, haplotypes=c('blue'))

mutation.plot(33576083, 33576083, text='', col='black', drop=-0.12101478860460824, haplotypes=c('blue'))

mutation.plot(33576106, 33576106, text='', col='black', drop=-0.15447060100097187, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.10881321542324768, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.10496290796821262, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.11241557736069491, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.10110544206899466, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.1717302891284758, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.8


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.8')

genemodel.plot(model=transcript, start=33574179, bpstop=33580935, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.8')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580935, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.10081920505869374, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.1418815956318346, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.10708276637771125, haplotypes=c('blue'))

mutation.plot(33576083, 33576083, text='', col='black', drop=-0.16626726618278223, haplotypes=c('blue'))

mutation.plot(33576106, 33576106, text='', col='black', drop=-0.1973462438704402, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.10846990055733317, haplotypes=c('blue'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.19627830031463755, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.1211312017587976, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.13732156495241926, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.1308885522500994, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123000.9


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123000.9')

genemodel.plot(model=transcript, start=33574179, bpstop=33580970, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123000.9')
genemodel.plot_domain(model=t_domains, start=33574179, bpstop=33580970, orientation='forward')

mutation.plot(33574220, 33574220, text='', col='black', drop=-0.19201827635887325, haplotypes=c('blue'))

mutation.plot(33574254, 33574254, text='', col='black', drop=-0.11265082130117955, haplotypes=c('blue'))

mutation.plot(33574321, 33574321, text='', col='black', drop=-0.14149357690183761, haplotypes=c('blue'))

mutation.plot(33577017, 33577017, text='', col='black', drop=-0.3916325699447231, haplotypes=c('orange'))

mutation.plot(33577727, 33577727, text='', col='black', drop=-0.10092829950309087, haplotypes=c('blue'))

mutation.plot(33577868, 33577868, text='', col='black', drop=-0.15140535672837355, haplotypes=c('blue'))

mutation.plot(33578683, 33578683, text='', col='black', drop=-0.14829851039021336, haplotypes=c('blue'))

mutation.plot(33580916, 33580916, text='', col='black', drop=-0.12821497880064509, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123490


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123490')

gene.transcript.model.plot(model=gene, gene_start=36856182, gene_bpstop=36859459, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123490')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36856182, gene_bpstop=36859459, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123490.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123490.1')

genemodel.plot(model=transcript, start=36856182, bpstop=36859459, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123490.1')
genemodel.plot_domain(model=t_domains, start=36856182, bpstop=36859459, orientation='forward')

mutation.plot(36857104, 36857104, text='', col='black', drop=-0.14916252221889498, haplotypes=c('blue'))

mutation.plot(36857114, 36857114, text='', col='black', drop=-0.17868200419389804, haplotypes=c('blue'))

mutation.plot(36859120, 36859120, text='', col='black', drop=-0.10384048013771335, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122890


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122890')

gene.transcript.model.plot(model=gene, gene_start=33164843, gene_bpstop=33166254, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122890')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33164843, gene_bpstop=33166254, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG122890.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122890.1')

genemodel.plot(model=transcript, start=33164843, bpstop=33166254, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122890.1')
genemodel.plot_domain(model=t_domains, start=33164843, bpstop=33166254, orientation='forward')

mutation.plot(33165496, 33165496, text='', col='black', drop=-0.14358890256469073, haplotypes=c('blue'))

mutation.plot(33165646, 33165646, text='', col='black', drop=-0.11215841779426369, haplotypes=c('blue'))

mutation.plot(33165734, 33165734, text='', col='black', drop=-0.3890116370472139, haplotypes=c('orange'))

mutation.plot(33166067, 33166067, text='', col='black', drop=-0.15222278513303172, haplotypes=c('blue'))

mutation.plot(33166082, 33166082, text='', col='black', drop=-0.16134590846270583, haplotypes=c('blue'))

mutation.plot(33166106, 33166106, text='', col='black', drop=-0.12918503715552182, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122890.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122890.2')

genemodel.plot(model=transcript, start=33165341, bpstop=33166061, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122890.2')
genemodel.plot_domain(model=t_domains, start=33165341, bpstop=33166061, orientation='forward')

mutation.plot(33165496, 33165496, text='', col='black', drop=-0.12140618062690815, haplotypes=c('blue'))

mutation.plot(33165646, 33165646, text='', col='black', drop=-0.1750345452451048, haplotypes=c('blue'))

mutation.plot(33165734, 33165734, text='', col='black', drop=-0.3456408665399776, haplotypes=c('orange'))

#Working on BaRT2v18chr3HG123450


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123450')

gene.transcript.model.plot(model=gene, gene_start=36436892, gene_bpstop=36437232, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123450')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36436892, gene_bpstop=36437232, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123450.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123450.1')

genemodel.plot(model=transcript, start=36436892, bpstop=36437232, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123450.1

#Working on BaRT2v18chr3HG123110


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123110')

gene.transcript.model.plot(model=gene, gene_start=34748337, gene_bpstop=34751773, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123110')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=34748337, gene_bpstop=34751773, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123110.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.1')

genemodel.plot(model=transcript, start=34748337, bpstop=34751695, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.1

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.13583949521170507, haplotypes=c('blue'))

mutation.plot(34748574, 34748574, text='', col='black', drop=-0.16524806362442856, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.13369709157667092, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.19225212945835543, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.11491880249419979, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.1822609039090357, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.14851756735926322, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.1676031324278658, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.1429284287677566, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.17766546843052922, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.15157464259262488, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.11522199053659227, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.11083236410584252, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.10284473048173703, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.1394411160131879, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.16799913189525434, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.12377637621399637, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.12882183456469937, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.19037281519411775, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.14186062858679913, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.1273824853889646, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.14111600749441183, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.15579493199857308, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.1334728519337252, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.19316620634940362, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.13435447921336635, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.10


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.10')

genemodel.plot(model=transcript, start=34748408, bpstop=34750957, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.10

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.1148581503571536, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.10873025889166064, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.1143199690481777, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1115454265261403, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.15719094431345346, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.18323132974359121, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.12361700454337841, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.1502837736791111, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.18719996443204, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.16261899925802148, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1792558528012791, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.12363339213309146, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.10184166845423172, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.13884408030411588, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.18366133268512558, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.10067963610545928, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.11


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.11')

genemodel.plot(model=transcript, start=34748531, bpstop=34751598, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123110.11')
genemodel.plot_domain(model=t_domains, start=34748531, bpstop=34751598, orientation='forward')

mutation.plot(34748574, 34748574, text='', col='black', drop=-0.12232716246468783, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.1323296349163281, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.10896244746230178, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.12275950093483752, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.11676154632071682, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.13371546683475188, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.14149885529826206, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.18510125013367557, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.1412117260654281, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.19660254649631057, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.10656630863836486, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.1925848664062321, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.12516666482745722, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.14278685847778852, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3781003973348599, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.395240409840085, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.13433819678952993, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5218028581325117, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.3941336235419743, haplotypes=c('orange'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.3339612324333489, haplotypes=c('orange'))

mutation.plot(34751473, 34751473, text='stop_lost', col='black', drop=-0.5733816910004887, haplotypes=c('red'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1891117763312253, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.12


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.12')

genemodel.plot(model=transcript, start=34748534, bpstop=34751706, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123110.12')
genemodel.plot_domain(model=t_domains, start=34748534, bpstop=34751706, orientation='forward')

mutation.plot(34748574, 34748574, text='', col='black', drop=-0.1676303851312833, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.14625559713781733, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.1180633636685161, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1568927152478381, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.17820253371268033, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.13243842243103185, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.1194605424268806, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.11590528073999186, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.12462641475733792, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.16892348258133838, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1707965746901574, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.10226617463945947, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.17122424562690852, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.19718339668721918, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.30080631299214944, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.327097663760971, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.10168429068185551, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5180136158854116, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.3578367929527826, haplotypes=c('orange'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.30841523851349223, haplotypes=c('orange'))

mutation.plot(34751473, 34751473, text='stop_lost', col='black', drop=-0.5510995284349465, haplotypes=c('red'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1898604797361081, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.13


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.13')

genemodel.plot(model=transcript, start=34748534, bpstop=34751768, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.13

mutation.plot(34748574, 34748574, text='', col='black', drop=-0.15966798682625855, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.17516543385994027, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.1888288169030825, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.14454717332025363, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.10681948542883477, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.11326861481855033, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.17285473058604317, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.15285151397587218, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.1980276298741574, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.16371903832663603, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.10897728071143582, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.18988567902840497, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.19388954584698082, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.10583249176655377, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3572108971349356, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.38389174453739816, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.3881006006947805, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.19600832136987173, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.15398834110707368, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.1421173658624605, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.10108968112406284, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.10502822326245968, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.14


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.14')

genemodel.plot(model=transcript, start=34748543, bpstop=34751768, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.14

mutation.plot(34748574, 34748574, text='', col='black', drop=-0.1394610516216525, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.15310632380702593, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.1532102417790716, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1952275285091014, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.15548038565099292, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.10913826333427906, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.17169459988325678, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.1893384076849057, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.1961672676555111, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.18341431012079415, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.18363278677299938, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.1222384334843508, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.11289616639683148, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.10596595618129662, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.16949167135025855, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.15086948319626908, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.17222800370427427, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.18223033549873374, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.15249574026429813, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.18771030710219908, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.11992665959234397, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.18573916008691543, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.14966844945371616, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.16187378838025246, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.18203856643607194, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.15


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.15')

genemodel.plot(model=transcript, start=34748578, bpstop=34751703, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.15

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.17205467022948312, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.10670758033856487, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.19873674516500126, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.16246958216136806, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.14938469570541377, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.12201101084330615, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.19959830415011812, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.13213138577042435, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.14643373935319762, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.15441738494658422, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.16466182897473536, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.12348172180643496, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.11381373844313408, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.16622563658391182, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.14940335610086158, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.17982821451993228, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.17900381821352807, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.186891047835345, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.16148527603262677, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.1078454355688849, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.16074155678674842, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.13240224069131082, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1263227219342596, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1922088211385308, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.16


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.16')

genemodel.plot(model=transcript, start=34748578, bpstop=34751666, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.16

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.11938975221202408, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.15305610598907798, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.12787350098000927, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.17084078743071873, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.18291406542904542, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.10916509220376298, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.1090158315148461, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.19053952447188338, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.1830938409663911, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1479956873756105, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.13886795797606222, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.17040852133354567, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.10192448178390352, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.10434746002899183, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.17642387805590298, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.121013454634537, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.15958571524600415, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.14932497817595788, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.18694681927427195, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.19820206554734374, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.15595747826670797, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.1004248634171064, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.12836456058763532, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.19551716135172115, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.17


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.17')

genemodel.plot(model=transcript, start=34748583, bpstop=34751773, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.17

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.1961690689973069, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.10911368026820799, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.17907706172571977, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.14838065285520163, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.14115868215779148, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.1464672452416008, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.15936241751655475, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.1578950254487937, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.1514572698232757, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.17280879055220094, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.18009566174265262, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.1202560948181833, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.1512700741341128, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.14507036206930354, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.11407620182837397, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.15286964268703, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1452651196041659, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1830480822719263, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.18


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.18')

genemodel.plot(model=transcript, start=34748583, bpstop=34751653, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.18

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.15020779124370776, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.19220393330216812, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1443570154065948, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.1820254806051902, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.11375095356732975, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.18778478733462642, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.10238242770188877, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.18294120245672688, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.1839573927935669, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.1361260848024841, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.17119131329069165, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.1642976951192851, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.11651096261402569, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.11317392382489871, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.18190794052420367, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.14436758609209718, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1972834055385415, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.11472186591158724, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.19


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.19')

genemodel.plot(model=transcript, start=34748583, bpstop=34751703, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.19

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.15035780785509473, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.18884251269851926, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.11541579866957119, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.1254108300709155, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.1781119300290378, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.10713118842039876, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.16394718621180776, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.13268218358116857, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.16881358576964156, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1699433013426756, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.1388373449842905, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.1476502481848333, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.18707385683725009, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.15304418966393196, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.18488592940672205, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.19694286423354151, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.19841755193452526, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.11253307753424807, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.17992999826596626, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.11247690350019504, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.13664431766961788, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.13144564298799685, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1604640867312176, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.13971920895309822, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.2')

genemodel.plot(model=transcript, start=34748340, bpstop=34751764, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.2

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.16604627520800383, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.13838461738553726, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.11564188063696154, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.17902827497754298, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.13462266932805877, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.16646133110715655, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.15531456628760296, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.12256258503169036, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.14112205679853462, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.15257752696907423, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.17310795484461092, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.1932727898656565, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.12887794681179268, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.18526716597819154, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.14014702569638574, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.1422059340006319, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.12075631942381257, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.15182045493346719, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.18136264064506885, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.10596244534872612, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.1649924386881727, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.11868284733910835, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.11842362429449838, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1777600387230946, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.16316984664221473, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.20


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.20')

genemodel.plot(model=transcript, start=34748591, bpstop=34751689, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.20

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.12350557982484353, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.16177774200910516, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.12470147983755028, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.18080981804947188, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.16471459097946145, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.17971859196086587, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.14110821810518953, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.16421394258475996, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.15366242873400007, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.17565177880502345, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.16169838719400764, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.13885882764691085, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.13838446037314225, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.39615322302609834, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.3984524174846439, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.3912943999186685, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5270200098972913, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.16822262845050645, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.12745003706592933, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.16704972163446674, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.19142904399104813, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.3')

genemodel.plot(model=transcript, start=34748346, bpstop=34751760, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.3

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.10147364807780848, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.18057639978128714, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.13032336615037654, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.11158594301201519, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.17329295591729527, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.12168871767795539, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.10131530251988968, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.17187414766402426, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.15400081690987494, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.12308267289853293, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.17772773808305084, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.17640979724019198, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.1772918906890411, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.16006754643828214, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.18180184783861275, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.11455994693880228, haplotypes=c('blue'))

mutation.plot(34750986, 34750986, text='', col='black', drop=-0.19613260114409403, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.11590706523232501, haplotypes=c('blue'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.13083118585046039, haplotypes=c('blue'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.11568351950761661, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.18108294720974627, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.1755756962374283, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.17713512139100623, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.14452893155594737, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.14154789708774718, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.4')

genemodel.plot(model=transcript, start=34748346, bpstop=34751760, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.4

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.17077559475593793, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.1731554030217825, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.16655091525177962, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1281064432298036, haplotypes=c('blue'))

mutation.plot(34749487, 34749487, text='', col='black', drop=-0.19137535484374357, haplotypes=c('blue'))

mutation.plot(34749525, 34749525, text='', col='black', drop=-0.14394043324807895, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.11114455117837826, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.1970775376848322, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.19011017050218157, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.19330601048122648, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.1679598840791714, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.15981563284175582, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.12193587114153562, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.16080349675790395, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.1929098099774139, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.15164921529173087, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3894849985301979, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.3900096036152407, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.38152776173137914, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.504902611374112, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.14007720771558127, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.11560982378090284, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.10026208216093949, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.12879409189830035, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.5')

genemodel.plot(model=transcript, start=34748375, bpstop=34751635, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123110.5')
genemodel.plot_domain(model=t_domains, start=34748375, bpstop=34751635, orientation='forward')

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.13378151606117067, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.16623227452811074, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.14381316471856198, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1904315339826546, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.13514501079321803, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.17665808511863612, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.11819151187153491, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.12871522153546844, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.17579491368295666, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.16380671960038015, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1729620925599311, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.11569125519747243, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.14910821994933782, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.19107687869561613, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3384891368994463, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.3083027865225844, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.18739140902322654, haplotypes=c('blue'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5881030036195057, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.33068131360959835, haplotypes=c('orange'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.31445157683388647, haplotypes=c('orange'))

mutation.plot(34751473, 34751473, text='stop_lost', col='black', drop=-0.5253739907449236, haplotypes=c('red'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1357147674835834, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.6')

genemodel.plot(model=transcript, start=34748376, bpstop=34751702, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.6

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.10362565364752846, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.18489573955859964, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.11418877690403832, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.13525030123063994, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.11352185139158638, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.1355550472445457, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.1944945547543573, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.12930941115885242, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.37079543032025664, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.34938435511409827, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.35873693272650836, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5575256229804753, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.19020315861075646, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.19474779868909362, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.10998961458315663, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.17810515243778796, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.7')

genemodel.plot(model=transcript, start=34748389, bpstop=34751769, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.7

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.15531123365949404, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.14664162938651384, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.180070564920982, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.13628179761203024, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.18105424984157054, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.11202285554512836, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.18337010425272626, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.13742010776402888, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3922548733064858, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.38142654922606906, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.3671033537980771, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='frameshift', col='black', drop=-0.5649642555200346, haplotypes=c('red'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.16632575063805236, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.1903625748780742, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.1833024033236465, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.1634629587572781, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.8


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.8')

genemodel.plot(model=transcript, start=34748389, bpstop=34751523, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.8

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.12570377294636, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.13692647556759996, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.13295392152148017, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1228841935551733, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.12503288939487758, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.12269708923812295, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.18608201691844034, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.12060970055840971, haplotypes=c('blue'))

mutation.plot(34751060, 34751060, text='', col='black', drop=-0.3925048192989061, haplotypes=c('orange'))

mutation.plot(34751106, 34751106, text='', col='black', drop=-0.34986350398262334, haplotypes=c('orange'))

mutation.plot(34751148, 34751148, text='', col='black', drop=-0.3685367994907396, haplotypes=c('orange'))

mutation.plot(34751252, 34751257, text='', col='black', drop=-0.123361417825989, haplotypes=c('blue'))

mutation.plot(34751342, 34751342, text='', col='black', drop=-0.12487622666066878, haplotypes=c('blue'))

mutation.plot(34751416, 34751416, text='', col='black', drop=-0.16719909233157326, haplotypes=c('blue'))

mutation.plot(34751473, 34751473, text='', col='black', drop=-0.16094319797741405, haplotypes=c('blue'))

mutation.plot(34751502, 34751502, text='', col='black', drop=-0.19584355418181848, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123110.9


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123110.9')

genemodel.plot(model=transcript, start=34748389, bpstop=34750957, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123110.9

mutation.plot(34748456, 34748456, text='', col='black', drop=-0.13783769306248878, haplotypes=c('blue'))

mutation.plot(34748648, 34748649, text='', col='black', drop=-0.15882343289132467, haplotypes=c('blue'))

mutation.plot(34749036, 34749036, text='', col='black', drop=-0.19364721654859773, haplotypes=c('blue'))

mutation.plot(34749044, 34749044, text='', col='black', drop=-0.1127971978794913, haplotypes=c('blue'))

mutation.plot(34749617, 34749617, text='', col='black', drop=-0.1593676256818312, haplotypes=c('blue'))

mutation.plot(34749683, 34749683, text='', col='black', drop=-0.15189785532949904, haplotypes=c('blue'))

mutation.plot(34749723, 34749723, text='', col='black', drop=-0.11565171228122512, haplotypes=c('blue'))

mutation.plot(34749786, 34749786, text='', col='black', drop=-0.15701655204572196, haplotypes=c('blue'))

mutation.plot(34750071, 34750071, text='', col='black', drop=-0.17270541883310925, haplotypes=c('blue'))

mutation.plot(34750157, 34750157, text='', col='black', drop=-0.14611336644001888, haplotypes=c('blue'))

mutation.plot(34750203, 34750203, text='', col='black', drop=-0.1560727349010105, haplotypes=c('blue'))

mutation.plot(34750217, 34750217, text='', col='black', drop=-0.19827020976801973, haplotypes=c('blue'))

mutation.plot(34750342, 34750342, text='', col='black', drop=-0.10839584053336307, haplotypes=c('blue'))

mutation.plot(34750354, 34750354, text='', col='black', drop=-0.1588008686925566, haplotypes=c('blue'))

mutation.plot(34750477, 34750477, text='', col='black', drop=-0.11595717559950644, haplotypes=c('blue'))

mutation.plot(34750705, 34750705, text='', col='black', drop=-0.15647782517046627, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123300


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123300')

gene.transcript.model.plot(model=gene, gene_start=36053203, gene_bpstop=36053479, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123300')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36053203, gene_bpstop=36053479, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123300.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123300.1')

genemodel.plot(model=transcript, start=36053203, bpstop=36053479, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123300.1

#Working on BaRT2v18chr3HG123020


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123020')

gene.transcript.model.plot(model=gene, gene_start=33585285, gene_bpstop=33589376, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123020')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33585285, gene_bpstop=33589376, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123020.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123020.1')

genemodel.plot(model=transcript, start=33585285, bpstop=33589376, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123020.1')
genemodel.plot_domain(model=t_domains, start=33585285, bpstop=33589376, orientation='forward')

mutation.plot(33585379, 33585379, text='', col='black', drop=-0.12035065503718498, haplotypes=c('blue'))

mutation.plot(33585425, 33585425, text='', col='black', drop=-0.17915889829824505, haplotypes=c('blue'))

mutation.plot(33585446, 33585446, text='', col='black', drop=-0.1316561413959426, haplotypes=c('blue'))

mutation.plot(33585484, 33585484, text='', col='black', drop=-0.18468137512839175, haplotypes=c('blue'))

mutation.plot(33585917, 33585917, text='', col='black', drop=-0.3664864581757753, haplotypes=c('orange'))

mutation.plot(33586222, 33586222, text='', col='black', drop=-0.18979923810410482, haplotypes=c('blue'))

mutation.plot(33586252, 33586252, text='', col='black', drop=-0.1286565098544924, haplotypes=c('blue'))

mutation.plot(33587087, 33587087, text='', col='black', drop=-0.30327847496879745, haplotypes=c('orange'))

mutation.plot(33587360, 33587360, text='', col='black', drop=-0.30850805693342187, haplotypes=c('orange'))

mutation.plot(33587459, 33587459, text='', col='black', drop=-0.36079710338164733, haplotypes=c('orange'))

mutation.plot(33588106, 33588106, text='', col='black', drop=-0.18222944401094301, haplotypes=c('blue'))

mutation.plot(33588226, 33588226, text='', col='black', drop=-0.1682504454208662, haplotypes=c('blue'))

mutation.plot(33589063, 33589063, text='', col='black', drop=-0.1980778496142659, haplotypes=c('blue'))

mutation.plot(33589215, 33589215, text='', col='black', drop=-0.15938981651759293, haplotypes=c('blue'))

mutation.plot(33589321, 33589322, text='', col='black', drop=-0.15646100187820672, haplotypes=c('blue'))

mutation.plot(33589322, 33589322, text='', col='black', drop=-0.13881210883638273, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123020.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123020.2')

genemodel.plot(model=transcript, start=33587281, bpstop=33589375, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123020.2')
genemodel.plot_domain(model=t_domains, start=33587281, bpstop=33589375, orientation='forward')

mutation.plot(33587360, 33587360, text='', col='black', drop=-0.197861287303099, haplotypes=c('blue'))

mutation.plot(33587459, 33587459, text='', col='black', drop=-0.34362439839596703, haplotypes=c('orange'))

mutation.plot(33588106, 33588106, text='', col='black', drop=-0.14905341973894592, haplotypes=c('blue'))

mutation.plot(33588226, 33588226, text='', col='black', drop=-0.10177776916904822, haplotypes=c('blue'))

mutation.plot(33589063, 33589063, text='', col='black', drop=-0.18555223270140572, haplotypes=c('blue'))

mutation.plot(33589215, 33589215, text='', col='black', drop=-0.12100232210823704, haplotypes=c('blue'))

mutation.plot(33589321, 33589322, text='', col='black', drop=-0.19329961766580236, haplotypes=c('blue'))

mutation.plot(33589322, 33589322, text='', col='black', drop=-0.12486832550949466, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122920')

gene.transcript.model.plot(model=gene, gene_start=33237820, gene_bpstop=33244461, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122920')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33237820, gene_bpstop=33244461, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG122920.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.1')

genemodel.plot(model=transcript, start=33237820, bpstop=33244358, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.1')
genemodel.plot_domain(model=t_domains, start=33237820, bpstop=33244358, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.10446575217560011, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.17806677193450024, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.3842298465954509, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.1451365877603223, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.1529000924622807, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.2')

genemodel.plot(model=transcript, start=33237821, bpstop=33244410, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.2')
genemodel.plot_domain(model=t_domains, start=33237821, bpstop=33244410, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.14372589897991514, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.14620465742997488, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.3022488421020753, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.13362908008674532, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.18664403281147385, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.3')

genemodel.plot(model=transcript, start=33237821, bpstop=33244415, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.3')
genemodel.plot_domain(model=t_domains, start=33237821, bpstop=33244415, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.10350450843032558, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.1356471483228658, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.30436118244030463, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.1676090903424091, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.13815091382395384, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.4')

genemodel.plot(model=transcript, start=33237824, bpstop=33244459, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.4')
genemodel.plot_domain(model=t_domains, start=33237824, bpstop=33244459, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.10576357612387673, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.15802776733765647, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.3989124532974264, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.164354998554315, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.16321029195299247, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.5')

genemodel.plot(model=transcript, start=33237826, bpstop=33244414, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.5')
genemodel.plot_domain(model=t_domains, start=33237826, bpstop=33244414, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.18624731930452817, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.1650435189443294, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.32241300988998406, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.18992115933100087, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.1257271389232049, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.6')

genemodel.plot(model=transcript, start=33237826, bpstop=33244461, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.6')
genemodel.plot_domain(model=t_domains, start=33237826, bpstop=33244461, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.12271161497668209, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.11039916327971663, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.3155376162494028, haplotypes=c('orange'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.1793292293708938, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.18480575350796632, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122920.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122920.7')

genemodel.plot(model=transcript, start=33237826, bpstop=33244414, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122920.7')
genemodel.plot_domain(model=t_domains, start=33237826, bpstop=33244414, orientation='forward')

mutation.plot(33237835, 33237835, text='', col='black', drop=-0.1430400991351629, haplotypes=c('blue'))

mutation.plot(33237849, 33237849, text='', col='black', drop=-0.13205940814492143, haplotypes=c('blue'))

mutation.plot(33241825, 33241825, text='', col='black', drop=-0.1366174540185566, haplotypes=c('blue'))

mutation.plot(33244171, 33244171, text='', col='black', drop=-0.18738576647722502, haplotypes=c('blue'))

mutation.plot(33244188, 33244188, text='', col='black', drop=-0.10916509743638804, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123100


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123100')

gene.transcript.model.plot(model=gene, gene_start=34646064, gene_bpstop=34648257, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123100')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=34646064, gene_bpstop=34648257, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123100.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123100.1')

genemodel.plot(model=transcript, start=34646064, bpstop=34648257, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123100.1')
genemodel.plot_domain(model=t_domains, start=34646064, bpstop=34648257, orientation='forward')

mutation.plot(34646132, 34646132, text='', col='black', drop=-0.10825234574052986, haplotypes=c('blue'))

mutation.plot(34646150, 34646151, text='', col='black', drop=-0.13704757033875212, haplotypes=c('blue'))

mutation.plot(34646226, 34646226, text='', col='black', drop=-0.1282460785421113, haplotypes=c('blue'))

mutation.plot(34646299, 34646299, text='', col='black', drop=-0.16974923534059588, haplotypes=c('blue'))

mutation.plot(34647607, 34647607, text='', col='black', drop=-0.17409716732377345, haplotypes=c('blue'))

mutation.plot(34647682, 34647682, text='', col='black', drop=-0.1231726796014646, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122880')

gene.transcript.model.plot(model=gene, gene_start=33138640, gene_bpstop=33145156, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122880')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33138640, gene_bpstop=33145156, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG122880.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.1')

genemodel.plot(model=transcript, start=33138640, bpstop=33145156, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.1')
genemodel.plot_domain(model=t_domains, start=33138640, bpstop=33145156, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.15465190953503735, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.33886960232470364, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3991122353713158, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.30068729944581, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.12524773655807558, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.13055536844226676, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.18016631942938097, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.1579605955058737, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.16659788172894696, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.15899053943706387, haplotypes=c('blue'))

mutation.plot(33145057, 33145057, text='', col='black', drop=-0.17162730802699272, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.10


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.10')

genemodel.plot(model=transcript, start=33139613, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.10')
genemodel.plot_domain(model=t_domains, start=33139613, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3124721457079805, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.30122899143934967, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3218014368417165, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.14713916482175224, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.12778790911312293, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.11476689070988513, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.11476097987234618, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.13464628196291112, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.19130694167378717, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.11


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.11')

genemodel.plot(model=transcript, start=33139616, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.11')
genemodel.plot_domain(model=t_domains, start=33139616, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.33567304938943604, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.31608673378386243, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3974327324210037, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.13755722029549017, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.17723150704940316, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.13434743484025918, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.10729449191606925, haplotypes=c('blue'))

mutation.plot(33143871, 33143871, text='', col='black', drop=-0.390925036515918, haplotypes=c('orange'))

mutation.plot(33143999, 33143999, text='', col='black', drop=-0.17652797324733632, haplotypes=c('blue'))

mutation.plot(33144113, 33144114, text='', col='black', drop=-0.12654851628809863, haplotypes=c('blue'))

mutation.plot(33144117, 33144117, text='', col='black', drop=-0.1982446375044009, haplotypes=c('blue'))

mutation.plot(33144257, 33144257, text='', col='black', drop=-0.19016476363902945, haplotypes=c('blue'))

mutation.plot(33144282, 33144282, text='', col='black', drop=-0.1630305916447224, haplotypes=c('blue'))

mutation.plot(33144324, 33144324, text='', col='black', drop=-0.11047386409025908, haplotypes=c('blue'))

mutation.plot(33144371, 33144371, text='', col='black', drop=-0.16974132983613643, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.16690721706751868, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.19021399459286736, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.12


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.12')

genemodel.plot(model=transcript, start=33139641, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.12')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3356783919962077, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.35476724939109117, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.30613492569592676, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1900660059772109, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.13176117050587915, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.1650991426448357, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.11859596039378255, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.10710655520358714, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.14380439797115938, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.13


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.13')

genemodel.plot(model=transcript, start=33139641, bpstop=33145012, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.13')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145012, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3738242374315611, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.36493585483088237, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3378334594611231, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.13397545760543897, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.113413040177347, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.19523433574535323, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.11406103695807542, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.12111842104339172, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.19304949016866246, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.14


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.14')

genemodel.plot(model=transcript, start=33139641, bpstop=33145025, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.14')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145025, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3548921127974503, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3970880316558719, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.32231767677927753, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.18445936785784484, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.12797828644696063, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.14030354406154716, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.146777966226888, haplotypes=c('blue'))

mutation.plot(33143871, 33143871, text='', col='black', drop=-0.30293810743902877, haplotypes=c('orange'))

mutation.plot(33143999, 33143999, text='', col='black', drop=-0.11099579568640229, haplotypes=c('blue'))

mutation.plot(33144113, 33144114, text='', col='black', drop=-0.11997171356263074, haplotypes=c('blue'))

mutation.plot(33144117, 33144117, text='', col='black', drop=-0.1447079105345986, haplotypes=c('blue'))

mutation.plot(33144257, 33144257, text='', col='black', drop=-0.14776255772207295, haplotypes=c('blue'))

mutation.plot(33144282, 33144282, text='', col='black', drop=-0.14959779398964893, haplotypes=c('blue'))

mutation.plot(33144324, 33144324, text='', col='black', drop=-0.17435111899497935, haplotypes=c('blue'))

mutation.plot(33144371, 33144371, text='', col='black', drop=-0.17690467720851447, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.16047439377033185, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.16806605241809658, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.15


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.15')

genemodel.plot(model=transcript, start=33139641, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.15')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3541887406521815, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3715396690215549, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.34608443287830337, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.19632849224307625, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.1967954745606974, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.19233992287423338, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.13606297509514934, haplotypes=c('blue'))

mutation.plot(33142163, 33142163, text='', col='black', drop=-0.13288263871381759, haplotypes=c('blue'))

mutation.plot(33142282, 33142282, text='', col='black', drop=-0.1863123676741895, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.1530152266082018, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.17747951982182558, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.16


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.16')

genemodel.plot(model=transcript, start=33139641, bpstop=33145016, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.16')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145016, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3527489138114898, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.39635040778330055, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.39523956316974573, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.17983814611557591, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.10638349952210188, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.16292876625066674, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.1527804174881097, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.18635787245959354, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.1966082921690507, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.17


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.17')

genemodel.plot(model=transcript, start=33139641, bpstop=33145016, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.17')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145016, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.39330886847261354, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.30955962139947124, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.36919686948942926, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.11152151635193665, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.10143931963383576, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.1628601156923713, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.1285022124532948, haplotypes=c('blue'))

mutation.plot(33142163, 33142163, text='', col='black', drop=-0.1629887404020144, haplotypes=c('blue'))

mutation.plot(33142282, 33142282, text='', col='black', drop=-0.18786293158956457, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.19998952638323467, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.13725128824478183, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.18


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.18')

genemodel.plot(model=transcript, start=33139641, bpstop=33145012, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.18')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145012, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3581232645074481, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.33384678204558493, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3438450873624879, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.17325586242922753, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.1912418402243352, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.18449541871937317, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.122619998187953, haplotypes=c('blue'))

mutation.plot(33143871, 33143871, text='', col='black', drop=-0.1712675467964469, haplotypes=c('blue'))

mutation.plot(33143999, 33143999, text='', col='black', drop=-0.11846760278650978, haplotypes=c('blue'))

mutation.plot(33144113, 33144114, text='', col='black', drop=-0.18797601578430578, haplotypes=c('blue'))

mutation.plot(33144117, 33144117, text='', col='black', drop=-0.12459821945401028, haplotypes=c('blue'))

mutation.plot(33144257, 33144257, text='', col='black', drop=-0.15788065076804997, haplotypes=c('blue'))

mutation.plot(33144282, 33144282, text='', col='black', drop=-0.10144441534170094, haplotypes=c('blue'))

mutation.plot(33144324, 33144324, text='', col='black', drop=-0.1140603150078549, haplotypes=c('blue'))

mutation.plot(33144371, 33144371, text='', col='black', drop=-0.17229192524679887, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.10641655678942669, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.15559729341944473, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.19


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.19')

genemodel.plot(model=transcript, start=33139641, bpstop=33145012, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.19')
genemodel.plot_domain(model=t_domains, start=33139641, bpstop=33145012, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3376620651276356, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.30392479295621927, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.33378269200600524, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.18351586847789095, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.13546978324205716, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.18771535072053855, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.14587239252498416, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.19487334262594425, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.18802065731111361, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.2')

genemodel.plot(model=transcript, start=33138701, bpstop=33145085, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.2')
genemodel.plot_domain(model=t_domains, start=33138701, bpstop=33145085, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.16171591559773316, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.30376621733199716, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3261239133698453, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.37878936562743193, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.13783426941436533, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.1770491453662563, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.1402689472453948, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.10106632589380882, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.10550970074429325, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.11775773930586303, haplotypes=c('blue'))

mutation.plot(33145057, 33145057, text='', col='black', drop=-0.141999748246199, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.20


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.20')

genemodel.plot(model=transcript, start=33139680, bpstop=33145012, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.20')
genemodel.plot_domain(model=t_domains, start=33139680, bpstop=33145012, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3061027181410383, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.39097502824571867, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.38280536161325907, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1600498661015283, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.14220197883210875, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.35933874023481444, haplotypes=c('orange'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.13192168605612456, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.11518037138549361, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.17283174285638558, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.21


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.21')

genemodel.plot(model=transcript, start=33139680, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.21')
genemodel.plot_domain(model=t_domains, start=33139680, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3912329182456221, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3976553137169333, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.31956212141362944, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.17750637999082736, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.18504138354489869, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.357114088831894, haplotypes=c('orange'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.17403168887549358, haplotypes=c('blue'))

mutation.plot(33142163, 33142163, text='', col='black', drop=-0.19289652688668207, haplotypes=c('blue'))

mutation.plot(33142282, 33142282, text='', col='black', drop=-0.10975285175510455, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.14084628160722465, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.10330385166738504, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.22


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.22')

genemodel.plot(model=transcript, start=33139680, bpstop=33145005, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.22')
genemodel.plot_domain(model=t_domains, start=33139680, bpstop=33145005, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.30845676825907326, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3080768338796302, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3878500740229943, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1718099624286924, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.1003375656027367, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.13301799375678153, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.13778477211185153, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.18129386580192863, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.10485048130918961, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.3')

genemodel.plot(model=transcript, start=33139186, bpstop=33144902, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.3')
genemodel.plot_domain(model=t_domains, start=33139186, bpstop=33144902, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.145215561179396, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.36696736835481414, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3384225220546378, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.37712459154717337, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.12541379514310763, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.13971293353909145, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.1584692326515095, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.1622074443323264, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.1956272988582257, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.1266314445361527, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.4')

genemodel.plot(model=transcript, start=33139187, bpstop=33145131, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.4')
genemodel.plot_domain(model=t_domains, start=33139187, bpstop=33145131, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.1006717576870476, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3033696662381948, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.33423559523870744, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.38102059463954546, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.12199696761508977, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.13694420158626286, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.3436948076664975, haplotypes=c('orange'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.19047890869048587, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.15646244480504495, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.18862344900192915, haplotypes=c('blue'))

mutation.plot(33145057, 33145057, text='', col='black', drop=-0.12243556401204897, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.5')

genemodel.plot(model=transcript, start=33139191, bpstop=33144887, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.5')
genemodel.plot_domain(model=t_domains, start=33139191, bpstop=33144887, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.12159619617615004, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.38958235634876315, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3462545655140755, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3129795886844234, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1032612638878894, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.12868027202474153, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.18617926044769417, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.170469471674359, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.15769095318580661, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.16061019073354754, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.6')

genemodel.plot(model=transcript, start=33139225, bpstop=33144853, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.6')
genemodel.plot_domain(model=t_domains, start=33139225, bpstop=33144853, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.13850686470278653, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3171451014283104, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.36286873218863674, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3874326682731535, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1277307515266849, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.16666086111938244, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.16697942333929128, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.16856279370802807, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.12164764448949184, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.1163692344107271, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.7')

genemodel.plot(model=transcript, start=33139560, bpstop=33145007, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.7')
genemodel.plot_domain(model=t_domains, start=33139560, bpstop=33145007, orientation='forward')

mutation.plot(33139539, 33139545, text='', col='black', drop=-0.161372285593248, haplotypes=c('blue'))

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3473365985835173, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3593981136905157, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3092859696638616, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.13429591445433153, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.1878876759370006, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.1835302691532405, haplotypes=c('blue'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.15426155726863897, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.16465344728084313, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.17620011666755914, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.8


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.8')

genemodel.plot(model=transcript, start=33139607, bpstop=33145016, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.8')
genemodel.plot_domain(model=t_domains, start=33139607, bpstop=33145016, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.3834812496305798, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3233211849690948, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.37102386149169847, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.1601904601358733, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.10304084269718072, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.11537490256996322, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.1402412792144236, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.15622157173825893, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122880.9


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122880.9')

genemodel.plot(model=transcript, start=33139607, bpstop=33145012, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122880.9')
genemodel.plot_domain(model=t_domains, start=33139607, bpstop=33145012, orientation='forward')

mutation.plot(33139937, 33139937, text='', col='black', drop=-0.30029090921059165, haplotypes=c('orange'))

mutation.plot(33139958, 33139958, text='', col='black', drop=-0.3440085636633914, haplotypes=c('orange'))

mutation.plot(33140588, 33140588, text='', col='black', drop=-0.3829358454577326, haplotypes=c('orange'))

mutation.plot(33140590, 33140590, text='', col='black', drop=-0.19548653680206635, haplotypes=c('blue'))

mutation.plot(33140868, 33140868, text='', col='black', drop=-0.10377844341410603, haplotypes=c('blue'))

mutation.plot(33141601, 33141601, text='', col='black', drop=-0.3118679490074867, haplotypes=c('orange'))

mutation.plot(33141652, 33141652, text='', col='black', drop=-0.11320962608920507, haplotypes=c('blue'))

mutation.plot(33144645, 33144645, text='', col='black', drop=-0.16862706081686157, haplotypes=c('blue'))

mutation.plot(33144803, 33144803, text='', col='black', drop=-0.14014429133634712, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123320


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123320')

gene.transcript.model.plot(model=gene, gene_start=36060037, gene_bpstop=36060240, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123320')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36060037, gene_bpstop=36060240, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123320.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123320.1')

genemodel.plot(model=transcript, start=36060037, bpstop=36060240, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123320.1

#Working on BaRT2v18chr3HG122930


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122930')

gene.transcript.model.plot(model=gene, gene_start=33248983, gene_bpstop=33250634, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122930')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33248983, gene_bpstop=33250634, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG122930.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122930.1')

genemodel.plot(model=transcript, start=33248983, bpstop=33250634, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122930.1')
genemodel.plot_domain(model=t_domains, start=33248983, bpstop=33250634, orientation='forward')

mutation.plot(33249270, 33249270, text='', col='black', drop=-0.18290030142120137, haplotypes=c('blue'))

mutation.plot(33249283, 33249283, text='', col='black', drop=-0.12106374938154012, haplotypes=c('blue'))

mutation.plot(33250161, 33250161, text='', col='black', drop=-0.17860309201044378, haplotypes=c('blue'))

mutation.plot(33250186, 33250186, text='', col='black', drop=-0.10312947976647245, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123330


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123330')

gene.transcript.model.plot(model=gene, gene_start=36065640, gene_bpstop=36066442, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123330')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36065640, gene_bpstop=36066442, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123330.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123330.1')

genemodel.plot(model=transcript, start=36065640, bpstop=36066442, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123330.1')
genemodel.plot_domain(model=t_domains, start=36065640, bpstop=36066442, orientation='forward')

#Working on BaRT2v18chr3HG123280


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123280')

gene.transcript.model.plot(model=gene, gene_start=36013370, gene_bpstop=36020675, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123280')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36013370, gene_bpstop=36020675, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123280.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123280.1')

genemodel.plot(model=transcript, start=36013370, bpstop=36020675, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123280.1')
genemodel.plot_domain(model=t_domains, start=36013370, bpstop=36020675, orientation='forward')

mutation.plot(36013467, 36013467, text='', col='black', drop=-0.17799027730527978, haplotypes=c('blue'))

mutation.plot(36013908, 36013908, text='', col='black', drop=-0.15127711156785478, haplotypes=c('blue'))

mutation.plot(36014042, 36014042, text='', col='black', drop=-0.36665965393960287, haplotypes=c('orange'))

mutation.plot(36019106, 36019106, text='', col='black', drop=-0.10475930927441507, haplotypes=c('blue'))

mutation.plot(36019155, 36019155, text='', col='black', drop=-0.32583489408958594, haplotypes=c('orange'))

mutation.plot(36019356, 36019356, text='', col='black', drop=-0.14019170557997607, haplotypes=c('blue'))

mutation.plot(36019370, 36019370, text='', col='black', drop=-0.31040802121212735, haplotypes=c('orange'))

mutation.plot(36019463, 36019466, text='', col='black', drop=-0.39615863954318414, haplotypes=c('orange'))

mutation.plot(36019506, 36019506, text='', col='black', drop=-0.18540030297889518, haplotypes=c('blue'))

mutation.plot(36019846, 36019846, text='', col='black', drop=-0.3527784415381251, haplotypes=c('orange'))

mutation.plot(36020006, 36020006, text='', col='black', drop=-0.3931117983426244, haplotypes=c('orange'))

mutation.plot(36020072, 36020072, text='', col='black', drop=-0.32530200964655676, haplotypes=c('orange'))

mutation.plot(36020089, 36020089, text='', col='black', drop=-0.31345164948708976, haplotypes=c('orange'))

mutation.plot(36020097, 36020097, text='', col='black', drop=-0.18036050804014783, haplotypes=c('blue'))

mutation.plot(36020510, 36020510, text='', col='black', drop=-0.1881106045076098, haplotypes=c('blue'))

mutation.plot(36020537, 36020537, text='', col='black', drop=-0.1447341503213047, haplotypes=c('blue'))

mutation.plot(36020550, 36020550, text='', col='black', drop=-0.1126135625179267, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123060


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123060')

gene.transcript.model.plot(model=gene, gene_start=33916015, gene_bpstop=33919967, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123060')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33916015, gene_bpstop=33919967, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123060.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123060.1')

genemodel.plot(model=transcript, start=33916015, bpstop=33919967, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123060.1')
genemodel.plot_domain(model=t_domains, start=33916015, bpstop=33919967, orientation='forward')

mutation.plot(33916152, 33916152, text='', col='black', drop=-0.32556035617224716, haplotypes=c('orange'))

mutation.plot(33918364, 33918364, text='', col='black', drop=-0.3797530752663033, haplotypes=c('orange'))

mutation.plot(33918809, 33918809, text='', col='black', drop=-0.11218343489332061, haplotypes=c('blue'))

mutation.plot(33918899, 33918899, text='', col='black', drop=-0.33245799414396393, haplotypes=c('orange'))

mutation.plot(33919449, 33919449, text='', col='black', drop=-0.11419299582730492, haplotypes=c('blue'))

mutation.plot(33919700, 33919700, text='', col='black', drop=-0.18431788886247824, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123060.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123060.2')

genemodel.plot(model=transcript, start=33916015, bpstop=33919967, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123060.2')
genemodel.plot_domain(model=t_domains, start=33916015, bpstop=33919967, orientation='forward')

mutation.plot(33916152, 33916152, text='', col='black', drop=-0.30073262444913956, haplotypes=c('orange'))

mutation.plot(33916709, 33916709, text='', col='black', drop=-0.11941676897409925, haplotypes=c('blue'))

mutation.plot(33916754, 33916754, text='', col='black', drop=-0.12601291115329982, haplotypes=c('blue'))

mutation.plot(33917015, 33917015, text='', col='black', drop=-0.16769762729437537, haplotypes=c('blue'))

mutation.plot(33918364, 33918364, text='', col='black', drop=-0.1264330880900023, haplotypes=c('blue'))

mutation.plot(33918809, 33918809, text='', col='black', drop=-0.17986937971155695, haplotypes=c('blue'))

mutation.plot(33918899, 33918899, text='', col='black', drop=-0.10037603702439844, haplotypes=c('blue'))

mutation.plot(33919449, 33919449, text='', col='black', drop=-0.18368701231117263, haplotypes=c('blue'))

mutation.plot(33919700, 33919700, text='', col='black', drop=-0.16181865285802735, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123120


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123120')

gene.transcript.model.plot(model=gene, gene_start=35079777, gene_bpstop=35082881, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123120')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35079777, gene_bpstop=35082881, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123120.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123120.1')

genemodel.plot(model=transcript, start=35079777, bpstop=35082881, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123120.1')
genemodel.plot_domain(model=t_domains, start=35079777, bpstop=35082881, orientation='forward')

mutation.plot(35079862, 35079862, text='', col='black', drop=-0.12413643363039904, haplotypes=c('blue'))

mutation.plot(35079920, 35079925, text='', col='black', drop=-0.19587421300441715, haplotypes=c('blue'))

mutation.plot(35079950, 35079950, text='', col='black', drop=-0.17423267738166845, haplotypes=c('blue'))

mutation.plot(35080651, 35080651, text='', col='black', drop=-0.15336085357572776, haplotypes=c('blue'))

mutation.plot(35082257, 35082257, text='', col='black', drop=-0.11366939178875114, haplotypes=c('blue'))

mutation.plot(35082460, 35082460, text='', col='black', drop=-0.13213184849036558, haplotypes=c('blue'))

mutation.plot(35082465, 35082465, text='', col='black', drop=-0.17069307151445873, haplotypes=c('blue'))

mutation.plot(35082674, 35082674, text='', col='black', drop=-0.10507563509438349, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123120.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123120.2')

genemodel.plot(model=transcript, start=35079795, bpstop=35082871, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123120.2')
genemodel.plot_domain(model=t_domains, start=35079795, bpstop=35082871, orientation='forward')

mutation.plot(35079862, 35079862, text='', col='black', drop=-0.19040686344134303, haplotypes=c('blue'))

mutation.plot(35079920, 35079925, text='', col='black', drop=-0.10284317295479434, haplotypes=c('blue'))

mutation.plot(35079950, 35079950, text='', col='black', drop=-0.18722752309692234, haplotypes=c('blue'))

mutation.plot(35080651, 35080651, text='', col='black', drop=-0.14416837359043932, haplotypes=c('blue'))

mutation.plot(35082257, 35082257, text='', col='black', drop=-0.1695772283045967, haplotypes=c('blue'))

mutation.plot(35082460, 35082460, text='', col='black', drop=-0.3365280352734655, haplotypes=c('orange'))

mutation.plot(35082465, 35082465, text='', col='black', drop=-0.3451035889428363, haplotypes=c('orange'))

mutation.plot(35082674, 35082674, text='', col='black', drop=-0.1857254275836049, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123120.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123120.3')

genemodel.plot(model=transcript, start=35079795, bpstop=35082835, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123120.3')
genemodel.plot_domain(model=t_domains, start=35079795, bpstop=35082835, orientation='forward')

mutation.plot(35079862, 35079862, text='', col='black', drop=-0.18214822547851167, haplotypes=c('blue'))

mutation.plot(35079920, 35079925, text='', col='black', drop=-0.13749075925672052, haplotypes=c('blue'))

mutation.plot(35079950, 35079950, text='', col='black', drop=-0.17538640577684048, haplotypes=c('blue'))

mutation.plot(35080651, 35080651, text='', col='black', drop=-0.10610292546887523, haplotypes=c('blue'))

mutation.plot(35082257, 35082257, text='', col='black', drop=-0.147554018916879, haplotypes=c('blue'))

mutation.plot(35082460, 35082460, text='', col='black', drop=-0.18468284676737498, haplotypes=c('blue'))

mutation.plot(35082465, 35082465, text='', col='black', drop=-0.13355764801784326, haplotypes=c('blue'))

mutation.plot(35082674, 35082674, text='', col='black', drop=-0.19738743113077667, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123120.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123120.4')

genemodel.plot(model=transcript, start=35079795, bpstop=35082835, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123120.4')
genemodel.plot_domain(model=t_domains, start=35079795, bpstop=35082835, orientation='forward')

mutation.plot(35079862, 35079862, text='', col='black', drop=-0.10223052337183955, haplotypes=c('blue'))

mutation.plot(35079920, 35079925, text='', col='black', drop=-0.1122532193097506, haplotypes=c('blue'))

mutation.plot(35079950, 35079950, text='', col='black', drop=-0.12898038220501984, haplotypes=c('blue'))

mutation.plot(35080651, 35080651, text='', col='black', drop=-0.1318747992082253, haplotypes=c('blue'))

mutation.plot(35082257, 35082257, text='', col='black', drop=-0.11263816484639332, haplotypes=c('blue'))

mutation.plot(35082460, 35082460, text='', col='black', drop=-0.16901912486677112, haplotypes=c('blue'))

mutation.plot(35082465, 35082465, text='', col='black', drop=-0.19455282337935811, haplotypes=c('blue'))

mutation.plot(35082674, 35082674, text='', col='black', drop=-0.1626815087882437, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123290


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123290')

gene.transcript.model.plot(model=gene, gene_start=36048615, gene_bpstop=36050173, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123290')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36048615, gene_bpstop=36050173, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123290.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123290.1')

genemodel.plot(model=transcript, start=36048615, bpstop=36050173, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123290.1')
genemodel.plot_domain(model=t_domains, start=36048615, bpstop=36050173, orientation='forward')

#Working on BaRT2v18chr3HG123160


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123160')

gene.transcript.model.plot(model=gene, gene_start=35151109, gene_bpstop=35152295, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123160')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35151109, gene_bpstop=35152295, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123160.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123160.1')

genemodel.plot(model=transcript, start=35151109, bpstop=35152295, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123160.1

mutation.plot(35151670, 35151670, text='', col='black', drop=-0.32649017548481607, haplotypes=c('orange'))

mutation.plot(35152022, 35152022, text='', col='black', drop=-0.10294811733112581, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123090')

gene.transcript.model.plot(model=gene, gene_start=34535739, gene_bpstop=34566620, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123090')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=34535739, gene_bpstop=34566620, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123090.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.1')

genemodel.plot(model=transcript, start=34535739, bpstop=34566620, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.1')
genemodel.plot_domain(model=t_domains, start=34535739, bpstop=34566620, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.3118447285911998, haplotypes=c('orange'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.11505461573741144, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.12378986689389077, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.16933882294508756, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.16230289383952134, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.1013510175448225, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.31549579325206323, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.1786727409582316, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.10


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.10')

genemodel.plot(model=transcript, start=34538444, bpstop=34566298, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.10')
genemodel.plot_domain(model=t_domains, start=34538444, bpstop=34566298, orientation='forward')

mutation.plot(34538680, 34538680, text='', col='black', drop=-0.18674235464381933, haplotypes=c('blue'))

mutation.plot(34538707, 34538707, text='', col='black', drop=-0.15618930699020844, haplotypes=c('blue'))

mutation.plot(34538771, 34538771, text='', col='black', drop=-0.13712483432636627, haplotypes=c('blue'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.10028164646953802, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.10886089379842045, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.16243162917119555, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.19874572054308165, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.1558638973410369, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.3718694695162416, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.1356647328217695, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.11


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.11')

genemodel.plot(model=transcript, start=34542454, bpstop=34566457, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.11')
genemodel.plot_domain(model=t_domains, start=34542454, bpstop=34566457, orientation='forward')

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.1306159650872372, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.1172311735741362, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.19829012171364604, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.3162433386818314, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.11743890160401275, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.12


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.12')

genemodel.plot(model=transcript, start=34542979, bpstop=34566355, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.12')
genemodel.plot_domain(model=t_domains, start=34542979, bpstop=34566355, orientation='forward')

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.17291887580669163, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.12696186665612988, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.10175115373285536, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.36633680357799614, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.15097116049285714, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.13


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.13')

genemodel.plot(model=transcript, start=34543350, bpstop=34566357, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.13')
genemodel.plot_domain(model=t_domains, start=34543350, bpstop=34566357, orientation='forward')

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.10887513537959743, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.1546514122076436, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.1137339755848783, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.3159646132780655, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.18415637254789222, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.2')

genemodel.plot(model=transcript, start=34535739, bpstop=34566620, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.2')
genemodel.plot_domain(model=t_domains, start=34535739, bpstop=34566620, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.36025570019031766, haplotypes=c('orange'))

mutation.plot(34538680, 34538680, text='', col='black', drop=-0.12339733042884227, haplotypes=c('blue'))

mutation.plot(34538707, 34538707, text='', col='black', drop=-0.1895795007786859, haplotypes=c('blue'))

mutation.plot(34538771, 34538771, text='', col='black', drop=-0.18940762156131521, haplotypes=c('blue'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.11560095248146984, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.11048878370343151, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.1831784331021708, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.11708262534404366, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.1683389128553568, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.12071117637760896, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.18697863058760636, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.3')

genemodel.plot(model=transcript, start=34535775, bpstop=34566310, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.3')
genemodel.plot_domain(model=t_domains, start=34535775, bpstop=34566310, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.11627968976772778, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.15730645634900053, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.15301727683438676, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.1838857726807756, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.1336947260436428, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.15563983967271777, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.18915310995560816, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.4')

genemodel.plot(model=transcript, start=34535791, bpstop=34566298, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.4')
genemodel.plot_domain(model=t_domains, start=34535791, bpstop=34566298, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.3843530315883748, haplotypes=c('orange'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.10115478716922133, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.11942086718273531, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.13924343353019, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.1616557554698318, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.19186560065861416, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.15401346256971782, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.128583330552901, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.5')

genemodel.plot(model=transcript, start=34535791, bpstop=34566298, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.5')
genemodel.plot_domain(model=t_domains, start=34535791, bpstop=34566298, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.38472958464164436, haplotypes=c('orange'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.13582599592201755, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.12118433044857434, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.13144992612334555, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.19324664524646737, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.19381006062798686, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.1546643731018157, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.10080621574573828, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.6')

genemodel.plot(model=transcript, start=34535792, bpstop=34566408, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.6')
genemodel.plot_domain(model=t_domains, start=34535792, bpstop=34566408, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.3603673882110628, haplotypes=c('orange'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.1785621878780342, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.15416787867326354, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.11806013051391298, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.10665726042469292, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.34759910578443665, haplotypes=c('orange'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.1662238106342112, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.7')

genemodel.plot(model=transcript, start=34535797, bpstop=34566298, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.7')
genemodel.plot_domain(model=t_domains, start=34535797, bpstop=34566298, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.3864292774978414, haplotypes=c('orange'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.10745799419027668, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.16728627160391385, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.18949941762645403, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.12428669766959836, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.17490333765485028, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.1419778448193174, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.19923054970535087, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.8


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.8')

genemodel.plot(model=transcript, start=34535797, bpstop=34563555, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123090.8')
genemodel.plot_domain(model=t_domains, start=34535797, bpstop=34563555, orientation='forward')

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.31340815850122106, haplotypes=c('orange'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.19444917847844018, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.10896467419197327, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.19597538835021538, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.19553235817785783, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.17629861408670136, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123090.9


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123090.9')

genemodel.plot(model=transcript, start=34536903, bpstop=34566408, orientation='forward', xaxis=T)

#No domains for BaRT2v18chr3HG123090.9

mutation.plot(34538301, 34538301, text='', col='black', drop=-0.1338986952543138, haplotypes=c('blue'))

mutation.plot(34542090, 34542090, text='', col='black', drop=-0.18071334608792167, haplotypes=c('blue'))

mutation.plot(34542313, 34542313, text='', col='black', drop=-0.10329945666831003, haplotypes=c('blue'))

mutation.plot(34560902, 34560902, text='', col='black', drop=-0.1747470307007703, haplotypes=c('blue'))

mutation.plot(34562922, 34562922, text='', col='black', drop=-0.15241505817929693, haplotypes=c('blue'))

mutation.plot(34562964, 34562964, text='', col='black', drop=-0.16127575225288882, haplotypes=c('blue'))

mutation.plot(34564241, 34564241, text='', col='black', drop=-0.19200434196354385, haplotypes=c('blue'))

mutation.plot(34564453, 34564453, text='', col='black', drop=-0.1356556159544936, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123360


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123360')

gene.transcript.model.plot(model=gene, gene_start=36113464, gene_bpstop=36114428, orientation='forward', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123360')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36113464, gene_bpstop=36114428, orientation='forward', gap=0.2)

#Working on BaRT2v18chr3HG123360.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123360.1')

genemodel.plot(model=transcript, start=36113464, bpstop=36114428, orientation='forward', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123360.1')
genemodel.plot_domain(model=t_domains, start=36113464, bpstop=36114428, orientation='forward')

#Working on BaRT2v18chr3HG123180


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123180')

gene.transcript.model.plot(model=gene, gene_start=35319467, gene_bpstop=35320226, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123180')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35319467, gene_bpstop=35320226, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123180.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123180.1')

genemodel.plot(model=transcript, start=35319467, bpstop=35320226, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123180.1

mutation.plot(35319569, 35319569, text='', col='black', drop=-0.19445552907408867, haplotypes=c('blue'))

mutation.plot(35319772, 35319772, text='stop_lost', col='black', drop=-0.5511651393211574, haplotypes=c('red'))

#Working on BaRT2v18chr3HG123340


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123340')

gene.transcript.model.plot(model=gene, gene_start=36068732, gene_bpstop=36069008, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123340')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36068732, gene_bpstop=36069008, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123340.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123340.1')

genemodel.plot(model=transcript, start=36068732, bpstop=36069008, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123340.1

#Working on BaRT2v18chr3HG123470


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123470')

gene.transcript.model.plot(model=gene, gene_start=36572423, gene_bpstop=36573695, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123470')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36572423, gene_bpstop=36573695, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123470.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123470.1')

genemodel.plot(model=transcript, start=36572423, bpstop=36573695, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123470.1

#Working on BaRT2v18chr3HG123030


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123030')

gene.transcript.model.plot(model=gene, gene_start=33648690, gene_bpstop=33652900, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123030')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33648690, gene_bpstop=33652900, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123030.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123030.1')

genemodel.plot(model=transcript, start=33648690, bpstop=33652571, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123030.1')
genemodel.plot_domain(model=t_domains, start=33648690, bpstop=33652571, orientation='reverse')

mutation.plot(33648878, 33648878, text='', col='black', drop=-0.1453930571534338, haplotypes=c('blue'))

mutation.plot(33648926, 33648926, text='', col='black', drop=-0.10104963048498845, haplotypes=c('blue'))

mutation.plot(33649037, 33649037, text='', col='black', drop=-0.11028599757624688, haplotypes=c('blue'))

mutation.plot(33649080, 33649080, text='', col='black', drop=-0.10179075831374233, haplotypes=c('blue'))

mutation.plot(33649146, 33649146, text='', col='black', drop=-0.10376292612464115, haplotypes=c('blue'))

mutation.plot(33649575, 33649575, text='', col='black', drop=-0.18549327684379985, haplotypes=c('blue'))

mutation.plot(33650866, 33650866, text='', col='black', drop=-0.10292066758319553, haplotypes=c('blue'))

mutation.plot(33652198, 33652198, text='', col='black', drop=-0.10199461268537074, haplotypes=c('blue'))

mutation.plot(33652370, 33652373, text='', col='black', drop=-0.18157887209592305, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123030.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123030.2')

genemodel.plot(model=transcript, start=33648690, bpstop=33652900, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123030.2')
genemodel.plot_domain(model=t_domains, start=33648690, bpstop=33652900, orientation='reverse')

mutation.plot(33648878, 33648878, text='', col='black', drop=-0.17865779947165333, haplotypes=c('blue'))

mutation.plot(33648926, 33648926, text='', col='black', drop=-0.15641316177696776, haplotypes=c('blue'))

mutation.plot(33649037, 33649037, text='', col='black', drop=-0.16683974371271978, haplotypes=c('blue'))

mutation.plot(33649080, 33649080, text='', col='black', drop=-0.12217530638551123, haplotypes=c('blue'))

mutation.plot(33649146, 33649146, text='', col='black', drop=-0.1646626269077483, haplotypes=c('blue'))

mutation.plot(33649575, 33649575, text='', col='black', drop=-0.13227622325229882, haplotypes=c('blue'))

mutation.plot(33650866, 33650866, text='', col='black', drop=-0.11409447676625334, haplotypes=c('blue'))

mutation.plot(33652198, 33652198, text='', col='black', drop=-0.179468567136912, haplotypes=c('blue'))

mutation.plot(33652370, 33652373, text='', col='black', drop=-0.15808125913920126, haplotypes=c('blue'))

mutation.plot(33652580, 33652581, text='', col='black', drop=-0.1881212769511871, haplotypes=c('blue'))

mutation.plot(33652830, 33652830, text='', col='black', drop=-0.16258601420043498, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123030.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123030.3')

genemodel.plot(model=transcript, start=33648691, bpstop=33652571, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123030.3')
genemodel.plot_domain(model=t_domains, start=33648691, bpstop=33652571, orientation='reverse')

mutation.plot(33648878, 33648878, text='', col='black', drop=-0.143267668818852, haplotypes=c('blue'))

mutation.plot(33648926, 33648926, text='', col='black', drop=-0.18823827350085606, haplotypes=c('blue'))

mutation.plot(33649037, 33649037, text='', col='black', drop=-0.1812533141437795, haplotypes=c('blue'))

mutation.plot(33649080, 33649080, text='', col='black', drop=-0.16973012702156548, haplotypes=c('blue'))

mutation.plot(33649146, 33649146, text='', col='black', drop=-0.19901104209586007, haplotypes=c('blue'))

mutation.plot(33649575, 33649575, text='', col='black', drop=-0.15141714981645876, haplotypes=c('blue'))

mutation.plot(33650866, 33650866, text='', col='black', drop=-0.18262431567343992, haplotypes=c('blue'))

mutation.plot(33651066, 33651066, text='', col='black', drop=-0.15156224513515215, haplotypes=c('blue'))

mutation.plot(33652198, 33652198, text='', col='black', drop=-0.12327041260139934, haplotypes=c('blue'))

mutation.plot(33652370, 33652373, text='', col='black', drop=-0.18740791228874207, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122980


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122980')

gene.transcript.model.plot(model=gene, gene_start=33300255, gene_bpstop=33301691, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122980')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33300255, gene_bpstop=33301691, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122980.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122980.1')

genemodel.plot(model=transcript, start=33300255, bpstop=33301691, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122980.1')
genemodel.plot_domain(model=t_domains, start=33300255, bpstop=33301691, orientation='reverse')

mutation.plot(33300272, 33300272, text='', col='black', drop=-0.18850058520539487, haplotypes=c('blue'))

mutation.plot(33300332, 33300332, text='', col='black', drop=-0.1803189312213642, haplotypes=c('blue'))

mutation.plot(33300334, 33300334, text='', col='black', drop=-0.16121415036753128, haplotypes=c('blue'))

mutation.plot(33300374, 33300374, text='', col='black', drop=-0.19694131714547594, haplotypes=c('blue'))

mutation.plot(33300655, 33300655, text='', col='black', drop=-0.18175976562046972, haplotypes=c('blue'))

mutation.plot(33300944, 33300944, text='', col='black', drop=-0.16196303977665588, haplotypes=c('blue'))

mutation.plot(33301668, 33301668, text='', col='black', drop=-0.16546171688968309, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123420')

gene.transcript.model.plot(model=gene, gene_start=36340164, gene_bpstop=36344124, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123420')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36340164, gene_bpstop=36344124, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123420.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.1')

genemodel.plot(model=transcript, start=36340164, bpstop=36344124, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.1')
genemodel.plot_domain(model=t_domains, start=36340164, bpstop=36344124, orientation='reverse')

mutation.plot(36340259, 36340259, text='', col='black', drop=-0.16476150372808085, haplotypes=c('blue'))

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.13061000294192004, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.14415701676941645, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1796976483972622, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.15636935973872892, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.13038269399301858, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.19598815367988925, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.11378684080490986, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.2')

genemodel.plot(model=transcript, start=36340251, bpstop=36344118, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.2')
genemodel.plot_domain(model=t_domains, start=36340251, bpstop=36344118, orientation='reverse')

mutation.plot(36340259, 36340259, text='', col='black', drop=-0.168580260312113, haplotypes=c('blue'))

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.18308788695990108, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.10892724429250453, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1702244103536497, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.19488823371469055, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.1817377432889074, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.1834670749417062, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.12335156426278963, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.3')

genemodel.plot(model=transcript, start=36340325, bpstop=36344102, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.3')
genemodel.plot_domain(model=t_domains, start=36340325, bpstop=36344102, orientation='reverse')

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.16778695899479457, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.1013665856562856, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1475123361634315, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.1838882664907598, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.15568866661883352, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.12204529511921719, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.11088666774661658, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.4')

genemodel.plot(model=transcript, start=36340349, bpstop=36344099, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.4')
genemodel.plot_domain(model=t_domains, start=36340349, bpstop=36344099, orientation='reverse')

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.1854045411628508, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.16323530982718143, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.12876452786583859, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.10831497959760364, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.11398656804927329, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.18111769054166063, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.12734830088783913, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.5')

genemodel.plot(model=transcript, start=36340349, bpstop=36344099, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.5')
genemodel.plot_domain(model=t_domains, start=36340349, bpstop=36344099, orientation='reverse')

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.1783171811640461, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.1279755410199696, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1460386899394618, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.13993035553185623, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.12302364854644804, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.12094351690316808, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.13943535729828269, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.6')

genemodel.plot(model=transcript, start=36340353, bpstop=36343960, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.6')
genemodel.plot_domain(model=t_domains, start=36340353, bpstop=36343960, orientation='reverse')

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.1186468364039821, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.13711002716358703, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1703629905935186, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.12225568386179239, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.1601108296830901, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.187244196212654, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.15675801320602958, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123420.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123420.7')

genemodel.plot(model=transcript, start=36340353, bpstop=36344099, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123420.7')
genemodel.plot_domain(model=t_domains, start=36340353, bpstop=36344099, orientation='reverse')

mutation.plot(36340373, 36340373, text='', col='black', drop=-0.16721134432495866, haplotypes=c('blue'))

mutation.plot(36341289, 36341289, text='', col='black', drop=-0.11366530516217452, haplotypes=c('blue'))

mutation.plot(36341534, 36341534, text='', col='black', drop=-0.1924733902505656, haplotypes=c('blue'))

mutation.plot(36341645, 36341645, text='', col='black', drop=-0.10118630365514268, haplotypes=c('blue'))

mutation.plot(36343629, 36343629, text='', col='black', drop=-0.1492838749744646, haplotypes=c('blue'))

mutation.plot(36343638, 36343638, text='', col='black', drop=-0.1787662250070973, haplotypes=c('blue'))

mutation.plot(36343792, 36343792, text='', col='black', drop=-0.1654593096768357, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123430


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123430')

gene.transcript.model.plot(model=gene, gene_start=36389182, gene_bpstop=36393618, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123430')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36389182, gene_bpstop=36393618, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123430.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123430.1')

genemodel.plot(model=transcript, start=36389182, bpstop=36393618, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123430.1')
genemodel.plot_domain(model=t_domains, start=36389182, bpstop=36393618, orientation='reverse')

mutation.plot(36389236, 36389236, text='', col='black', drop=-0.12050400934013679, haplotypes=c('blue'))

mutation.plot(36389359, 36389359, text='', col='black', drop=-0.19136082855555822, haplotypes=c('blue'))

mutation.plot(36389403, 36389406, text='', col='black', drop=-0.17278624033663442, haplotypes=c('blue'))

mutation.plot(36389475, 36389476, text='', col='black', drop=-0.17224222332846784, haplotypes=c('blue'))

mutation.plot(36389539, 36389539, text='', col='black', drop=-0.14828740526070333, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123430.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123430.2')

genemodel.plot(model=transcript, start=36389182, bpstop=36393618, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123430.2')
genemodel.plot_domain(model=t_domains, start=36389182, bpstop=36393618, orientation='reverse')

mutation.plot(36389236, 36389236, text='', col='black', drop=-0.1442768335645739, haplotypes=c('blue'))

mutation.plot(36389359, 36389359, text='', col='black', drop=-0.125781303646053, haplotypes=c('blue'))

mutation.plot(36389403, 36389406, text='', col='black', drop=-0.1525251892272236, haplotypes=c('blue'))

mutation.plot(36389475, 36389476, text='', col='black', drop=-0.17084212723211273, haplotypes=c('blue'))

mutation.plot(36389539, 36389539, text='', col='black', drop=-0.16676010463840196, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123250


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123250')

gene.transcript.model.plot(model=gene, gene_start=35489477, gene_bpstop=35490317, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123250')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35489477, gene_bpstop=35490317, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123250.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123250.1')

genemodel.plot(model=transcript, start=35489477, bpstop=35490317, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123250.1')
genemodel.plot_domain(model=t_domains, start=35489477, bpstop=35490317, orientation='reverse')

mutation.plot(35489658, 35489658, text='', col='black', drop=-0.3441764932074989, haplotypes=c('orange'))

mutation.plot(35489686, 35489686, text='', col='black', drop=-0.30973886138793927, haplotypes=c('orange'))

mutation.plot(35489826, 35489826, text='', col='black', drop=-0.3726665654272744, haplotypes=c('orange'))

mutation.plot(35490018, 35490018, text='', col='black', drop=-0.1257612364091464, haplotypes=c('blue'))

mutation.plot(35490067, 35490067, text='', col='black', drop=-0.10850807178544891, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123390


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123390')

gene.transcript.model.plot(model=gene, gene_start=36175890, gene_bpstop=36176983, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123390')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36175890, gene_bpstop=36176983, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123390.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123390.1')

genemodel.plot(model=transcript, start=36175890, bpstop=36176983, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123390.1

#Working on BaRT2v18chr3HG123150


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123150')

gene.transcript.model.plot(model=gene, gene_start=35151052, gene_bpstop=35152311, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123150')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35151052, gene_bpstop=35152311, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123150.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123150.1')

genemodel.plot(model=transcript, start=35151052, bpstop=35152311, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123150.1')
genemodel.plot_domain(model=t_domains, start=35151052, bpstop=35152311, orientation='reverse')

mutation.plot(35151670, 35151670, text='', col='black', drop=-0.17682900274852364, haplotypes=c('blue'))

mutation.plot(35152022, 35152022, text='', col='black', drop=-0.16197667958002204, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122960


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122960')

gene.transcript.model.plot(model=gene, gene_start=33290198, gene_bpstop=33290442, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122960')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33290198, gene_bpstop=33290442, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122960.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122960.1')

genemodel.plot(model=transcript, start=33290198, bpstop=33290442, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG122960.1

#Working on BaRT2v18chr3HG122900


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122900')

gene.transcript.model.plot(model=gene, gene_start=33170409, gene_bpstop=33174103, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122900')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33170409, gene_bpstop=33174103, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122900.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122900.1')

genemodel.plot(model=transcript, start=33170409, bpstop=33174103, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122900.1')
genemodel.plot_domain(model=t_domains, start=33170409, bpstop=33174103, orientation='reverse')

mutation.plot(33170465, 33170465, text='', col='black', drop=-0.16423222257409648, haplotypes=c('blue'))

mutation.plot(33170524, 33170524, text='', col='black', drop=-0.12023719923788595, haplotypes=c('blue'))

mutation.plot(33170569, 33170569, text='', col='black', drop=-0.16139002462517404, haplotypes=c('blue'))

mutation.plot(33170631, 33170631, text='', col='black', drop=-0.11491839460578385, haplotypes=c('blue'))

mutation.plot(33170675, 33170675, text='', col='black', drop=-0.1554535810660761, haplotypes=c('blue'))

mutation.plot(33170692, 33170692, text='', col='black', drop=-0.14283413347956192, haplotypes=c('blue'))

mutation.plot(33170724, 33170725, text='', col='black', drop=-0.17918089566420536, haplotypes=c('blue'))

mutation.plot(33170813, 33170813, text='', col='black', drop=-0.16353539636363648, haplotypes=c('blue'))

mutation.plot(33170844, 33170844, text='', col='black', drop=-0.10035765633330829, haplotypes=c('blue'))

mutation.plot(33171315, 33171315, text='', col='black', drop=-0.12161762257778048, haplotypes=c('blue'))

mutation.plot(33172481, 33172481, text='', col='black', drop=-0.12580449711127495, haplotypes=c('blue'))

mutation.plot(33172803, 33172803, text='', col='black', drop=-0.3875775958195137, haplotypes=c('orange'))

mutation.plot(33172883, 33172883, text='', col='black', drop=-0.14066246571855653, haplotypes=c('blue'))

mutation.plot(33173072, 33173072, text='', col='black', drop=-0.1658390761117391, haplotypes=c('blue'))

mutation.plot(33173930, 33173930, text='', col='black', drop=-0.11490432395203487, haplotypes=c('blue'))

mutation.plot(33173931, 33173931, text='', col='black', drop=-0.19220223434423603, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122900.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122900.2')

genemodel.plot(model=transcript, start=33170492, bpstop=33174006, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122900.2')
genemodel.plot_domain(model=t_domains, start=33170492, bpstop=33174006, orientation='reverse')

mutation.plot(33170524, 33170524, text='', col='black', drop=-0.10870893365667449, haplotypes=c('blue'))

mutation.plot(33170569, 33170569, text='', col='black', drop=-0.1105673673763998, haplotypes=c('blue'))

mutation.plot(33170631, 33170631, text='', col='black', drop=-0.11061839221153774, haplotypes=c('blue'))

mutation.plot(33170675, 33170675, text='', col='black', drop=-0.13438911280035234, haplotypes=c('blue'))

mutation.plot(33170692, 33170692, text='', col='black', drop=-0.13663152660108507, haplotypes=c('blue'))

mutation.plot(33170724, 33170725, text='', col='black', drop=-0.18044113677529716, haplotypes=c('blue'))

mutation.plot(33170813, 33170813, text='', col='black', drop=-0.3775984272371696, haplotypes=c('orange'))

mutation.plot(33170844, 33170844, text='', col='black', drop=-0.1476916567718874, haplotypes=c('blue'))

mutation.plot(33171315, 33171315, text='', col='black', drop=-0.10747671358346843, haplotypes=c('blue'))

mutation.plot(33172481, 33172481, text='', col='black', drop=-0.13190467677959178, haplotypes=c('blue'))

mutation.plot(33172803, 33172803, text='', col='black', drop=-0.34323811022357775, haplotypes=c('orange'))

mutation.plot(33172883, 33172883, text='', col='black', drop=-0.15591612376762637, haplotypes=c('blue'))

mutation.plot(33173072, 33173072, text='', col='black', drop=-0.14734038133573096, haplotypes=c('blue'))

mutation.plot(33173930, 33173930, text='', col='black', drop=-0.15349721014398965, haplotypes=c('blue'))

mutation.plot(33173931, 33173931, text='', col='black', drop=-0.1705943625107899, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122900.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122900.3')

genemodel.plot(model=transcript, start=33170492, bpstop=33173983, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122900.3')
genemodel.plot_domain(model=t_domains, start=33170492, bpstop=33173983, orientation='reverse')

mutation.plot(33170524, 33170524, text='', col='black', drop=-0.15034060692406984, haplotypes=c('blue'))

mutation.plot(33170569, 33170569, text='', col='black', drop=-0.18822191103703623, haplotypes=c('blue'))

mutation.plot(33170631, 33170631, text='', col='black', drop=-0.11305184438344851, haplotypes=c('blue'))

mutation.plot(33170675, 33170675, text='', col='black', drop=-0.12528058936121061, haplotypes=c('blue'))

mutation.plot(33170692, 33170692, text='', col='black', drop=-0.161978365015488, haplotypes=c('blue'))

mutation.plot(33170724, 33170725, text='', col='black', drop=-0.1043108885887035, haplotypes=c('blue'))

mutation.plot(33170813, 33170813, text='', col='black', drop=-0.17546893029283686, haplotypes=c('blue'))

mutation.plot(33170844, 33170844, text='', col='black', drop=-0.13464495078634764, haplotypes=c('blue'))

mutation.plot(33171315, 33171315, text='', col='black', drop=-0.14958412083882047, haplotypes=c('blue'))

mutation.plot(33172481, 33172481, text='', col='black', drop=-0.10711441540191013, haplotypes=c('blue'))

mutation.plot(33172803, 33172803, text='', col='black', drop=-0.3516806761303657, haplotypes=c('orange'))

mutation.plot(33172883, 33172883, text='', col='black', drop=-0.1674333868473237, haplotypes=c('blue'))

mutation.plot(33173072, 33173072, text='', col='black', drop=-0.17070206743478689, haplotypes=c('blue'))

mutation.plot(33173930, 33173930, text='', col='black', drop=-0.17548899351724806, haplotypes=c('blue'))

mutation.plot(33173931, 33173931, text='', col='black', drop=-0.1512538806168558, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122900.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122900.4')

genemodel.plot(model=transcript, start=33170641, bpstop=33174006, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122900.4')
genemodel.plot_domain(model=t_domains, start=33170641, bpstop=33174006, orientation='reverse')

mutation.plot(33170675, 33170675, text='', col='black', drop=-0.15199597278977592, haplotypes=c('blue'))

mutation.plot(33170692, 33170692, text='', col='black', drop=-0.1346794625888469, haplotypes=c('blue'))

mutation.plot(33170724, 33170725, text='', col='black', drop=-0.19293893655246314, haplotypes=c('blue'))

mutation.plot(33170813, 33170813, text='', col='black', drop=-0.3467951942957521, haplotypes=c('orange'))

mutation.plot(33170844, 33170844, text='', col='black', drop=-0.1946357925650829, haplotypes=c('blue'))

mutation.plot(33171315, 33171315, text='', col='black', drop=-0.1252849258647431, haplotypes=c('blue'))

mutation.plot(33172481, 33172481, text='', col='black', drop=-0.13506869030695715, haplotypes=c('blue'))

mutation.plot(33172803, 33172803, text='', col='black', drop=-0.35862079693975224, haplotypes=c('orange'))

mutation.plot(33172883, 33172883, text='', col='black', drop=-0.12076275805901808, haplotypes=c('blue'))

mutation.plot(33173072, 33173072, text='', col='black', drop=-0.17476652954135255, haplotypes=c('blue'))

mutation.plot(33173930, 33173930, text='', col='black', drop=-0.1569559404593775, haplotypes=c('blue'))

mutation.plot(33173931, 33173931, text='', col='black', drop=-0.18807429801966652, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123380


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123380')

gene.transcript.model.plot(model=gene, gene_start=36160347, gene_bpstop=36167544, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123380')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36160347, gene_bpstop=36167544, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123380.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.1')

genemodel.plot(model=transcript, start=36160347, bpstop=36167539, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.1')
genemodel.plot_domain(model=t_domains, start=36160347, bpstop=36167539, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.18569130779327647, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.1338381812409576, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.14895441911867463, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.18589602957824036, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.3859730482174878, haplotypes=c('orange'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.11599899065687883, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.3373481241985502, haplotypes=c('orange'))

mutation.plot(36166773, 36166773, text='', col='black', drop=-0.10211643359261714, haplotypes=c('blue'))

mutation.plot(36166948, 36166948, text='', col='black', drop=-0.1083497199516148, haplotypes=c('blue'))

mutation.plot(36167173, 36167195, text='', col='black', drop=-0.10508270321402349, haplotypes=c('blue'))

mutation.plot(36167252, 36167252, text='', col='black', drop=-0.1948726423667495, haplotypes=c('blue'))

mutation.plot(36167404, 36167404, text='', col='black', drop=-0.14077632106118163, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123380.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.2')

genemodel.plot(model=transcript, start=36160347, bpstop=36167543, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.2')
genemodel.plot_domain(model=t_domains, start=36160347, bpstop=36167543, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.14038017299078437, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.12645350625496427, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.16897528193703015, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.12775094497172534, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.35396614640630614, haplotypes=c('orange'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.10416386326792298, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.31273465858144434, haplotypes=c('orange'))

mutation.plot(36166492, 36166492, text='', col='black', drop=-0.14176429539478247, haplotypes=c('blue'))

mutation.plot(36166773, 36166773, text='', col='black', drop=-0.15956966449737459, haplotypes=c('blue'))

mutation.plot(36166948, 36166948, text='', col='black', drop=-0.11955654328390031, haplotypes=c('blue'))

mutation.plot(36167173, 36167195, text='', col='black', drop=-0.14817882758721715, haplotypes=c('blue'))

mutation.plot(36167252, 36167252, text='', col='black', drop=-0.12029461338671375, haplotypes=c('blue'))

mutation.plot(36167404, 36167404, text='', col='black', drop=-0.1360382483134755, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123380.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.3')

genemodel.plot(model=transcript, start=36160361, bpstop=36167539, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.3')
genemodel.plot_domain(model=t_domains, start=36160361, bpstop=36167539, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.13699412771561978, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.17628553953103282, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.19443587524698253, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.16450528197473904, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.14799208447957993, haplotypes=c('blue'))

mutation.plot(36163244, 36163244, text='', col='black', drop=-0.19463618085431353, haplotypes=c('blue'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.10814230425251684, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.37033988823254793, haplotypes=c('orange'))

mutation.plot(36166492, 36166492, text='', col='black', drop=-0.19633260719163526, haplotypes=c('blue'))

mutation.plot(36166773, 36166773, text='', col='black', drop=-0.13906341034351488, haplotypes=c('blue'))

mutation.plot(36166948, 36166948, text='', col='black', drop=-0.137553614817355, haplotypes=c('blue'))

mutation.plot(36167173, 36167195, text='', col='black', drop=-0.1982018580096953, haplotypes=c('blue'))

mutation.plot(36167252, 36167252, text='', col='black', drop=-0.13173624604609374, haplotypes=c('blue'))

mutation.plot(36167404, 36167404, text='', col='black', drop=-0.1790729437011476, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123380.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.4')

genemodel.plot(model=transcript, start=36160488, bpstop=36167544, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.4')
genemodel.plot_domain(model=t_domains, start=36160488, bpstop=36167544, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.17346639694538568, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.17332597435768896, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.19180760974768124, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.12889558060045225, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.1187158352952245, haplotypes=c('blue'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.1911526502968376, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.31088256826672367, haplotypes=c('orange'))

mutation.plot(36166773, 36166773, text='', col='black', drop=-0.1527970885478882, haplotypes=c('blue'))

mutation.plot(36166948, 36166948, text='', col='black', drop=-0.11393042061934516, haplotypes=c('blue'))

mutation.plot(36167173, 36167195, text='', col='black', drop=-0.10424212833387701, haplotypes=c('blue'))

mutation.plot(36167252, 36167252, text='', col='black', drop=-0.12025192339126191, haplotypes=c('blue'))

mutation.plot(36167404, 36167404, text='', col='black', drop=-0.1871874586677452, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123380.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.5')

genemodel.plot(model=transcript, start=36160617, bpstop=36165907, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.5')
genemodel.plot_domain(model=t_domains, start=36160617, bpstop=36165907, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.142763554340061, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.19783186564532976, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.10887822111156703, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.17742762542073573, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.13230162054420058, haplotypes=c('blue'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.1247149551406086, haplotypes=c('blue'))

mutation.plot(36164764, 36164764, text='', col='black', drop=-0.13695560797573036, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.38646504570863166, haplotypes=c('orange'))

#Working on BaRT2v18chr3HG123380.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123380.6')

genemodel.plot(model=transcript, start=36160634, bpstop=36165938, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123380.6')
genemodel.plot_domain(model=t_domains, start=36160634, bpstop=36165938, orientation='reverse')

mutation.plot(36160756, 36160758, text='', col='black', drop=-0.18543264574165827, haplotypes=c('blue'))

mutation.plot(36161744, 36161744, text='', col='black', drop=-0.12263780129908201, haplotypes=c('blue'))

mutation.plot(36161849, 36161849, text='', col='black', drop=-0.17297628333357506, haplotypes=c('blue'))

mutation.plot(36161852, 36161852, text='', col='black', drop=-0.10408654525755204, haplotypes=c('blue'))

mutation.plot(36162341, 36162341, text='', col='black', drop=-0.39739387467365156, haplotypes=c('orange'))

mutation.plot(36164336, 36164336, text='', col='black', drop=-0.1771723822606212, haplotypes=c('blue'))

mutation.plot(36165670, 36165670, text='', col='black', drop=-0.38625322998260553, haplotypes=c('orange'))

#Working on BaRT2v18chr3HG122970


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122970')

gene.transcript.model.plot(model=gene, gene_start=33294772, gene_bpstop=33299957, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122970')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33294772, gene_bpstop=33299957, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122970.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122970.1')

genemodel.plot(model=transcript, start=33294772, bpstop=33299932, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122970.1')
genemodel.plot_domain(model=t_domains, start=33294772, bpstop=33299932, orientation='reverse')

mutation.plot(33294871, 33294871, text='', col='black', drop=-0.15275997498483204, haplotypes=c('blue'))

mutation.plot(33294891, 33294891, text='', col='black', drop=-0.18858121986123838, haplotypes=c('blue'))

mutation.plot(33294929, 33294929, text='', col='black', drop=-0.16711221489292433, haplotypes=c('blue'))

mutation.plot(33295153, 33295153, text='', col='black', drop=-0.19739918383756466, haplotypes=c('blue'))

mutation.plot(33295177, 33295178, text='', col='black', drop=-0.1492735697341675, haplotypes=c('blue'))

mutation.plot(33295236, 33295236, text='', col='black', drop=-0.182723177690241, haplotypes=c('blue'))

mutation.plot(33296098, 33296098, text='', col='black', drop=-0.11812023101166405, haplotypes=c('blue'))

mutation.plot(33296420, 33296420, text='', col='black', drop=-0.11692903139198359, haplotypes=c('blue'))

mutation.plot(33297796, 33297796, text='', col='black', drop=-0.12692015667466944, haplotypes=c('blue'))

mutation.plot(33297999, 33297999, text='', col='black', drop=-0.10426615911189835, haplotypes=c('blue'))

mutation.plot(33298005, 33298005, text='', col='black', drop=-0.16798813892109327, haplotypes=c('blue'))

mutation.plot(33298006, 33298006, text='', col='black', drop=-0.16775435204876876, haplotypes=c('blue'))

mutation.plot(33298503, 33298503, text='', col='black', drop=-0.36641283806980013, haplotypes=c('orange'))

mutation.plot(33298672, 33298672, text='', col='black', drop=-0.31867245206854006, haplotypes=c('orange'))

mutation.plot(33298677, 33298677, text='', col='black', drop=-0.3332045100394102, haplotypes=c('orange'))

mutation.plot(33298718, 33298718, text='', col='black', drop=-0.10169790647453086, haplotypes=c('blue'))

mutation.plot(33298968, 33298968, text='', col='black', drop=-0.33027999167156935, haplotypes=c('orange'))

mutation.plot(33299021, 33299021, text='', col='black', drop=-0.17039726116833795, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122970.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122970.2')

genemodel.plot(model=transcript, start=33294814, bpstop=33299884, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122970.2')
genemodel.plot_domain(model=t_domains, start=33294814, bpstop=33299884, orientation='reverse')

mutation.plot(33294871, 33294871, text='', col='black', drop=-0.18915416698228638, haplotypes=c('blue'))

mutation.plot(33294891, 33294891, text='', col='black', drop=-0.1324317829445553, haplotypes=c('blue'))

mutation.plot(33294929, 33294929, text='', col='black', drop=-0.11419477959741134, haplotypes=c('blue'))

mutation.plot(33295153, 33295153, text='', col='black', drop=-0.1881263855245342, haplotypes=c('blue'))

mutation.plot(33295177, 33295178, text='', col='black', drop=-0.1059958073338636, haplotypes=c('blue'))

mutation.plot(33295236, 33295236, text='', col='black', drop=-0.13928751613441842, haplotypes=c('blue'))

mutation.plot(33296098, 33296098, text='', col='black', drop=-0.18858322537941696, haplotypes=c('blue'))

mutation.plot(33296420, 33296420, text='', col='black', drop=-0.17464065910875617, haplotypes=c('blue'))

mutation.plot(33296572, 33296572, text='', col='black', drop=-0.1788054340430557, haplotypes=c('blue'))

mutation.plot(33297796, 33297796, text='', col='black', drop=-0.11436861872941466, haplotypes=c('blue'))

mutation.plot(33297999, 33297999, text='', col='black', drop=-0.11213332061255521, haplotypes=c('blue'))

mutation.plot(33298005, 33298005, text='', col='black', drop=-0.1839711936971556, haplotypes=c('blue'))

mutation.plot(33298006, 33298006, text='', col='black', drop=-0.1498007908254324, haplotypes=c('blue'))

mutation.plot(33298503, 33298503, text='', col='black', drop=-0.39198627611876335, haplotypes=c('orange'))

mutation.plot(33298672, 33298672, text='', col='black', drop=-0.39471022164763303, haplotypes=c('orange'))

mutation.plot(33298677, 33298677, text='', col='black', drop=-0.34950187536809096, haplotypes=c('orange'))

mutation.plot(33298718, 33298718, text='', col='black', drop=-0.15506215985136465, haplotypes=c('blue'))

mutation.plot(33298968, 33298968, text='', col='black', drop=-0.3629705411230278, haplotypes=c('orange'))

mutation.plot(33299021, 33299021, text='', col='black', drop=-0.17682948008000365, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122970.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122970.3')

genemodel.plot(model=transcript, start=33295650, bpstop=33299957, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG122970.3')
genemodel.plot_domain(model=t_domains, start=33295650, bpstop=33299957, orientation='reverse')

mutation.plot(33296098, 33296098, text='', col='black', drop=-0.11532442761487278, haplotypes=c('blue'))

mutation.plot(33296420, 33296420, text='', col='black', drop=-0.14099200727150007, haplotypes=c('blue'))

mutation.plot(33297796, 33297796, text='', col='black', drop=-0.19431990693860557, haplotypes=c('blue'))

mutation.plot(33297999, 33297999, text='', col='black', drop=-0.1378044300437779, haplotypes=c('blue'))

mutation.plot(33298005, 33298005, text='', col='black', drop=-0.13982324642970434, haplotypes=c('blue'))

mutation.plot(33298006, 33298006, text='', col='black', drop=-0.10914973059265481, haplotypes=c('blue'))

mutation.plot(33298503, 33298503, text='', col='black', drop=-0.3603439590716656, haplotypes=c('orange'))

mutation.plot(33298672, 33298672, text='', col='black', drop=-0.3006535337176106, haplotypes=c('orange'))

mutation.plot(33298677, 33298677, text='', col='black', drop=-0.332192035425495, haplotypes=c('orange'))

mutation.plot(33298718, 33298718, text='', col='black', drop=-0.19301739436445078, haplotypes=c('blue'))

mutation.plot(33298968, 33298968, text='', col='black', drop=-0.3616059411514629, haplotypes=c('orange'))

mutation.plot(33299021, 33299021, text='', col='black', drop=-0.1667793914049024, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123050


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123050')

gene.transcript.model.plot(model=gene, gene_start=33873281, gene_bpstop=33876342, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123050')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33873281, gene_bpstop=33876342, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123050.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123050.1')

genemodel.plot(model=transcript, start=33873281, bpstop=33876342, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123050.1')
genemodel.plot_domain(model=t_domains, start=33873281, bpstop=33876342, orientation='reverse')

mutation.plot(33873335, 33873335, text='', col='black', drop=-0.18899447822961268, haplotypes=c('blue'))

mutation.plot(33873337, 33873337, text='', col='black', drop=-0.33060449690693283, haplotypes=c('orange'))

mutation.plot(33873343, 33873343, text='', col='black', drop=-0.329069932362595, haplotypes=c('orange'))

mutation.plot(33873492, 33873492, text='', col='black', drop=-0.3212106825435438, haplotypes=c('orange'))

mutation.plot(33873694, 33873694, text='', col='black', drop=-0.39876557206752966, haplotypes=c('orange'))

mutation.plot(33873715, 33873715, text='', col='black', drop=-0.38892237287954456, haplotypes=c('orange'))

mutation.plot(33873802, 33873802, text='', col='black', drop=-0.35691848604556414, haplotypes=c('orange'))

mutation.plot(33873988, 33873988, text='', col='black', drop=-0.3497477660567593, haplotypes=c('orange'))

mutation.plot(33875039, 33875039, text='', col='black', drop=-0.16463289111888213, haplotypes=c('blue'))

mutation.plot(33876153, 33876153, text='', col='black', drop=-0.12916734621698045, haplotypes=c('blue'))

mutation.plot(33876216, 33876217, text='', col='black', drop=-0.15730299034159662, haplotypes=c('blue'))

mutation.plot(33876241, 33876241, text='', col='black', drop=-0.13929480727075838, haplotypes=c('blue'))

mutation.plot(33876282, 33876282, text='', col='black', drop=-0.14380191880324578, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123480')

gene.transcript.model.plot(model=gene, gene_start=36576714, gene_bpstop=36583003, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123480')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36576714, gene_bpstop=36583003, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123480.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.1')

genemodel.plot(model=transcript, start=36576714, bpstop=36582957, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.1')
genemodel.plot_domain(model=t_domains, start=36576714, bpstop=36582957, orientation='reverse')

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.11624501458586961, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5814264050888577, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.13393025211749876, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.10


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.10')

genemodel.plot(model=transcript, start=36577190, bpstop=36582946, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.10')
genemodel.plot_domain(model=t_domains, start=36577190, bpstop=36582946, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.17130687502835965, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.13924152472711493, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.19311999492609316, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.1991748367493988, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.34318874088768087, haplotypes=c('orange'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5024520478320218, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.12143004219839626, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.11


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.11')

genemodel.plot(model=transcript, start=36577194, bpstop=36582946, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.11')
genemodel.plot_domain(model=t_domains, start=36577194, bpstop=36582946, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.11399603552151186, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.13058105275984527, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.18746605162160662, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.18634903427013869, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.13790809339683155, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5068382671329202, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.11752364049667677, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.12


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.12')

genemodel.plot(model=transcript, start=36577194, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.12')
genemodel.plot_domain(model=t_domains, start=36577194, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.18028996997394514, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.10982801098155562, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.16931488579693027, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.16881191958423103, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.13463822957221347, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5930310471310006, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.16995906511844736, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.13


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.13')

genemodel.plot(model=transcript, start=36577197, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.13')
genemodel.plot_domain(model=t_domains, start=36577197, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.13690730201484272, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.19358528348233456, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.11262668006224738, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.1749830072857782, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.1805184206798195, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5117518190919669, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.154453363992975, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.14


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.14')

genemodel.plot(model=transcript, start=36577197, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.14')
genemodel.plot_domain(model=t_domains, start=36577197, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.1202849878916532, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.10890586298925692, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.11792478726311247, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.10916441702048611, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.3773597596574741, haplotypes=c('orange'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5076316267370257, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.16922621494298612, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.15


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.15')

genemodel.plot(model=transcript, start=36577208, bpstop=36582961, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123480.15

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.18057690930159132, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.10085899260667816, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.1848087376287531, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.17221675970334419, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.16467329074969606, haplotypes=c('blue'))

mutation.plot(36580924, 36580924, text='', col='black', drop=-0.1449842673585365, haplotypes=c('blue'))

mutation.plot(36581027, 36581027, text='', col='black', drop=-0.184381339466619, haplotypes=c('blue'))

mutation.plot(36581207, 36581207, text='', col='black', drop=-0.1363083542754441, haplotypes=c('blue'))

mutation.plot(36581689, 36581689, text='', col='black', drop=-0.17704505213996274, haplotypes=c('blue'))

mutation.plot(36581814, 36581814, text='', col='black', drop=-0.1436485112411668, haplotypes=c('blue'))

mutation.plot(36581923, 36581923, text='', col='black', drop=-0.1425268512553187, haplotypes=c('blue'))

mutation.plot(36582042, 36582042, text='', col='black', drop=-0.33407737971559576, haplotypes=c('orange'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5236171437731161, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.11635587822287435, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.2')

genemodel.plot(model=transcript, start=36576991, bpstop=36583003, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123480.2

mutation.plot(36577111, 36577111, text='', col='black', drop=-0.16126261295358413, haplotypes=c('blue'))

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.174555200886927, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.16943850766885482, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.1994252587888371, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.178291963205941, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.14261560970214748, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5198033533762595, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.13376418622848812, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.3')

genemodel.plot(model=transcript, start=36577102, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.3')
genemodel.plot_domain(model=t_domains, start=36577102, bpstop=36582961, orientation='reverse')

mutation.plot(36577111, 36577111, text='', col='black', drop=-0.16424101638226576, haplotypes=c('blue'))

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.1618816176369109, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.16342153679734167, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.11987683715878734, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.11825099033942951, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.17136998905046097, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5792059191472125, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.11131019780567215, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.4


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.4')

genemodel.plot(model=transcript, start=36577131, bpstop=36582950, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.4')
genemodel.plot_domain(model=t_domains, start=36577131, bpstop=36582950, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.19416553270718356, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.111976808104204, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.19617661821952423, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.19747790387705957, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.19175166571429364, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5443839709059323, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.18219025464217292, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.5


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.5')

genemodel.plot(model=transcript, start=36577144, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.5')
genemodel.plot_domain(model=t_domains, start=36577144, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.1348158996146733, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.14363313484011375, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.18013958098506466, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.198521518313525, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.12620989730143714, haplotypes=c('blue'))

mutation.plot(36578975, 36578975, text='', col='black', drop=-0.12270130010004641, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5598220103850566, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.1803282047093398, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.6


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.6')

genemodel.plot(model=transcript, start=36577155, bpstop=36582950, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.6')
genemodel.plot_domain(model=t_domains, start=36577155, bpstop=36582950, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.1809115133713996, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.15527386028959628, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.17054263572920053, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.19087362463199575, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.13149501348694984, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5177664646801907, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.17481526781082932, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.7


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.7')

genemodel.plot(model=transcript, start=36577176, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.7')
genemodel.plot_domain(model=t_domains, start=36577176, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.17883314588604585, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.17480853231031884, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.17130301654770547, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.11100901007408547, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.11609246830003887, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5808920208821946, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.1401044520928067, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.8


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.8')

genemodel.plot(model=transcript, start=36577179, bpstop=36582961, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.8')
genemodel.plot_domain(model=t_domains, start=36577179, bpstop=36582961, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.19117335687429635, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.1212419294937919, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.13833587009330414, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.10854754366608124, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.12663015999170765, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5188516144686348, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.11479134968625647, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123480.9


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123480.9')

genemodel.plot(model=transcript, start=36577181, bpstop=36582950, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123480.9')
genemodel.plot_domain(model=t_domains, start=36577181, bpstop=36582950, orientation='reverse')

mutation.plot(36577209, 36577209, text='', col='black', drop=-0.17887360772455768, haplotypes=c('blue'))

mutation.plot(36577257, 36577260, text='', col='black', drop=-0.11864527994924474, haplotypes=c('blue'))

mutation.plot(36577301, 36577301, text='', col='black', drop=-0.10362865663859229, haplotypes=c('blue'))

mutation.plot(36577464, 36577464, text='', col='black', drop=-0.16452095493895053, haplotypes=c('blue'))

mutation.plot(36578235, 36578235, text='', col='black', drop=-0.16707556122468517, haplotypes=c('blue'))

mutation.plot(36578975, 36578975, text='', col='black', drop=-0.14604813674614434, haplotypes=c('blue'))

mutation.plot(36582714, 36582714, text='splice_donor&intron', col='black', drop=-0.5795843268307754, haplotypes=c('red'))

mutation.plot(36582792, 36582798, text='', col='black', drop=-0.18282951327264355, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123140


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123140')

gene.transcript.model.plot(model=gene, gene_start=35131940, gene_bpstop=35133679, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123140')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35131940, gene_bpstop=35133679, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123140.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123140.1')

genemodel.plot(model=transcript, start=35131940, bpstop=35133679, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123140.1')
genemodel.plot_domain(model=t_domains, start=35131940, bpstop=35133679, orientation='reverse')

mutation.plot(35132007, 35132008, text='', col='black', drop=-0.1050415101298875, haplotypes=c('blue'))

mutation.plot(35132279, 35132279, text='', col='black', drop=-0.1390371516224383, haplotypes=c('blue'))

mutation.plot(35132380, 35132381, text='', col='black', drop=-0.12151512919502346, haplotypes=c('blue'))

mutation.plot(35132434, 35132434, text='', col='black', drop=-0.19931463453648776, haplotypes=c('blue'))

mutation.plot(35133408, 35133408, text='', col='black', drop=-0.1862342817678904, haplotypes=c('blue'))

mutation.plot(35133421, 35133421, text='', col='black', drop=-0.3589790562654053, haplotypes=c('orange'))

mutation.plot(35133481, 35133482, text='frameshift', col='black', drop=-0.5906840510636303, haplotypes=c('red'))

mutation.plot(35133509, 35133509, text='', col='black', drop=-0.18680068751406545, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123140.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123140.2')

genemodel.plot(model=transcript, start=35132339, bpstop=35133618, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123140.2')
genemodel.plot_domain(model=t_domains, start=35132339, bpstop=35133618, orientation='reverse')

mutation.plot(35132380, 35132381, text='', col='black', drop=-0.15817906421164765, haplotypes=c('blue'))

mutation.plot(35132434, 35132434, text='', col='black', drop=-0.10434790336998305, haplotypes=c('blue'))

mutation.plot(35133408, 35133408, text='', col='black', drop=-0.11257383215431682, haplotypes=c('blue'))

mutation.plot(35133421, 35133421, text='', col='black', drop=-0.3310227320788934, haplotypes=c('orange'))

mutation.plot(35133481, 35133482, text='frameshift', col='black', drop=-0.5092688238878718, haplotypes=c('red'))

mutation.plot(35133509, 35133509, text='', col='black', drop=-0.15391879795430302, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123370


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123370')

gene.transcript.model.plot(model=gene, gene_start=36117124, gene_bpstop=36120695, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123370')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36117124, gene_bpstop=36120695, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123370.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123370.1')

genemodel.plot(model=transcript, start=36117124, bpstop=36120695, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123370.1')
genemodel.plot_domain(model=t_domains, start=36117124, bpstop=36120695, orientation='reverse')

#Working on BaRT2v18chr3HG123370.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123370.2')

genemodel.plot(model=transcript, start=36117220, bpstop=36120691, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123370.2')
genemodel.plot_domain(model=t_domains, start=36117220, bpstop=36120691, orientation='reverse')

#Working on BaRT2v18chr3HG122940


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122940')

gene.transcript.model.plot(model=gene, gene_start=33251832, gene_bpstop=33252065, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122940')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33251832, gene_bpstop=33252065, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122940.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122940.1')

genemodel.plot(model=transcript, start=33251832, bpstop=33252065, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG122940.1

#Working on BaRT2v18chr3HG123440


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123440')

gene.transcript.model.plot(model=gene, gene_start=36397204, gene_bpstop=36400472, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123440')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36397204, gene_bpstop=36400472, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123440.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123440.1')

genemodel.plot(model=transcript, start=36397204, bpstop=36400472, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123440.1')
genemodel.plot_domain(model=t_domains, start=36397204, bpstop=36400472, orientation='reverse')

mutation.plot(36397547, 36397555, text='', col='black', drop=-0.15154716247767247, haplotypes=c('blue'))

mutation.plot(36397563, 36397563, text='', col='black', drop=-0.16241615134138213, haplotypes=c('blue'))

mutation.plot(36397586, 36397587, text='', col='black', drop=-0.1858797257669309, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123440.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123440.2')

genemodel.plot(model=transcript, start=36397206, bpstop=36400472, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123440.2')
genemodel.plot_domain(model=t_domains, start=36397206, bpstop=36400472, orientation='reverse')

mutation.plot(36397547, 36397555, text='', col='black', drop=-0.16687226530550814, haplotypes=c('blue'))

mutation.plot(36397563, 36397563, text='', col='black', drop=-0.1316372372008316, haplotypes=c('blue'))

mutation.plot(36397586, 36397587, text='', col='black', drop=-0.1697142350801793, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123440.3


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123440.3')

genemodel.plot(model=transcript, start=36397217, bpstop=36400472, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123440.3')
genemodel.plot_domain(model=t_domains, start=36397217, bpstop=36400472, orientation='reverse')

mutation.plot(36397547, 36397555, text='', col='black', drop=-0.1927520063727298, haplotypes=c('blue'))

mutation.plot(36397563, 36397563, text='', col='black', drop=-0.192774609581771, haplotypes=c('blue'))

mutation.plot(36397586, 36397587, text='', col='black', drop=-0.11375501393149166, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG122990


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122990')

gene.transcript.model.plot(model=gene, gene_start=33302619, gene_bpstop=33306323, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122990')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33302619, gene_bpstop=33306323, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122990.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122990.1')

genemodel.plot(model=transcript, start=33302619, bpstop=33306323, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG122990.1

mutation.plot(33302658, 33302658, text='', col='black', drop=-0.16425003627726742, haplotypes=c('blue'))

mutation.plot(33302831, 33302831, text='', col='black', drop=-0.1421731381516533, haplotypes=c('blue'))

mutation.plot(33305702, 33305702, text='', col='black', drop=-0.3086661268754809, haplotypes=c('orange'))

mutation.plot(33305765, 33305765, text='', col='black', drop=-0.3831110880840275, haplotypes=c('orange'))

mutation.plot(33305766, 33305766, text='', col='black', drop=-0.17150567652303597, haplotypes=c('blue'))

mutation.plot(33305985, 33305985, text='', col='black', drop=-0.3969931054390413, haplotypes=c('orange'))

mutation.plot(33305990, 33305990, text='', col='black', drop=-0.13541608409203437, haplotypes=c('blue'))

mutation.plot(33306244, 33306244, text='', col='black', drop=-0.1941650289164058, haplotypes=c('blue'))

mutation.plot(33306245, 33306245, text='', col='black', drop=-0.3549854666979515, haplotypes=c('orange'))

mutation.plot(33306298, 33306298, text='', col='black', drop=-0.12843872187959757, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123240


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123240')

gene.transcript.model.plot(model=gene, gene_start=35414746, gene_bpstop=35415475, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123240')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35414746, gene_bpstop=35415475, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123240.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123240.1')

genemodel.plot(model=transcript, start=35414746, bpstop=35415475, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123240.1

#Working on BaRT2v18chr3HG123270


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123270')

gene.transcript.model.plot(model=gene, gene_start=36013350, gene_bpstop=36014191, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123270')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36013350, gene_bpstop=36014191, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123270.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123270.1')

genemodel.plot(model=transcript, start=36013350, bpstop=36014191, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123270.1

mutation.plot(36013467, 36013467, text='', col='black', drop=-0.1089548967545305, haplotypes=c('blue'))

mutation.plot(36013908, 36013908, text='', col='black', drop=-0.11978328698420941, haplotypes=c('blue'))

mutation.plot(36014042, 36014042, text='', col='black', drop=-0.1366573395679101, haplotypes=c('blue'))

mutation.plot(36014079, 36014079, text='', col='black', drop=-0.13902003814949185, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123400


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123400')

gene.transcript.model.plot(model=gene, gene_start=36310060, gene_bpstop=36311939, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123400')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=36310060, gene_bpstop=36311939, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123400.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123400.1')

genemodel.plot(model=transcript, start=36310060, bpstop=36311939, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123400.1')
genemodel.plot_domain(model=t_domains, start=36310060, bpstop=36311939, orientation='reverse')

mutation.plot(36310248, 36310248, text='', col='black', drop=-0.3028152961329473, haplotypes=c('orange'))

mutation.plot(36310302, 36310302, text='', col='black', drop=-0.3068883981642948, haplotypes=c('orange'))

mutation.plot(36310353, 36310353, text='', col='black', drop=-0.3768546118141769, haplotypes=c('orange'))

mutation.plot(36310405, 36310405, text='', col='black', drop=-0.17835877365091307, haplotypes=c('blue'))

mutation.plot(36310551, 36310551, text='', col='black', drop=-0.31900884944429536, haplotypes=c('orange'))

mutation.plot(36311012, 36311012, text='', col='black', drop=-0.3177999039530279, haplotypes=c('orange'))

#Working on BaRT2v18chr3HG123070


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123070')

gene.transcript.model.plot(model=gene, gene_start=33924212, gene_bpstop=33925756, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123070')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33924212, gene_bpstop=33925756, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123070.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123070.1')

genemodel.plot(model=transcript, start=33924212, bpstop=33925756, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123070.1

#Working on BaRT2v18chr3HG123070.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123070.2')

genemodel.plot(model=transcript, start=33924212, bpstop=33925742, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123070.2

#Working on BaRT2v18chr3HG122950


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG122950')

gene.transcript.model.plot(model=gene, gene_start=33260503, gene_bpstop=33260714, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG122950')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33260503, gene_bpstop=33260714, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG122950.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG122950.1')

genemodel.plot(model=transcript, start=33260503, bpstop=33260714, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG122950.1

#Working on BaRT2v18chr3HG123170


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123170')

gene.transcript.model.plot(model=gene, gene_start=35210478, gene_bpstop=35212372, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123170')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35210478, gene_bpstop=35212372, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123170.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123170.1')

genemodel.plot(model=transcript, start=35210478, bpstop=35212372, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123170.1')
genemodel.plot_domain(model=t_domains, start=35210478, bpstop=35212372, orientation='reverse')

#Working on BaRT2v18chr3HG123010


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123010')

gene.transcript.model.plot(model=gene, gene_start=33583827, gene_bpstop=33589403, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123010')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=33583827, gene_bpstop=33589403, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123010.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123010.1')

genemodel.plot(model=transcript, start=33583827, bpstop=33589403, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123010.1

mutation.plot(33583952, 33583958, text='', col='black', drop=-0.11276598888068722, haplotypes=c('blue'))

mutation.plot(33584530, 33584535, text='', col='black', drop=-0.17216965240941562, haplotypes=c('blue'))

mutation.plot(33584650, 33584650, text='', col='black', drop=-0.19288323586370978, haplotypes=c('blue'))

mutation.plot(33585221, 33585222, text='', col='black', drop=-0.1480872463583245, haplotypes=c('blue'))

mutation.plot(33585276, 33585276, text='', col='black', drop=-0.14339864347062153, haplotypes=c('blue'))

mutation.plot(33585379, 33585379, text='', col='black', drop=-0.10913820146467312, haplotypes=c('blue'))

mutation.plot(33585425, 33585425, text='', col='black', drop=-0.14188924185541146, haplotypes=c('blue'))

mutation.plot(33585446, 33585446, text='', col='black', drop=-0.18412076245017156, haplotypes=c('blue'))

mutation.plot(33585484, 33585484, text='', col='black', drop=-0.10489589020349009, haplotypes=c('blue'))

mutation.plot(33585917, 33585917, text='', col='black', drop=-0.17762260747009218, haplotypes=c('blue'))

mutation.plot(33586222, 33586222, text='', col='black', drop=-0.12103602836746274, haplotypes=c('blue'))

mutation.plot(33586252, 33586252, text='', col='black', drop=-0.16938237702123204, haplotypes=c('blue'))

mutation.plot(33587087, 33587087, text='', col='black', drop=-0.376782602950268, haplotypes=c('orange'))

mutation.plot(33587360, 33587360, text='', col='black', drop=-0.19681470711923083, haplotypes=c('blue'))

mutation.plot(33587459, 33587459, text='', col='black', drop=-0.13428564360053208, haplotypes=c('blue'))

mutation.plot(33588106, 33588106, text='', col='black', drop=-0.1305536610065578, haplotypes=c('blue'))

mutation.plot(33588226, 33588226, text='', col='black', drop=-0.14944642389810203, haplotypes=c('blue'))

mutation.plot(33589063, 33589063, text='', col='black', drop=-0.11109772923967765, haplotypes=c('blue'))

mutation.plot(33589215, 33589215, text='', col='black', drop=-0.15032941486123572, haplotypes=c('blue'))

mutation.plot(33589321, 33589322, text='', col='black', drop=-0.1892753245491747, haplotypes=c('blue'))

mutation.plot(33589322, 33589322, text='', col='black', drop=-0.16753277080737733, haplotypes=c('blue'))

mutation.plot(33589388, 33589388, text='', col='black', drop=-0.11161753609634731, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123010.2


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123010.2')

genemodel.plot(model=transcript, start=33584598, bpstop=33589345, orientation='reverse', xaxis=T)

#No domains for BaRT2v18chr3HG123010.2

mutation.plot(33584650, 33584650, text='', col='black', drop=-0.1572560434201225, haplotypes=c('blue'))

mutation.plot(33585221, 33585222, text='', col='black', drop=-0.1872064809500231, haplotypes=c('blue'))

mutation.plot(33585276, 33585276, text='', col='black', drop=-0.1279618884122628, haplotypes=c('blue'))

mutation.plot(33585379, 33585379, text='', col='black', drop=-0.12320919061654187, haplotypes=c('blue'))

mutation.plot(33585425, 33585425, text='', col='black', drop=-0.15822326867105677, haplotypes=c('blue'))

mutation.plot(33585446, 33585446, text='', col='black', drop=-0.11157541798527504, haplotypes=c('blue'))

mutation.plot(33585484, 33585484, text='', col='black', drop=-0.11620009635487173, haplotypes=c('blue'))

mutation.plot(33585917, 33585917, text='', col='black', drop=-0.16181027162964418, haplotypes=c('blue'))

mutation.plot(33586222, 33586222, text='', col='black', drop=-0.15204376104692208, haplotypes=c('blue'))

mutation.plot(33586252, 33586252, text='', col='black', drop=-0.16607821798918995, haplotypes=c('blue'))

mutation.plot(33587087, 33587087, text='', col='black', drop=-0.3679200643238334, haplotypes=c('orange'))

mutation.plot(33587360, 33587360, text='', col='black', drop=-0.10395373256861974, haplotypes=c('blue'))

mutation.plot(33587459, 33587459, text='', col='black', drop=-0.19280481803421554, haplotypes=c('blue'))

mutation.plot(33588106, 33588106, text='', col='black', drop=-0.15800855630198973, haplotypes=c('blue'))

mutation.plot(33589063, 33589063, text='', col='black', drop=-0.18685258223995072, haplotypes=c('blue'))

mutation.plot(33589215, 33589215, text='', col='black', drop=-0.197738056241758, haplotypes=c('blue'))

mutation.plot(33589321, 33589322, text='', col='black', drop=-0.13423748128769514, haplotypes=c('blue'))

mutation.plot(33589322, 33589322, text='', col='black', drop=-0.11424163642875895, haplotypes=c('blue'))

#Working on BaRT2v18chr3HG123220


gene <- transcripts %>% filter(gene == 'BaRT2v18chr3HG123220')

gene.transcript.model.plot(model=gene, gene_start=35388678, gene_bpstop=35389490, orientation='reverse', xaxis=T, gap=0.2)

#Note: gap parameter must be the same for both gene and g_domains

g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123220')
transcript_vector <- as.vector(unique(gene$transcript))
gene.transcript.model.plot_domain(model=g_domains, transcript_vector=transcript_vector, gene_start=35388678, gene_bpstop=35389490, orientation='reverse', gap=0.2)

#Working on BaRT2v18chr3HG123220.1


transcript <- transcripts %>% filter(transcript == 'BaRT2v18chr3HG123220.1')

genemodel.plot(model=transcript, start=35388678, bpstop=35389490, orientation='reverse', xaxis=T)

t_domains <- domains %>% filter(transcript == 'BaRT2v18chr3HG123220.1')
genemodel.plot_domain(model=t_domains, start=35388678, bpstop=35389490, orientation='reverse')

mutation.plot(35389119, 35389119, text='', col='black', drop=-0.3464586766017554, haplotypes=c('orange'))

mutation.plot(35389331, 35389331, text='', col='black', drop=-0.38451088559987534, haplotypes=c('orange'))

