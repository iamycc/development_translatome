## Figure1
load('data1.RData') 
########## correlation analysis
library(PerformanceAnalytics)
chart.Correlation(one(dat = rna[,1:3],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = rna[,4:6],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = ribo[,1:3],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = ribo[,4:6],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')

##### 3nt periodicity analysis
library(ggplot2)
ggplot(data = esc_3nt[esc_3nt$from=='ribo',], aes(distance, per))+
  geom_bar(stat="identity", aes(fill=col))+
  geom_line(esc_3nt[esc_3nt$from=='rna',], mapping = aes(distance, per), 
            colour = 'purple')+
  facet_wrap(~type, scales = 'free_x')+
  scale_fill_brewer(palette = "Dark2")+
  geom_vline(xintercept = seq(-30,20,3)-0.5, colour = "grey80")+
  scale_y_continuous(expand = c(0,0.1))+
  scale_x_continuous(expand = c(0,1))+
  labs(x='CDS Positions [nt from start]',
       y="Read 5' ends [A.U.]",title = 'ES')+
  theme(axis.text = element_text(size = 14), 
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_line(colour = "grey80"), 
        axis.title = element_text(size = 14, face = "bold"), 
        panel.background = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        plot.margin=unit(rep(3,4),'lines'))

ggplot(data = npc_3nt[npc_3nt$from=='ribo',], aes(distance, per))+
  geom_bar(stat="identity", aes(fill=col))+
  geom_line(npc_3nt[npc_3nt$from=='rna',], mapping = aes(distance, per), 
            colour = 'purple')+
  facet_wrap(~type, scales = 'free_x')+
  scale_fill_brewer(palette = "Dark2")+
  geom_vline(xintercept = seq(-30,20,3)-0.5, colour = "grey80")+
  scale_y_continuous(expand = c(0,0.1))+
  scale_x_continuous(expand = c(0,1))+
  labs(x='CDS Positions [nt from start]',
       y="Read 5' ends [A.U.]",title = 'NPC')+
  theme(axis.text = element_text(size = 14), 
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_line(colour = "grey80"), 
        axis.title = element_text(size = 14, face = "bold"), 
        panel.background = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        plot.margin=unit(rep(3,4),'lines'))

## Figure2
load('data2.RData')
rna_tpm=rna_tpm[rowMeans(rna_tpm)>1,]
ribo_cpm=ribo_cpm[rowMeans(ribo_cpm)>1,]
gene=intersect(rownames(ribo_cpm),rownames(rna_tpm))

rna_deg=rna_deg[gene,c(1,2,4,8)]
ribo_deg=ribo_deg[gene,c(1,2,4,8)]

## Scatter plot of differentially expressed genes
ribo_up_down = get_up_down(data = ribo_deg, cutoff = 1.8)
rna_up_down = get_up_down(data = rna_deg, cutoff = 1.1)

library(ggrepel)
library(ggplot2)

ggplot(rna_up_down, aes(x=rank, y=log2FoldChange , color=log2FoldChange ))+
  geom_point()+
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  geom_vline(xintercept = c(sum(rna_up_down$la=='down')+0.5), linetype=2, size=0.25)+
  scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint=0)+
  geom_point(inherit.aes = F, data = rna_up_down[c('POU5F1','NANOG','PAX6','SOX1'),],
             aes(x=rank, y=log2FoldChange), size=3, color='black')+
  geom_point(inherit.aes = F, data = rna_up_down[c('POU5F1','NANOG','PAX6','SOX1'),],
             aes(x=rank, y=log2FoldChange), size=2, color='yellow')+
  geom_text_repel(inherit.aes = F, data = rna_up_down[c('POU5F1','NANOG','PAX6','SOX1'),], 
                  aes(x=rank, y=log2FoldChange, label=name), size=3)+
  xlab('rank of differentially expressed genes')+
  theme_bw() + 
  theme(panel.grid = element_line(color='white'), legend.title.align = 0.5,
        plot.margin=unit(rep(3,4),'lines'))

ggplot(ribo_up_down, aes(x=rank, y=log2FoldChange , color=log2FoldChange ))+
  geom_point()+
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  geom_vline(xintercept = c(sum(ribo_up_down$la=='down')+0.5), linetype=2, size=0.25)+
  scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint=0)+
  geom_point(inherit.aes = F, data = ribo_up_down[c('POU5F1','NANOG','PAX6','SOX1'),],
             aes(x=rank, y=log2FoldChange), size=3, color='black')+
  geom_point(inherit.aes = F, data = ribo_up_down[c('POU5F1','NANOG','PAX6','SOX1'),],
             aes(x=rank, y=log2FoldChange), size=2, color='yellow')+
  geom_text_repel(inherit.aes = F, data = ribo_up_down[c('POU5F1','NANOG','PAX6','SOX1'),], 
                  aes(x=rank, y=log2FoldChange, label=name), size=3)+
  xlab('rank of differentially expressed genes')+
  theme_bw() + 
  theme(panel.grid = element_line(color='white'), legend.title.align = 0.5,
        plot.margin=unit(rep(3,4),'lines'))


## combine analysis: RNA-seq & Ribo-seq
ggplot(data = comb, 
       mapping = aes(x=rna_log2FoldChange,y=ribo_log2FoldChange,color=la))+
  geom_point(size=4,alpha = 0.5)+
  scale_color_manual(values=c("red", "blue","purple","grey"),
                     name="Experimental\nCondition",
                     breaks=c("both", "ribo","rna",'no'))+
  ggtitle('RNA + RIBO')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        plot.title = element_text(size=18, face="bold", hjust = 0.5)
        #legend.position="none"
  )

## Figure3
load('data3.RData')
library(psych)
library(ggpubr)
library(ggrepel)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)

out=data.frame()
for (i in 1:nrow(all)) {
  term_id=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$term_id
  cor_r=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$cor_r
  cor_r2=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$cor_r2
  cor_fdr=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$cor_fdr
  meanFCrna=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$meanFCrna
  meanFCribo=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$meanFCribo
  meanTe=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$meanTe
  te_up=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$te_up
  te_down=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$te_down
  te_no=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$te_no
  geneN=one(n = i,all = all,rna_deg = rna_deg,ribo_deg = ribo_deg,te = te,f = F)$n
  out[i,'term_id']=term_id
  out[i,'cor_r']=cor_r
  out[i,'cor_r2']=cor_r2
  out[i,'cor_fdr']=cor_fdr
  out[i,'meanFCrna']=meanFCrna
  out[i,'meanFCribo']=meanFCribo
  out[i,'meanTe']=meanTe
  out[i,'te_up']=te_up
  out[i,'te_down']=te_down
  out[i,'te_no']=te_no
  out[i,'genes']=geneN
}

t=cbind(out,all[match(out$term_id,all$term_id),c(1,2,4,8)])
rownames(t)=t$term_name

ggplot(t,aes(x=te_no,y=cor_r2))+
  geom_point()+geom_smooth(method = 'lm')+
  geom_point(data = t[order(t$cor_r2)[c(1,4,8,14,17)],], size=4, alpha=0.3, color='blue')+
  geom_text_repel(data = t[order(t$cor_r2)[c(1,4,8,14,17)],],aes(label=term_name))+
  geom_point(data = t[order(t$cor_r2,decreasing = T)[c(3,6,7,11,15)],], size=4, alpha=0.3, color='blue')+
  geom_text_repel(data = t[order(t$cor_r2,decreasing = T)[c(3,6,7,11,15)],],aes(label=term_name))+
  theme_pubr()



pca <- PCA(t[,5:10], scale.unit = T, graph = F)
pcaPlot=pca$ind$coord
pcaPlot=cbind(pcaPlot[,1:2], t[match(rownames(pcaPlot), rownames(t)), 1:15] )
ggplot(data = pcaPlot,aes(x=Dim.1, y=Dim.2,color=meanTe,size=te_no) )+
  geom_point()+scale_colour_gradient2(low = 'blue',mid = 'orange',high = 'black')+
  scale_size_area()+
  geom_text_repel(data = pcaPlot[order(pcaPlot$cor_r2)[c(1,3)],],aes(label=term_name))+
  geom_text_repel(data = pcaPlot[order(pcaPlot$cor_r2,decreasing = T)[c(3)],],aes(label=term_name))+
  theme_bw()+
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=16),
    legend.text = element_text(size =16),
    legend.title = element_text(size =16 ,face="bold"),
    plot.margin=unit(rep(3,4),'lines'),
    plot.title = element_text(size=18, face="bold", hjust = 0.5),
    aspect.ratio=1)



tePlot(gene = neuralGene$V1, te = te, term = 'neural leading edge genes')
tePlot(gene = reGene(all=all,term = 'dendrite morphogenesis'), te = te, term = 'dendrite morphogenesis')
tePlot(gene = reGene(all=all,term = 'Wnt signaling pathway'), te = te, term = 'Wnt signaling pathway')
tePlot(gene = reGene(all=all,term = 'Cholesterol biosynthesis'), te = te, term = 'Cholesterol biosynthesis')
tePlot(gene = reGene(all=all,term = 'regulation of MAP kinase activity'), te = te, term = 'regulation of MAP kinase activity')
tePlot(gene = reGene(all=all,term = 'regulation of canonical Wnt signaling pathway'), te = te, term = 'regulation of canonical Wnt signaling pathway')
tePlot(gene = reGene(all=all,term = 'positive regulation of MAP kinase activity'), te = te, term = 'positive regulation of MAP kinase activity')
tePlot(gene = reGene(all=all,term = 'ECM-receptor interaction'), te = te, term = 'ECM-receptor interaction')
tePlot(gene = reGene(all=all,term = 'mesenchymal cell development'), te = te, term = 'mesenchymal cell development')
tePlot(gene = reGene(all=all,term = 'neural crest cell development'), te = te, term = 'neural crest cell development')
tePlot(gene = reGene(all=all,term = 'neural crest cell migration'), te = te, term = 'neural crest cell migration')
tePlot(gene = reGene(all=all,term = 'Transcriptional regulation of pluripotent stem cells'), te = te, term = 'Transcriptional regulation of pluripotent stem cells')

## Figure4
library(ggplot2)
library(ggpubr)
library(cowplot)
load('data4.RData')
p1=UTR_fun(data=utr_5, con='length', cell='esc', pre="5'UTR_length_es_te_class")
p2=UTR_fun(data=utr_5, con='length', cell='npc', pre="5'UTR_length_npc_te_class")

p3=UTR_fun(data=utr_5, con='GC', cell='esc', pre="5'UTR_GC%_es_te_class")
p4=UTR_fun(data=utr_5, con='GC', cell='npc', pre="5'UTR_GC%_npc_te_class")

p5=UTR_fun(data=utr_3, con='length', cell='esc', pre="3'UTR_length_es_te_class")
p6=UTR_fun(data=utr_3, con='length', cell='npc', pre="3'UTR_length_npc_te_class")

p7=UTR_fun(data=utr_3, con='GC', cell='esc', pre="3'UTR_GC%_es_te_class")
p8=UTR_fun(data=utr_3, con='GC', cell='npc', pre="3'UTR_GC%_npc_te_class")

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, labels=LETTERS[1:8])


## Pvalue
p1=UTR_P(data=utr_5, con='length', cell='esc', pre="5'UTR_length_es_te_class")
p2=UTR_P(data=utr_5, con='length', cell='npc', pre="5'UTR_length_npc_te_class")

p3=UTR_P(data=utr_5, con='GC', cell='esc', pre="5'UTR_GC%_es_te_class")
p4=UTR_P(data=utr_5, con='GC', cell='npc', pre="5'UTR_GC%_npc_te_class")

p5=UTR_P(data=utr_3, con='length', cell='esc', pre="3'UTR_length_es_te_class")
p6=UTR_P(data=utr_3, con='length', cell='npc', pre="3'UTR_length_npc_te_class")

p7=UTR_P(data=utr_3, con='GC', cell='esc', pre="3'UTR_GC%_es_te_class")
p8=UTR_P(data=utr_3, con='GC', cell='npc', pre="3'UTR_GC%_npc_te_class")
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, labels=LETTERS[1:8])



######################### The relationship between Kozak sequences and TE

p1=UTR_fun2(data=kozak[kozak$V2>0.8,], cell='esc', pre="kozak_es_te_class")
p2=UTR_fun2(data=kozak[kozak$V2>0.8,], cell='npc', pre="kozak_npc_te_class")
p3=UTR_P2(data=kozak[kozak$V2>0.8,], cell='esc', pre="5'UTR_length_es_te_class")
p4=UTR_P2(data=kozak[kozak$V2>0.8,], cell='npc', pre="5'UTR_length_npc_te_class")

plot_grid(p1, p2, p3, p4, nrow = 2, labels=LETTERS[1:4])


## codon analysis 
library(gridExtra)
library(ggpubr)
load('data5.RData')

cod_one=function(es_site, npc_site, codon, colours=c('red','blue')){
  one = function(dat){ dat=dat[dat$count>100,]
  dat=dat[,!colnames(dat) %in% c('count','TAA','TAG','TGA') ]
  print(nrow(dat))
  dat=dat/apply(dat,1,sum)
  return(dat)
  }
  es_site=one(es_site)
  npc_site=one(npc_site)
  es_con=data.frame(gene=rownames(es_site),
                    clusters='es',
                    freq=as.numeric(es_site[,colnames(es_site)==codon])
  )
  remove_extr=function(data,left_cutoff=0.1,right_cutoff=0.95){
    out=(data > quantile(data, left_cutoff)) & (data < quantile(data, right_cutoff))
    return(out)
  }
  es_con=es_con[remove_extr(es_con$freq),]
  npc_con=data.frame(gene=rownames(npc_site),
                     clusters='npc',
                     freq=as.numeric(npc_site[,colnames(npc_site)==codon])
  )
  npc_con=npc_con[remove_extr(npc_con$freq),]
  es_npc=rbind(es_con,npc_con)
  fc = mean(es_npc$freq[es_npc$clusters=='npc'])/mean(es_npc$freq[es_npc$clusters=='es'])
  fc = round(fc,2)
  p <- ggboxplot(es_npc, x = "clusters", y = "freq", color = "white")+
    geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
    stat_summary(fun=mean,orientation = 'x', geom="point", shape=21, size=3,  stroke = 1, fill="white")+
    geom_hline(yintercept = mean(es_npc$freq[es_npc$clusters=='es']), linetype = 2)+
    stat_compare_means(method = "kruskal.test", label.y = max(es_npc$freq)+0.03)+
    theme_bw()+expand_limits(y = c(0, max(es_npc$freq)+0.05)) +
    ggtitle(paste0(codon,' fc:',as.character(fc)))+scale_fill_manual(values=colours)+
    theme(
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      axis.text=element_text(size=16),
      # axis.title=element_text(size=16),
      axis.title=element_blank(),
      legend.text = element_text(size =16),
      legend.title=element_blank(),
      aspect.ratio=0.5,
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  return(p)
}


codonslist=c('GAC','GAT','TAT','TTA','ATT','CAT','AAA','AAG','ATG','TTT','AGG','AGA')
p <- list()
for (codons in codonslist) {
  p[[codons]] <- cod_one(
    es_site = es_asite,
    npc_site = npc_asite,
    codon = codons,
    colours = c('red','purple')
  )
}

do.call(grid.arrange,c(p, ncol=2))


ggbarplot(diff_codon, x="codon", y="NPC/ESC", fill = "log10P", color = "black",
          sort.by.groups = FALSE, x.text.angle=60, 
          ylab = "NPC / ESC", xlab = FALSE, legend.title=expression(-log[10](Pvalue)))+
  scale_fill_gradient(low=rgb(1,1,1,alpha=1),high = rgb(1,0,0,alpha = 0.5))+
  scale_y_log10(breaks = c(0.8, 1, 1.2))+
  geom_hline(yintercept = c(0.8,1.2),lty=2,color=c('blue','red'))








