library(ggrepel)
library(ggplot2)
library(DESeq2)

rna_count = read.table('rna_rawCount.txt')
ribo_count = read.table('ribo_rawCount.txt')


run_DEG_deseq2 = function(dat, group_list, sample_A, sample_B){
  #library(DESeq2)
  colData = data.frame(row.names=colnames(dat), 
                       group_list=group_list)
  colData$group_list = factor(colData$group_list)
  
  #group_list需要变成factor, 不然会报错
  dds = DESeqDataSetFromMatrix(countData = dat,
                               colData = colData,
                               design = ~group_list)
  dds = DESeq(dds)
  
  normalized_counts = counts(dds, normalized=TRUE)
  res = results(dds,contrast = c("group_list", sample_A, sample_B))
  baseA <-counts(dds, normalized=TRUE)[, colData(dds)$group_list == sample_A]
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <- sample_A
  baseB <- counts(dds, normalized=TRUE)[, colData(dds)$group_list == sample_B]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- sample_B
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
  res <- cbind(name=rownames(res), as.data.frame(res))
  res$padj[is.na(res$padj)] <- 1
  res <- res[order(res$pvalue),]
  return(list(deg=res, normalized_counts=normalized_counts))
}

ribo_deg = run_DEG_deseq2(dat = ribo_count[,-1], group_list = rep(c('es', 'npc'), each=3)
                          , sample_A = 'npc', sample_B = 'es')
ribo_deseq_count=data.frame(ribo_deg$normalized_counts)
identical(rownames(ribo_count),rownames(ribo_deseq_count))
ribo_deseq_count$len=ribo_count$len
ribo_deseq_count=ribo_deseq_count[,c(7,1:6)]
ribo_deg=ribo_deg$deg
ribo_deg=ribo_deg[,-1]


rna_deg = run_DEG_deseq2(dat = rna_count[,-1], group_list = rep(c('es', 'npc'), each=3)
                         , sample_A = 'npc', sample_B = 'es')
rna_deseq_count=data.frame(rna_deg$normalized_counts)
identical(rownames(rna_count),rownames(rna_deseq_count))
rna_deseq_count$len=rna_count$len
rna_deseq_count=rna_deseq_count[,c(7,1:6)]
rna_deg=rna_deg$deg
rna_deg=rna_deg[,-1]

getTPM=function(count,name=T){
  kb <- count[,1] / 1000
  count=count[,-1]
  rpk <- count / kb
  
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  cpm <- t(t(count)/colSums(count) * 1000000)
  fpkm <- t(t(rpk)/colSums(count) * 1000000)
  if(name){
    tpm=as.data.frame(tpm)
    cpm=as.data.frame(cpm)
    fpkm=as.data.frame(fpkm)
    len=length(colnames(tpm))
    tpm$name=rownames(tpm)
    cpm$name=rownames(cpm)
    fpkm$name=rownames(fpkm)
    tpm=tpm[,c(len+1,1:len)]
    cpm=cpm[,c(len+1,1:len)]
    fpkm=fpkm[,c(len+1,1:len)]
  }
  return(list(tpm=tpm,fpkm=fpkm,cpm=cpm))
}
ribo_cpm=getTPM(count = ribo_deseq_count)$cpm
ribo_cpm=ribo_cpm[,-1]
rna_tpm=getTPM(count = rna_deseq_count)$tpm
rna_tpm=rna_tpm[,-1]

write.table(ribo_cpm, 'out/ribo_cpm.txt', quote = F, sep = '\t')
write.table(rna_tpm, 'out/rna_tpm.txt', quote = F, sep = '\t')

rna_tpm=rna_tpm[rowMeans(rna_tpm)>1,]
ribo_cpm=ribo_cpm[rowMeans(ribo_cpm)>1,]
gene=intersect(rownames(ribo_cpm),rownames(rna_tpm))



rna_deg=rna_deg[gene,c(1,2,4,8)]
ribo_deg=ribo_deg[gene,c(1,2,4,8)]
write.table(rna_deg, 'out/rna_deg.txt', sep = '\t', quote = F)
write.table(ribo_deg, 'out/ribo_deg.txt', sep = '\t', quote = F)

## Scatter plot of differentially expressed genes
get_up_down = function(data,cutoff){
  data$la = ifelse(data$padj>0.05,
                   'no',
                   ifelse(data$log2FoldChange>=cutoff,
                          'up',
                          ifelse(data$log2FoldChange <= -cutoff,
                                 'down',
                                 'no'
                          )))
  data = data[order(data$log2FoldChange),]
  data = data[data$la!='no',]
  data$rank = 1:nrow(data)
  data$name=rownames(data)
  return(data)}
ribo_up_down = get_up_down(data = ribo_deg, cutoff = 1.8)
rna_up_down = get_up_down(data = rna_deg, cutoff = 1.1)



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

rna_deg= read.table('out/rna_deg.txt')
ribo_deg= read.table('out/ribo_deg.txt')
ven_coding =function(gene){
  cod=read.table('coding_list.txt')
  return(intersect(gene,cod$V1))
}

gene=ven_coding(rownames(rna_deg))
rna_deg=rna_deg[gene,]
ribo_deg=ribo_deg[gene,]

colnames(rna_deg)=paste0('rna','_',colnames(rna_deg))
colnames(ribo_deg)=paste0('ribo','_',colnames(ribo_deg))
identical(rownames(rna_deg),rownames(ribo_deg))


rna_deg=rna_deg[rna_deg$rna_padj<0.05,]
ribo_deg=ribo_deg[ribo_deg$ribo_padj<0.05,]
ribo_diffGENE=rownames(ribo_deg[abs(ribo_deg$ribo_log2FoldChange)>1.8,])
rna_diffGENE=rownames(rna_deg[abs(rna_deg$rna_log2FoldChange)>1.1,])
both=intersect(ribo_diffGENE,rna_diffGENE) 
ribo=setdiff(ribo_diffGENE,rna_diffGENE)   
rna=setdiff(rna_diffGENE,ribo_diffGENE)



rna_deg= read.table('out/rna_deg.txt')
ribo_deg= read.table('out/ribo_deg.txt')

colnames(rna_deg)=paste0('rna','_',colnames(rna_deg))
colnames(ribo_deg)=paste0('ribo','_',colnames(ribo_deg))

gene=ven_coding(rownames(rna_deg))
comb=cbind(rna_deg[gene,],ribo_deg[gene,])
comb$la='no'
comb$la[rownames(comb)%in%rna]='rna'
comb$la[rownames(comb)%in%ribo]='ribo'
comb$la[rownames(comb)%in%both]='both'
write.table(comb, 'out/comb.txt', sep = '\t', quote = F)

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
