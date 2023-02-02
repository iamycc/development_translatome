library(DESeq2)
library(ggplot2)
library(psych)
library(ggrepel)
library(ggpubr)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)

ribo=read.table('ribo_rawCount.txt')
ribo=ribo[,-1]

rna= read.table('rna_rawCount.txt')
rna=rna[,-1]




merge <- cbind(ribo,rna)
coldata=data.frame(Condition=rep(c(1,2,1,2), each=3),
                   SeqType=rep(c('RIBO','RNA'), each=6) )
rownames(coldata)=c(colnames(ribo), colnames(rna))
coldata$Condition <- factor(coldata$Condition)
coldata$SeqType <- factor(coldata$SeqType)

ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                 colData = coldata, 
                                 design =~ Condition + SeqType + Condition:SeqType)

ddsMat$SeqType = relevel(ddsMat$SeqType,"RNA")
ddsMat <- DESeq(ddsMat)
resultsNames(ddsMat)



# Choose the term you want to look at from resultsNames(ddsMat) 
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs
# Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(ddsMat, contrast=list("Condition2.SeqTypeRIBO"))
summary(res)
length(which(res$padj < 0.05))
write.table(res,"out/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)

all=read.csv('ALL_GSEA.csv',header = T)
rna_deg=read.table('out/rna_deg.txt')
ribo_deg=read.table('out/ribo_deg.txt')

te=read.table('out/deltaTE.txt')
te=te[!is.na(te$padj),]
one=function(n,all,rna_deg,ribo_deg,te,f=F,for_gene=F,s=NA){
  if (!for_gene) {gene=unlist(strsplit(all$intersections[n],split = ','))
  term_title=all$term_name[n]
  term_id=all$term_id[n]} else {gene=s
  term_title='user self'
  term_id='00000'}
  
  term_l=data.frame(rnaFC=rna_deg[gene,]$log2FoldChange,
                    riboFC=ribo_deg[gene,]$log2FoldChange)
  rownames(term_l)=gene
  term_l=na.omit(term_l)
  
  cor_r=corr.test(term_l$rnaFC, term_l$riboFC , method = "pearson",adjust = "fdr")$r
  cor_fdr=corr.test(term_l$rnaFC, term_l$riboFC , method = "pearson",adjust = "fdr")$p
  cor_r2=cor_r^2
  meanFCrna=mean(term_l$rnaFC)
  meanFCribo=mean(term_l$riboFC)
  te2=na.omit(te[gene,])
  meanTe=mean(te2$log2FoldChange)
  te_up=round(nrow(te2[te2$log2FoldChange >= 1,])/nrow(te2),2)
  te_down=round(nrow(te2[te2$log2FoldChange <= -1,])/nrow(te2),2)
  te_no=round(nrow(te2[abs(te2$log2FoldChange) < 1,])/nrow(te2),2)
  t=paste0('r2:',round(cor_r2,2),', fdr:', signif(cor_fdr,3))
  
  if (f) {ggplot(data = term_l,aes(x=rnaFC,y=riboFC))+
      geom_point()+
      ggtitle(label = paste0(term_title,'(',t,')' ))+
      geom_smooth(method = 'lm')+
      theme_pubr()} else {return(list(term_id = term_id, cor_r=cor_r, cor_fdr=cor_fdr, 
                                      cor_r2=cor_r2, meanFCrna=meanFCrna, meanFCribo=meanFCribo,
                                      meanTe=meanTe,te_up=te_up, te_down=te_down,
                                      te_no=te_no, n=nrow(te2)) )}
  
}
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

tePlot=function(gene,te,term){
  gene=te[gene,]
  gene=na.omit(gene)
  gene$la= ifelse( abs(gene$log2FoldChange) >= 1, ifelse(gene$log2FoldChange >= 1,'up','down'), 'no' )
  gene=gene[order(gene$log2FoldChange),]
  ggplot(gene,aes(x=1:nrow(gene),y=log2FoldChange,fill=la))+
    geom_area(alpha=0.35)+
    ggtitle(term)+
    scale_fill_manual(values=c("darkblue", "orange", "darkred"),
                      name="Experimental\nCondition",
                      breaks=c("down", "no", "up"),
                      labels=c("ES", "D2", "D5"))+theme_pubr()
}

reGene=function(all,term){
  tmp=all[all$term_name==term,]
  gene=unlist(strsplit(tmp$intersections,split = ','))
  return(gene)
}


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
