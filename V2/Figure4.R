library(ggplot2)
library(ggpubr)
library(cowplot)

comb = read.table('out/comb.txt')
rna = read.table('out/rna_tpm.txt')

ribo = read.table('out/ribo_cpm.txt')

rna = rna[rownames(comb),]
ribo = ribo[rownames(comb),]
te=(ribo+1)/(rna+1)
colnames(te) = paste0(sapply(strsplit(colnames(te),split = '_'),'[',1),'_','te')

rna_ribo_te = data.frame(rna,ribo,te)
rna_ribo_te$name = rownames(rna_ribo_te)
rna_ribo_te$es_te_mean = rowMeans(rna_ribo_te[,13:15])
rna_ribo_te$npc_te_mean = rowMeans(rna_ribo_te[,16:18])
rna_ribo_te$es_te_class = ifelse(log2(rna_ribo_te$es_te_mean) < log2(quantile(rna_ribo_te$es_te_mean,0.1)),'low',
                                 ifelse(log2(rna_ribo_te$es_te_mean) > log2(quantile(rna_ribo_te$es_te_mean,0.9)),'high','none'))
rna_ribo_te$npc_te_class = ifelse(log2(rna_ribo_te$npc_te_mean) < log2(quantile(rna_ribo_te$npc_te_mean,0.1)),'low',
                                  ifelse(log2(rna_ribo_te$npc_te_mean) > log2(quantile(rna_ribo_te$npc_te_mean,0.9)),'high','none'))
rna_ribo_te$es_rna_tpmPlus1 = rowMeans(rna_ribo_te[,1:3])+1
rna_ribo_te$npc_rna_tpmPlus1 = rowMeans(rna_ribo_te[,4:6])+1
rna_ribo_te$es_ribo_cpmPlus1 = rowMeans(rna_ribo_te[,7:9])+1
rna_ribo_te$npc_ribo_cpmPlus1 = rowMeans(rna_ribo_te[,10:12])+1
te = rna_ribo_te
rm(comb, rna, ribo, rna_ribo_te)

utr_5 = read.table('utr_5.txt', header = T, check.names = F)
utr_3 = read.table('utr_3.txt', header = T, check.names = F)
kozak = read.table('kozak.txt', header = T, check.names = F)

UTR_fun = function(data, con, cell, pre=''){
  if (con == 'length') p = ggplot(data,aes(x=length))+scale_x_log10()
  if (con == 'GC') p = ggplot(data,aes(x=`GC%`))
  if (cell == 'esc') {high=te$name[te$es_te_class == 'high'];low=te$name[te$es_te_class == 'low']}
  if (cell == 'npc') {high=te$name[te$npc_te_class == 'high'];low=te$name[te$npc_te_class == 'low']}
  p=p+
    geom_histogram(aes(y = ..density..),bins = 50,fill='gray')+
    geom_histogram(data = data[data$gene %in% high,],aes(y = ..density..),bins = 50,fill='red',alpha=0.3)+
    geom_histogram(data = data[data$gene %in% low,],aes(y=..density..),bins = 50,fill='blue',alpha=0.3)+
    geom_density()+
    geom_density(data = data[data$gene %in% high,], color='red')+
    geom_density(data = data[data$gene %in% low,], color='blue')+
    theme_pubr()+
    ggtitle(pre)+
    theme(plot.margin=unit(c(1,2,1,2),'lines'))
  return(p)
}
UTR_fun2 = function(data, cell, pre=''){
  if (cell == 'esc') {high=te$name[te$es_te_class == 'high'];low=te$name[te$es_te_class == 'low']}
  if (cell == 'npc') {high=te$name[te$npc_te_class == 'high'];low=te$name[te$npc_te_class == 'low']}
  p=ggplot(data,aes(V2))+
    geom_histogram(aes(y = ..density..),bins = 50,fill='gray')+
    geom_histogram(data = data[data$gene %in% high,],aes(y = ..density..),bins = 50,fill='red',alpha=0.3)+
    geom_histogram(data = data[data$gene %in% low,],aes(y=..density..),bins = 50,fill='blue',alpha=0.3)+
    geom_density()+
    geom_density(data = data[data$gene %in% high,], color='red')+
    geom_density(data = data[data$gene %in% low,], color='blue')+
    theme_pubr()+
    ggtitle(pre)+
    scale_x_log10()+
    theme(plot.margin=unit(c(1,2,1,2),'lines'))
  return(p)
}
UTR_P = function(data, con, cell, pre=''){
  if (cell == 'esc') {high=te$name[te$es_te_class == 'high'];low=te$name[te$es_te_class == 'low']}
  if (cell == 'npc') {high=te$name[te$npc_te_class == 'high'];low=te$name[te$npc_te_class == 'low']}
  data$type=ifelse(data$gene %in% high,'high',
                   ifelse(data$gene %in% low,'low','no'))
  my_comparisons <- list( c("no", "high"), c("no", "low"), c("high", "low") )
  if (con == 'length') p=ggboxplot(data,x="type",y="length",color="type",
                                   palette=c(rgb(0,0,0,alpha=1),rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5)))+
    scale_y_log10()+
    stat_compare_means(comparisons = my_comparisons,label = "p.format")+
    geom_hline(yintercept = median(data$length), linetype = 2)+
    theme(plot.margin=unit(rep(1,4),'lines'))+
    ggtitle(pre)
  if (con == 'GC') p=ggboxplot(data, x = "type", y = "GC%",color ="type",
                               palette=c(rgb(0,0,0,alpha=1),rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5)))+
    stat_compare_means(comparisons = my_comparisons,label = "p.format")+
    geom_hline(yintercept = median(data$`GC%`), linetype = 2)+
    theme(plot.margin=unit(rep(1,4),'lines'))+
    ggtitle(pre)
  return(p)
}
UTR_P2 =function(data, cell, pre=''){
  if (cell == 'esc') {high=te$name[te$es_te_class == 'high'];low=te$name[te$es_te_class == 'low']}
  if (cell == 'npc') {high=te$name[te$npc_te_class == 'high'];low=te$name[te$npc_te_class == 'low']}
  data$type=ifelse(data$gene %in% high,'high',
                   ifelse(data$gene %in% low,'low','no'))
  my_comparisons <- list( c("no", "high"), c("no", "low"), c("high", "low") )
  p=ggboxplot(data, x = "type", y = "V2",color ="type",
              palette=c(rgb(0,0,0,alpha=1),rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5)))+
    stat_compare_means(comparisons = my_comparisons,label = "p.format")+
    geom_hline(yintercept = median(data$V2), linetype = 2)+
    theme(plot.margin=unit(rep(1,4),'lines'))+
    ggtitle(pre)
  return(p)
}

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


