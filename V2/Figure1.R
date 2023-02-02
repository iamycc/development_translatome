library(PerformanceAnalytics)
library(ggplot2)

one = function(dat, cut_mean=5, logg=F){
  dat$mean=apply(dat,1,mean)
  dat=dat[dat$mean>cut_mean,]
  if (logg) dat=log(dat+1,logg)
  dat=dat[,!colnames(dat) == "mean"]
  return(dat)
}


rna = read.table('rna_rawCount.txt')
rna=rna[,-1]
ribo = read.table('ribo_rawCount.txt')
ribo=ribo[,-1]
esc_3nt = read.table('esc_3nt.txt', header = T)
esc_3nt$col = factor(esc_3nt$col, levels = c(2, 1, 0))

npc_3nt = read.table('npc_3nt.txt', header = T)
npc_3nt$col = factor(npc_3nt$col, levels = c(2, 1, 0))
########## correlation analysis

chart.Correlation(one(dat = rna[,1:3],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = rna[,4:6],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = ribo[,1:3],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')
chart.Correlation(one(dat = ribo[,4:6],cut_mean = 5,logg = 2), 
                  histogram=TRUE, pch=19, method = 'spearman')

##### 3nt periodicity analysis

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
