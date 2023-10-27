library(ggplot2)

pData = data.frame(count = c(79,9,2464,1007,82,9),
                   test = factor(rep(c("Presence/absence","Sequence conservation","Stop codon conservation"),each = 2),
                                 levels = c("Presence/absence","Sequence conservation","Stop codon conservation")),
                   class = factor(rep(c("perennial","annual"),3),levels = c("perennial","annual")))

png("/workdir/sh2246/p_panAndOGASR/output/fig/FigB.png",width = 8.7*1.25,height = 8.7,units = "cm",res = 600,pointsize = 8)
ggplot(pData, aes(fill=class, y=count, x=class)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~test,scales = 'free') +
  scale_fill_grey(start = .3,end = .7) +
  theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 6),axis.title = element_text(size = 8),strip.text = element_text(size = 5)) +
  xlab("Presence/Conservation in") +
  ylab("Orthogroup counts")
dev.off()