library(ape)
tree1 <- read.tree("p_panAndOGASR/data/paspalum_anchors_ASTRALout.2023-08-23_noTdactm_SKHmodified_newLabeled.nwk")
tree1$edge.length = NULL

tree2 = drop.tip(tree1,c("tdacn2","tdacs2"))

tree3 = tree2
tree3$tip.label = c("Zea mays ssp. parviglumis - TIL11","Zea mays ssp. mays - B73","Zea mays ssp. parviglumis - TIL01",
                    "Zea mays ssp. mexicana - TIL25","Zea mays ssp. mexicana - TIL18","Zea mays ssp. huehuetenangensis",
                    "Zea nicaragensis", "Zea diploperennis - momo","Zea diploperennis - gigi",
                    "Tripsacum dactyloides - Southern","Tripsacum dactyloides - Northern","Tripsacum zopolitense",
                    "Urelytrum digitatum","Vossia cuspidata","Rhytacne rottboelloides","Rottboellia tuberculosa",
                    "Hemarthria compressa","Elionorus tripsacoides","Schizachyrium scoparium","Schizachyrium microstachyum",
                    "Andropogon virginicus","Andropogon chinensis","Andropogon gerardii","Cymbopogon refractus",
                    "Cymbopogon citratus","Heteropogon contortus","Themeda triandra","Bothriochloa laguroides",
                    "Pogonatherum paniceum","Sorghum bicolor","Ischaemum rugosum","Sorghastrum nutans",
                    "Andropogon tenuifolius","Thelepogon elegans","Chrysopogon serrulatus","Paspalum vaginatum","Setaria viridis")
edgeCol = apply(tree2$edge,1,function(x) any(x%in%c(grep("zm|zn|telega|sb|zTIL",tree2$tip.label),(8:13)+37)))
plot(tree3,tip.color = ifelse(grepl("zm|zn|telega|sb|zTIL",tree2$tip.label),"red","black"),edge.color = ifelse(edgeCol,"red","black"))

png("~/Dropbox/postDoc/projects/p_panAndOGASR/output/figure/FigA.png",width = 8.7,height = 8.7*0.75,units = "cm",res = 600,pointsize = 6)
plot(tree3,tip.color = ifelse(grepl("zm|zn|telega|sb|zTIL",tree2$tip.label),"orangered","black"),
     edge.color = ifelse(edgeCol,"orangered","black"),cex = 0.75,no.margin = T)
dev.off()
