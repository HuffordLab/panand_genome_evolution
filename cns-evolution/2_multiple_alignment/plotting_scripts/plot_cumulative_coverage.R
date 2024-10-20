library(ggplot2)

setwd("C:\\Users\\Armin\\ws\\panand\\analysis\\panand_cns_analysis\\2_multiple_alignment\\plotting_scripts")

cov_stats <- read.delim("..\\data\\panand_all_lt2.subfilt_10052023_clean_coverage_stats.csv",sep=',',header=TRUE)
cov_stats <- cov_stats[cov_stats$species!=1,]
cov_stats$species <- cov_stats$species -1 

cov_stats_cum <- cov_stats %>% map_df(rev) %>% mutate(csum=cumsum(sum_species_bp))

ggplot(cov_stats_cum ,aes(species,csum/1000000)) + 
  geom_area(fill="grey") +
  #geom_vline(xintercept=c(9), linetype="dotted")+
  #geom_hline(yintercept=c(271.686786), linetype="dotted")+
  #xlim(1,34)+
  scale_x_continuous(limits=c(1, 34), expand = c(0, 0))+
  scale_y_continuous(limits=c(0, 700), expand = c(0, 0))+
  ylab("Cumulative aligned sequence (Mbp)")+
  xlab("Minimum number of species")+
  theme_bw()+
  theme(text = element_text(size = 18))
ggsave("C:\\Users\\Armin\\ws\\panand\\analysis\\panand_cns_analysis\\2_multiple_alignment\\figures\\Fig_cumulative_aligned_mbp.pdf")

