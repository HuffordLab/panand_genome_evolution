setwd(
  "C:/Users/arun/OneDrive - Iowa State University/OrganizedDocuments/HuffordLab/PanAnd/annotation-stats"
)
library(tidyverse)
library(ggalt)
lai <-
  read.csv("newLAI.tsv",
           sep = "\t",
           stringsAsFactors = TRUE)
lai$diff <- lai$RawLAI - lai$LAI
lai$Tech <- paste0(lai$data, "(", lai$Type, ")")
lai$newSpp <-
  paste0 (lai$spp, " (", format(round(lai$genomeSize / 1000, 2), nsmall = 2), " Gb)")




g <- ggplot() + geom_segment(
  data = lai,
  aes(
    y = newSpp,
    yend = newSpp,
    x = 1,
    xend = 46
  ),
  color = "#b2b2b2",
  size = 0.15
) +
  geom_dumbbell(
    data = lai,
    aes(y = newSpp, x = LAI, xend = RawLAI),
    size = 1.5,
    color = "#b2b2b2",
    size_x = 3,
    size_xend = 3,
    colour_x = "darkorchid",
    colour_xend = "royalblue"
  ) +
  geom_text(
    data = filter(lai, spp == "Zea-mays-ssp.mays-B73.v5"),
    aes(x = LAI, y = newSpp, label = "LAI"),
    color = "darkorchid",
    vjust = -2.5, hjust = 2,
    fontface = "bold"
  ) +
  geom_text(
    data = filter(lai, spp == "Zea-mays-ssp.mays-B73.v5"),
    aes(x = RawLAI, y = newSpp, label = "LAI (raw)"),
    color = "royalblue",
    vjust = -2.5,hjust = -0.5,
    fontface = "bold"
  ) +
  geom_text(
    data = lai,
    aes(x = RawLAI, y = newSpp, label = RawLAI),
    color = "royalblue",
    size = 2.75,
    hjust = -0.5,
    
  ) +
  geom_text(
    data = lai,
    aes(x = LAI, y = newSpp, label = LAI),
    color = "darkorchid",
    size = 2.75,
    hjust = 1.5
  )   + geom_text(
    data = lai,
    aes(
      x = (LAI + RawLAI) / 2,
      y = newSpp,
      label = Tech
    ),
    color = "slateblue",
    size = 2.75,
    vjust = -0.75
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 46)) +
  scale_y_discrete(expand = c(0.075, 0))  +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic")
  )
ggsave(
  "LAI_scores_panand_genomes_supplement.pdf",
  g,
  width = 8,
  height = 9,
  dpi = 300
)