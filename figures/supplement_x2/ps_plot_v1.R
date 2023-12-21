library(data.table)
library(tidyverse)
library(scales)
library(ggtext)
setwd("C:/Users/arun/OneDrive - Iowa State University/OrganizedDocuments/HuffordLab/PanAnd/phylostratr/phylostratr_results")
psdata <- fread("combined_ps.tsv")

colnames(psdata) <- c("Species", "Strata", "Genes", "Order")

ps.collapsed <- psdata %>% group_by(Species) %>%
  mutate(SppTotal = sum(Genes)) %>%
  mutate(
    NewStrata =
      case_when(
        Strata == "cellular_organisms" ~ "Cellular Organisms",
        Strata == "Eukaryota" ~ "Eukaryota",
        Strata == "Viridiplantae" ~ "Viridiplantae",
        Strata == "Streptophyta" ~ "Streptophyta",
        Strata == "Streptophytina" ~ "Streptophyta",
        Strata == "Embryophyta" ~ "Embryophyta",
        Strata == "Tracheophyta" ~ "Embryophyta",
        Strata == "Euphyllophyta" ~ "Embryophyta",
        Strata == "Spermatophyta"   ~ "Embryophyta",
        Strata == "Magnoliopsida" ~ "Embryophyta",
        Strata == "Mesangiospermae" ~ "Embryophyta",
        Strata == "Liliopsida" ~ "Liliopsida",
        Strata == "Petrosaviidae"  ~ "Liliopsida",
        Strata == "commelinids" ~ "Liliopsida",
        Strata == "Poales" ~ "Liliopsida",
        Strata == "crown_grasses" ~ "Crown Grasses",
        Strata == "PACMAD_clade" ~ "Crown Grasses",
        Strata == "BOP_clade" ~ "Crown Grasses",
        Strata == "Panicoideae" ~ "Crown Grasses",
        Strata == "Androp+Paspal" ~ "Crown Grasses",
        Strata == "Andropogoneae" ~ "Andropogoneae",
        Strata == "awned_Andropogoneae" ~ "Andropogoneae",
        Strata == "burmanicus+the_rest" ~ "Andropogoneae",
        Strata == "Apludinae+the_rest" ~ "Andropogoneae",
        Strata == "Ischaemininae+the_rest" ~ "Andropogoneae",
        Strata == "Saccharinae+the_rest" ~ "Andropogoneae",
        Strata == "Germainiinae+the_rest" ~ "Andropogoneae",
        Strata == "core_Andropogoneae" ~ "Andropogoneae",
        Strata == "Cymbo+DASH_clade" ~ "Andropogoneae",
        Strata == "BBS_&_chinensis" ~ "Andropogoneae",
        Strata == "Andropogoninae" ~ "Andropogoneae",
        Strata == "Anatherum_&_Schizac" ~ "Andropogoneae",
        Strata == "Schizachyrium" ~ "Andropogoneae",
        Strata == "Cymbopogon" ~ "Andropogoneae",
        Strata == "incertae_sedis" ~ "Andropogoneae",
        Strata == "BCD+Them_&_Het" ~ "Andropogoneae",
        Strata == "Them_&_Het" ~ "Andropogoneae",
        Strata == "awnless_Andropogoneae" ~ "Andropogoneae",
        Strata == "Trips_Rhyt" ~ "Andropogoneae",
        Strata == "Ratz-Elio" ~ "Andropogoneae",
        Strata == "Ratzeburgiinae" ~ "Andropogoneae",
        Strata == "Voss-Urely" ~ "Andropogoneae",
        Strata == "Rhytachninae" ~ "Andropogoneae",
        Strata == "Tripsacinae" ~ "Andropogoneae",
        Strata == "Zea" ~ "Tripsacinae",
        Strata == "Central_America" ~ "Tripsacinae",
        Strata == "Zea_mays" ~ "Tripsacinae",
        Strata == "parv-mays-mex" ~ "Tripsacinae",
        TRUE                      ~ "Species Specific"
      )
  ) %>%
  group_by(Species, NewStrata) %>%
  mutate(GroupedGenes = sum(Genes)) %>%
  select(c(Species, NewStrata, SppTotal, GroupedGenes)) %>%
  mutate(Percentage = round(100 * GroupedGenes / SppTotal, 2)) %>%
  distinct()

ps.collapsed$NewStrata <-
  factor(
    ps.collapsed$NewStrata,
    levels = c(
      "Cellular Organisms",
      "Eukaryota",
      "Viridiplantae",
      "Streptophyta",
      "Embryophyta",
      "Liliopsida",
      "Crown Grasses",
      "Andropogoneae",
      "Tripsacinae",
      "Species Specific"
    )
  )
ps.collapsed$Species <-
  factor(
    ps.collapsed$Species,
    levels = c(
      "aburma",
      "achine",
      "agerar",
      "avirgi",
      "blagur",
      "ccitra",
      "crefra",
      "cserru",
      "hconto",
      "irugos",
      "ppanic",
      "smicro",
      "snutan",
      "sscopa",
      "telega",
      "ttrian",
      "sbicol",
      "hcompr",
      "etrips",
      "rtuber",
      "rrottb",
      "udigit",
      "vcuspi",
      "tdacn1",
      "tdacn2",
      "tdacs1",
      "tdacs2",
      "zTIL01",
      "zTIL11",
      "zTIL18",
      "zTIL25",
      "zdgigi",
      "zdmomo",
      "zmhuet",
      "znicar",
      "zB73v5"
    )
  )
stratacolors <- c(
  "Cellular Organisms" = "#3dbabe",
  "Eukaryota" = "#d35238",
  "Viridiplantae" = "#6cb543",
  "Streptophyta" = "#8361cc",
  "Embryophyta" = "#d89a35",
  "Liliopsida" = "#6a88cb",
  "Crown Grasses" = "#aba74b",
  "Andropogoneae" = "#c867b5",
  "Tripsacinae" = "#59722e",
  "Species Specific" = "#b77549"
)

ggplot(ps.collapsed, aes(x = Species, y = GroupedGenes, fill = NewStrata)) +
  geom_bar(stat = 'identity') +
  labs(x = "", y = "") + theme_minimal() +
  scale_y_continuous(labels = label_comma()) + scale_fill_manual(values = stratacolors) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 12,
      face = "italic"
    ),
    strip.text = element_text(
      face = "bold",
      color = "gray35",
      hjust = 0,
      size = 10
    ),
    strip.background = element_rect(fill = "white", linetype = "blank"),
    legend.position = "none"
  ) +
  facet_wrap("NewStrata", scales = "free_y", ncol = 3) 

ggsave(
  "phylostrata-counts_barplot_supplement_option3.pdf",
  width = 20,
  height = 15,
  dpi = 300
)

ggplot(ps.collapsed, aes(x = Species, y = Percentage, fill = NewStrata)) +
  geom_bar(stat = 'identity') +
  labs(x = "", y = "") + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = stratacolors) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 12,
      face = "italic"
    ),
    strip.text = element_text(
      face = "bold",
      color = "gray35",
      hjust = 0,
      size = 12
    ),
    strip.background = element_rect(fill = "white", linetype = "blank"),
    legend.position = "none"
  ) +
  facet_wrap("NewStrata", scales = "free_y", ncol = 3) 

ggsave(
  "phylostrata_barplot_supplement_option2.pdf",
  width = 20,
  height = 15,
  dpi = 300
)
# distint colors for sorghum, and B73.v5 (not sequenced in this project)
labelColors <- c(rep("#454840", 16), "#ac0052", rep("#454840", 18), "#ac0052") 

# plot starta
names <- as.data.frame(levels(ps.collapsed$NewStrata))
colnames(names) <- "strata"
p <- ggplot() +
  geom_bar(
    aes(x = Species, y = Percentage, fill = NewStrata),
    position = position_stack(reverse = TRUE),
    data = ps.collapsed,
    stat = "identity",
    width = 1
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, .1))
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = stratacolors) +
  xlab("") +
  ylab("Genes percent") +
  theme(
    legend.position = "none",
    legend.text = element_markdown(size = rel(1.2)),
    panel.background = element_rect(color = "#FFFFFF", fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y = element_markdown(size = 12, face = "bold"),
    axis.text.y = element_markdown(size = 12, face = "bold"),
    axis.text.x = element_markdown(
      angle = 90,
      vjust = 1,
      hjust = 1,
      size = 12,
      color = labelColors
    ),
    plot.margin = unit(c(0, 2, 0, 0), "in"),
    axis.ticks.y = element_line(colour = "white", size = 0),
    axis.ticks.x = element_line(colour = "#222222", size = 0),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = ggtext::element_markdown(size = rel(1.2))
  ) +
  geom_text(
    data = filter(ps.collapsed, Species == "zB73v5"),
    aes(
      x = 36.5,
      y = Percentage,
      label = names$strata
    ),
    position = position_stack(vjust = 0.5),
    color = stratacolors,
    size = 4,
    face = "bold",
    hjust = 0
  )
ggsave(
  "phylostrata_layer_supplement_option1.pdf",
  plot = p,
  dpi = 900,
  width = 14,
  height = 8,
  limitsize = F
)
