# Literature Review
# Hengxing Zou

# packages

library(tidyverse)
library(magrittr)
library(readxl)
library(janitor)
library(patchwork)

# color scheme from https://jfly.uni-koeln.de/color/
color_scheme = c("#444444", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#bbbbbb")

##### Read File and Clean Up #####

lit = read_excel("SupportingData.xlsx")

##### Focal Taxa #####

# Change all input to lower case

lit %<>% mutate(`Focal Taxa` = tolower(`Focal Taxa`), `Secondary Taxa` = tolower(`Secondary Taxa`))

# Parse multiple focal taxa

lit_all_taxa = separate(lit, col = `Focal Taxa`, into = c("Focal Taxa 1", "Focal Taxa 2", "Focal Taxa 3"), sep = "&")
lit_all_taxa$`Focal Taxa 1` = str_trim(lit_all_taxa$`Focal Taxa 1`)
lit_all_taxa$`Focal Taxa 2` = str_trim(lit_all_taxa$`Focal Taxa 2`)
lit_all_taxa$`Focal Taxa 3` = str_trim(lit_all_taxa$`Focal Taxa 3`)

# Tally Taxa

focal_taxa_freq = 
  lit_all_taxa %>% 
  pivot_longer(cols = 4:6, names_to = "Type", values_to = "Focal Taxa") %>% 
  tabyl(`Focal Taxa`, show_na = F) %>% 
  arrange(n)
focal_taxa_freq

##### Experimental Design #####

# Temporal treatments

lit %<>% 
  mutate(`# Manipulated Species Binned` = case_when(`# Manipulated Species` %in% c("1", "2") ~ "1-2", 
                                                    `# Manipulated Species` %in% as.character(seq(3, 10, 1)) ~ "3-10", 
                                                    `# Manipulated Species` == ">10" ~ ">10",
                                                    TRUE ~ "No\ntemporal\ntreatment")) %>% 
  mutate(`# Manipulated Species Binned` = factor(`# Manipulated Species Binned`, 
                                                 levels = c("No\ntemporal\ntreatment", "1-2", "3-10", ">10")))

tabyl(lit, `# Manipulated Species Binned`)

temp = 
  lit %>% 
  ggplot(aes(x = `# Manipulated Species Binned`)) + 
  geom_bar(stat = "count", width = 0.6, fill = "#444444") + 
  geom_text(stat = "count", aes(label = ..count..), hjust = 2, color = "white", size = 6) + 
  labs(x = "Number of\nmanipulated species", y = "Number of studies") + 
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
  ggtitle("A") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15), axis.title.x = element_blank(), axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20), panel.grid.minor = element_blank())

# Temporal gradient

tabyl(lit, `Temporal Gradient?`)

##### Durations #####

durations = lit %>% 
  filter(!(`Duration (Generations)` %in% c("Mixed","Not Applicable"))) %>% 
  mutate(`Duration (Generations)` = factor(`Duration (Generations)`, 
                                           levels = c("<=1", "2-5", ">5")))

margin.table(table(durations$`Duration (Generations)`, durations$`Generation Time`), margin = c(1, 2))
prop.table(table(durations$`Duration (Generations)`, durations$`Generation Time`))

dura = 
  durations %>% 
  ggplot(aes(x = `Duration (Generations)`, fill = `Generation Time`)) + 
  geom_bar() + 
  geom_text(stat = "count", aes(label = ..count..), hjust = 2.5, color = "white", size = 6) + 
  # some hard coded values for the extra label
  geom_text(aes(x = 2, y = 32, label = 4), color = "#E69F00", size = 6) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  ylab("Number of studies") + 
  coord_flip() + 
  ggtitle("B") + 
  theme_bw() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), legend.position = "bottom", 
        plot.title = element_text(size = 20), panel.grid.minor = element_blank())

##### Taxa and Durations #####

# Including all taxa in the same study

dura_taxa_df = lit_all_taxa %>% 
  pivot_longer(cols = `Focal Taxa 1`: `Focal Taxa 3`, names_to = "# Taxa", values_to = "Focal Taxa") %>% 
  drop_na("Focal Taxa") %>% 
  filter(!(`Duration (Generations)` %in% c("Mixed","Not Applicable"))) %>% 
  mutate(`Duration (Generations)` = factor(`Duration (Generations)`, 
                                           levels = c("<=1", "2-5", ">5")))
  
# Broader categorization of taxa

dura_taxa_df %<>% 
  mutate(`Broad Taxa` = case_when(`Focal Taxa` %in% c("fish", "amphibians") ~ "Fish and Amphibians", 
                                  `Focal Taxa` %in% c("nematodes", "trematodes") ~ "Nematodes and Trematodes", 
                                  `Focal Taxa` %in% c("insects", "arachnids", "crustacean", "zooplankton", "crustaceans") ~ "Crustaceans", 
                                  `Focal Taxa` %in% c("terrestrial plants", "aquatic plants", "phytoplankton") ~ "Plants and Algae", 
                                  `Focal Taxa` %in% c("mycorrhizal/wood-decomposing fungi", "fungi") ~ "Fungi", 
                                  `Focal Taxa` %in% c("bacteria", "archaea") ~ "Prokaryotes", 
                                  `Focal Taxa` == "protists" ~ "Protists", 
                                  `Focal Taxa` == "virus" ~ "Virus", 
                                  TRUE ~ "Other taxa")) %>% 
  mutate(`Broad Taxa` = factor(`Broad Taxa`, 
                               levels = c("Other taxa", "Fish and Amphibians", "Plants and Algae", "Fungi", 
                                          "Crustaceans", "Protists", "Prokaryotes", "Nematodes and Trematodes", "Virus")))


dura_taxa = 
  dura_taxa_df %>% 
  ggplot(aes(x = `Duration (Generations)`, fill = `Broad Taxa`)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = color_scheme) +
  scale_y_continuous(name = "Proportion of studies", labels = c("0", "25%", "50%", "75%", "100%")) + 
  coord_flip() + 
  ggtitle("C") + 
  theme_bw() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15), legend.title = element_text(size = 20), legend.position = "bottom", 
        plot.title = element_text(size = 20), panel.grid.minor = element_blank())

##### Measurements #####

lit %<>% mutate(`Primary Measurement` = tolower(`Primary Measurement`))

lit_all_meas = separate(lit, col = `Primary Measurement`, into = c("Measurement 1", "Measurement 2", "Measurement 3"), sep = "&")
lit_all_meas$`Measurement 1` = str_trim(lit_all_meas$`Measurement 1`)
lit_all_meas$`Measurement 2` = str_trim(lit_all_meas$`Measurement 2`)
lit_all_meas$`Measurement 3` = str_trim(lit_all_meas$`Measurement 3`)

meas_freq = 
  lit_all_meas %>% 
  pivot_longer(cols = 15:17, names_to = "Type", values_to = "Measurements") %>% 
  tabyl(`Measurements`, show_na = F)
meas_freq

##### Figures #####

# Create directory

fig_dir = "Figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

# Result figure

ggsave("Figures/Fig1I.pdf", 
       plot = temp / dura / dura_taxa, 
       width = 3000, height = 4000, units = "px")
