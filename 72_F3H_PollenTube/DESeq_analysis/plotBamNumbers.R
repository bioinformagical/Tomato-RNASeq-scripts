
suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(readr)
  library(tidyverse)
})

read_counts <- read_tsv("~/Downloads/bamnumbers-Sheet.tsv", col_types = "cnn") %>%
  rename(Aligned = ALIGNED, Unaligned = UNALIGNED) %>%
  drop_na() %>%
  mutate(Sample = str_remove_all(Sample, "[ │·]")) %>%
  pivot_longer(cols = 2:3, names_to = "Type", values_to = "Reads", names_transform = list(Type = ~ readr::parse_factor(.x, levels = c("Aligned", "Unaligned")))) %>%
  separate_wider_delim(Sample, "_", names = c("Junk1", "Sample")) %>%
  mutate(Sample = as.numeric(str_remove_all(Sample, "S")))

plot <- ggplot(read_counts, aes(x= Sample, y = Reads, fill = Type)) +
  geom_col(position = "stack") +
  scale_x_discrete() +
  scale_y_continuous(labels = scales::comma, expand = c(0,0)) +
  scale_fill_manual(values = c("orange","yellow" )) + 
  theme_minimal_hgrid()
#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot
