library(tidyverse)
library(ggforce)
source("geom_hpline.R")
library(ggpubr)

rename_nested_columns <- function(.df) {
  .df %>% 
    rename_with(tolower) %>% 
    rename(cell = `...1`)
}

# Import and join data ----
# data_path <- "data/1_straight/"
data_path <- "data/2_background_corrected/"

dat <-
  tibble(
    files = list.files(data_path, full.names = T, pattern = ".csv"),
    data = map(files, read_csv)) %>%
  mutate(files = basename(files)) %>%
  separate(col = files,
           into = c("donor_and_location", "staining_1", "staining_2"),
           sep = '_', extra = 'drop') %>% 
  mutate(donor = parse_number(donor_and_location),
         location = str_sub(donor_and_location, -4L, -1L) %>% tolower(),
         staining_2 = str_remove(staining_2, ".csv")) %>% 
  select(-donor_and_location) %>% 
  relocate(donor, location, staining_1, staining_2, data) %>% 
  mutate(data = map(data, rename_nested_columns)) %>% 
  rowid_to_column("id")
  
dat  
summary(unnest(dat, data))


# Plots ----
datl <- unnest(dat, data) %>% 
  mutate(location = as_factor(location) %>% fct_relevel(c("prox", "dist"))) %>% 
  mutate(mean = abs(mean))
custom_colors <- c("black", "dodgerblue")


## Nucleus shape ----
### Size ----
nuc_size <- ggplot(datl, aes(location, area, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = bquote('Nuclaer size ['*μm^2*']'))
nuc_size


### Aspect ratio ----
nuc_ar <- ggplot(datl, aes(location, ar, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = "Nuclear aspect ratio")
nuc_ar


### Roundness ----
nuc_round <- ggplot(datl, aes(location, round, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = "Nuclear roundness")
nuc_round


### Cell count ----
cell_count <- datl %>% 
  group_by(id, donor, location) %>%
  summarise(n_cells = n())

nuc_count <- ggplot(cell_count, aes(location, n_cells, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = "Cell count [n]")
nuc_count


## Nucleus signal ----
### HIF1a concentration ----
nuc_hif_conc <- ggplot(datl, aes(location, mean, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  # scale_y_continuous(limits = y_limits) +
  labs(x = "Donor site",
       y = "Nuclear HIF1a")
nuc_hif_conc


### HIF1a amount ----
# This is not corrected for background signal; it's trash.

# nuc_hif_amount <- ggplot(datl, aes(location, intden, 
#                       group = location,
#                       color = location)) +
#   geom_boxplot(show.legend = FALSE) +
#   theme_classic() +
#   scale_color_manual(values = custom_colors) +
#   # scale_y_continuous(limits = y_limits) +
#   labs(x = "Donor site",
#        y = "Nuclear HIF1a amount [Σsignal]")
# nuc_hif_amount


### HIF1a distribution ----
nuc_hif_distr <- ggplot(datl, aes(mean, 
                      group = location,
                      color = location,
                      fill = location)) +
  geom_freqpoly(show.legend = FALSE) +
  # geom_histogram(aes(y=..density..), position = "identity", 
                 # alpha = 0.5, show.legend = FALSE) +
  # geom_density(fill = NA, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Nuclear HIF1a",
       y = "Count [n]")
nuc_hif_distr


# Correlation: roundness~HIF ----
cor_round_hif <- cor(datl$round, datl$mean)
cor.test(datl$round, datl$mean, method = "pearson")
shapiro.test(datl$round)
shapiro.test(datl$mean)

nuc_cor <- ggplot(datl, aes(round, mean)) +
  geom_point(aes(color = location), pch = 16, alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear roundness",
       y = "Nuclear HIF1a") +
  annotate(geom = "text", x = 0.2, y = 250, cex = 3,
           label = expression("r=0.48\np<0.0001"))
nuc_cor


# Combine plots ----
fig <- ggarrange(ggarrange(nuc_count, nuc_size, nrow = 1), 
                 ggarrange(nuc_hif_conc, nuc_hif_distr, nrow = 1),
                 ggarrange(nuc_round, nuc_cor, nrow = 1), 
                 nrow = 3)
annotate_figure(fig, 
                top = text_grob("Proximal vs. distal overview", 
                                face = "bold", 
                                size = 16))

# ggsave("results/nuclear_shape_vs_HIF_overview.svg",
       # width = 350, height = 400, dpi = 72, units = "px")


# Export results ----
dat_final <- select(datl, id, donor, location, cell, area, mean, round) %>% 
  rename(image_id = id, nuclear_size = area, nuclear_hif = mean, nuclear_roundness = round)

# write_csv(dat_final, "results/nuclear_shape_vs_HIF_results.csv")


