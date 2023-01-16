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
data_path <- "data/3_batch_of_221218/csv/"


dat <-
tibble(files = list.files(data_path, full.names = T, pattern = ".csv"),
       data = map(files, read_csv)) %>%
  mutate(files = basename(files),
         files = str_remove(files, ".csv")) %>%
  separate(col = files,
           into = c("donor", "location", "image"),
           sep = '_', extra = 'drop') %>% 
  mutate(donor = parse_number(donor),
         location = tolower(location),
         image = ifelse(is.na(image), "0", image) %>% parse_number()) %>% 
  relocate(donor, location, image) %>% 
  mutate(data = map(data, rename_nested_columns),
         location = case_when(location == "prox" ~ "proximal",
                              location == "dist" ~ "distal",
                              TRUE ~ location)) %>% 
  rowid_to_column("id")

unique(dat$location)
dat
summary(unnest(dat, data))


## Data export for Greta ----
### by cell
# write_csv(datl, "results/230116_prox_vs_dist_by_cell.csv")

## by sample
dat_sample <- datl %>% 
  group_by(across(1:4)) %>% 
  summarise_all(mean) %>% 
  select(-cell)
dat_sample

# add cell count
cell_count <- datl %>% 
  group_by(across(1:4)) %>% 
  summarise(cell_count = n())
cell_count

dat_sample_with_count <- full_join(dat_sample, cell_count)

# write_csv(dat_sample_with_count, "results/230116_prox_vs_dist_by_sample.csv")

# Plots ----
datl <- unnest(dat, data) %>% 
  mutate(location = as_factor(location) %>% fct_relevel(c("proximal", "distal"))) %>% 
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
#### Boxplot ----
nuc_ar <- ggplot(datl, aes(location, ar, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = "Nuclear aspect ratio")
nuc_ar

#### Density ----
nuc_ar_dens <- ggplot(datl, aes(ar, 
                      group = location,
                      color = location)) +
  geom_density(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear aspect ratio",
       y = "Density")
nuc_ar_dens


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
#### Boxplot ----
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

##### Sample-wise ----
nuc_hif_conc_samplewise <- ggplot(datl, aes(location, mean, 
                      group = id,
                      color = location)) +
  # geom_boxplot(show.legend = FALSE) +
  stat_summary(geom = "point",
               position = position_jitter(width = 0.15),
               show.legend = FALSE,
               alpha = 0.5,
               pch = 16) +
  stat_summary(aes(group = location), 
               geom = "hpline",
               # fun = "median",
               show.legend = FALSE,
               # alpha = 0.7,
               width = 0.35) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  # scale_y_continuous(limits = y_limits) +
  labs(x = "Donor site",
       # y = bquote("Nuclear HIF1a [F"[nucleus]-"F"[background]~"]"))
       y = "Nuclear HIF1a")
nuc_hif_conc_samplewise


#### Density ----
nuc_hif_conc_dens <- ggplot(datl, aes(mean, 
                      group = location,
                      color = location)) +
  geom_density(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear HIF1a",
       y = "Density")
nuc_hif_conc_dens

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
  geom_point(aes(color = location), pch = 16, alpha = 0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear roundness",
       y = "Nuclear HIF1a")
  # annotate(geom = "text", x = 0.2, y = 250, cex = 3,
           # label = expression("r=0.48\np<0.0001"))
nuc_cor

# Correlation: ar~HIF ----
cor_round_hif <- cor(datl$ar, datl$mean)
cor.test(datl$ar, datl$mean, method = "pearson")
# shapiro.test(datl$round)
# shapiro.test(datl$mean)

hif_ar_cor <- ggplot(datl, aes(ar, mean)) +
  geom_point(aes(color = location), pch = 16, alpha = 0.3, show.legend = FALSE) +
  # geom_smooth(formula = mean~749*exp(-0.2*ar), se = FALSE, color = "darkred") +
  # geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear AR",
       y = "Nuclear HIF1a")
  # annotate(geom = "text", x = 0.2, y = 250, cex = 3,
           # label = expression("r=0.48\np<0.0001"))
hif_ar_cor


## Fit exponential model ----
### Estimate starting values ----
tt <- tibble(
  x = seq(0, 20),
  y = 2^(-x+2)
)
plot(tt)
lines(tt)
abline(v = 1)


### Fit model ----
# https://stackoverflow.com/questions/31851936/exponential-curve-fitting-in-r

#multiplicative error: we should use lm() on log-transformed data, because the error is constant on that scale
dat_mod <- select(datl, ar, mean) %>% mutate(log_mean = log(mean))

hif_model <- lm(log_mean~ar, dat_mod)
lm_coef <- coef(hif_model)
lm_coef

# formula with nls: nls(mean ~ a*exp(r*ar))
# linearized model: log(y) = log(a) + r*t, which is equal to: 
# Y = β0 + β1 * X, where β0 is our intercept and β1 our slope
# Therefore: intercept = log(a) & 

plot(dat_mod$ar, dat_mod$mean, col = alpha("black", 0.2))

hif_lm <- lm(dat_mod$log_mean ~ dat_mod$ar)
hif_lm_coef <- coef(hif_lm)
hif_lm_coef


lines(dat_mod$ar, exp(hif_lm_coef[1])*exp(hif_lm_coef[2]*dat_mod$ar), 
      col = "dodgerblue", lwd = 1)

hif_nls <- nls(mean ~ a*exp(r*ar), data = dat_mod, start = list(a = 3000, r = -0.2))
hif_nls_coef <- coef(hif_nls)
hif_nls_coef

lines(dat_mod$ar, nls_coef[1]*exp(nls_coef[2]*dat_mod$ar), 
      col = "darkred", lwd = 1)

# Conclusion: this data has too much noise to fit an exponential model


# Combine plots ----
fig <- ggarrange(ggarrange(nuc_count, nuc_size, nrow = 1), 
                 ggarrange(nuc_round, nuc_ar, nrow = 1), 
                 ggarrange(nuc_hif_conc, nuc_hif_distr, nuc_hif_conc_samplewise, nrow = 1),
                 nuc_cor,
                 nrow = 4)

annotate_figure(fig, 
                top = text_grob("Proximal vs. distal overview", 
                                face = "bold", 
                                size = 16))

# ggsave("results/221228_nuclear_shape_vs_HIF_overview.svg",
       # width = 350, height = 600, dpi = 72, units = "px")


# Export results ----
dat_final <- select(datl, id, donor, location, cell, area, mean, round) %>% 
  rename(image_id = id, nuclear_size = area, nuclear_hif = mean, nuclear_roundness = round)

# write_csv(dat_final, "results/nuclear_shape_vs_HIF_results.csv")


# Summmary ----
datl %>% 
  group_by(id, donor, image) %>% 
  summarise(mean = mean(mean)) %>% print(n = nrow(.))


ggplot(datl, aes(round, 
                 group = location,
                 color = location)) +
  geom_density(show.legend = FALSE) +
  # geom_histogram(aes(y=..density..), position = "identity", 
  # alpha = 0.5, show.legend = FALSE) +
  # geom_density(fill = NA, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Roundness",
       y = "Count [n]")

select(datl, -id, -donor, -location, -image, -cell) %>% 
  cor() %>% 
  corrplot::corrplot()


# HIF~AR normalized -------------------------------------------------------

dat2 <- datl %>% mutate(hif_score = mean/ar)

ggplot(dat2, aes(hif_score, 
                 group = location,
                 color = location)) +
  geom_density(show.legend = FALSE) +
  # geom_histogram(aes(y=..density..), position = "identity", 
  # alpha = 0.5, show.legend = FALSE) +
  # geom_density(fill = NA, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Roundness",
       y = "Count [n]")


ggplot(dat2, aes(location, hif_score, 
                      group = location,
                      color = location)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Donor site",
       y = "the fuck is this?")
