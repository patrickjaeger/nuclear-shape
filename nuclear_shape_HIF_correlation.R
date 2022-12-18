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
shapiro.test(datl$round)
shapiro.test(datl$mean)

hif_ar_cor <- ggplot(datl, aes(ar, mean)) +
  geom_point(aes(color = location), pch = 16, alpha = 0.3, show.legend = FALSE) +
  geom_smooth(formula = mean~749*exp(-0.2*ar), se = FALSE, color = "darkred") +
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
  x = seq(1, 20),
  y = 8000*(exp(-0.5*x+1))
)
plot(tt)
lines(tt)

### Fit model ----
hif_model <- nls(mean~I(a*exp(-b*ar)), datl, start = list(a = 4000, b=2), algorithm = "port", lower = c(3000, 0.5))
hif_model <- nls(mean~I(a*exp(-b*ar+1)), datl, start = list(a = 4000, b=2))
hif_model
exp_hif <- tibble(
  ar = seq(0,25, 0.2),
  mean = predict(hif_model, newdata = data.frame(ar = ar))
)
plot(c(0,25), c(0, 3000), type = "n")
lines(exp_hif)

ggplot(datl, aes(ar, mean)) +
  geom_point(aes(color = location), pch = 16, alpha = 0.3, show.legend = FALSE) +
  # geom_line(data = exp_hif, aes(ar, mean), color = "darkred") +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = custom_colors) +
  labs(x = "Nuclear AR",
       y = "Nuclear HIF1a")

## lm fit with log transform
hif_lm <- lm(log(mean)~ar, data = datl)
lm_coef <- coef(hif_lm)


dat_hif <- tibble(
  ar = seq(0,25, 0.2),
  # mean = predict(hif_lm, newdata = data.frame(ar = ar))
  mean = exp(lm_coef[1])*exp(lm_coef[2]*ar)
)
plot(c(0,25), c(0, 3000), type = "n")
lines(dat_hif)


# Combine plots ----
fig <- ggarrange(ggarrange(nuc_count, nuc_size, nrow = 1), 
                 ggarrange(nuc_round, nuc_ar, nrow = 1), 
                 ggarrange(nuc_hif_conc, nuc_hif_distr, nrow = 1),
                 nuc_cor,
                 nrow = 4)

annotate_figure(fig, 
                top = text_grob("Proximal vs. distal overview", 
                                face = "bold", 
                                size = 16))

# ggsave("results/221218_nuclear_shape_vs_HIF_overview.svg",
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
  labs(x = "Nuclear HIF1a",
       y = "Count [n]")

select(datl, -id, -donor, -location, -image, -cell) %>% 
  cor() %>% 
  corrplot::corrplot()

