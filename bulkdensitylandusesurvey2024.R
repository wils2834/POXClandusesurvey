# Bulk Density — 2024 CRP Project
library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)

# ---- Factor levels ------------------------------------------------------
type_levels <- c("Conventional","Soil_Health","Rotational_Grazing","New_CRP","Old_CRP")

type_labels <- c(
  "Conventional"       = "Conventional Cropping & Tillage",
  "Soil_Health"        = "No Till and/or Cover Crops",
  "Rotational_Grazing" = "Rotational Grazing",
  "New_CRP"            = "0-10 years CRP",
  "Old_CRP"            = "10+ years CRP"
)

fill_colors <- c(
  "Conventional"       = "#E69F00",
  "Soil_Health"        = "#56B4E9",
  "Rotational_Grazing" = "#009E73",
  "New_CRP"            = "#F0E442",
  "Old_CRP"            = "#0072B2"
)

# ---- Data import --------------------------------------------------------
bulk_density_driftless <- read.csv(
  "C:\\Users\\swils\\Downloads\\01 Projects\\2024 CRP\\Wilson Bulk Density 2024 - Sheet1.csv"
)

## ---- Data cleaning ------------------------------------------------------
bulk_density_driftless <- bulk_density_driftless %>%
  dplyr::rename(Bulk_Density = Bulk.Density) %>%
  filter(!is.na(Bulk_Density), !is.na(Type), !is.na(Farmer)) %>%
  mutate(
    Farmer = dplyr::recode(Farmer,
                           "Melissa_Wiest2"        = "Melissa_Wiest_2",
                           "Maurice_McLean_second" = "Maurice_McLean",
                           "Steve_Collins_second"  = "Steve_Collins",
                           "RLH_Inc_2_second"      = "RLH_Inc_2"
    ),
    Type = factor(Type, levels = type_levels),
    Farmer = as.factor(Farmer),
    County = as.factor(County)
  )
bulk_density_driftless <- bulk_density_driftless %>%
  mutate(Field = paste(Farmer, Type, sep = "_"))

n_distinct(bulk_density_driftless$Field)
cat("Rows after cleaning:", nrow(bulk_density_driftless), "\n")
cat("Distinct fields:", n_distinct(bulk_density_driftless$Field), "(expect 31)\n")
cat("Distinct UIDs:", n_distinct(bulk_density_driftless$UID), "\n\n")
stopifnot(n_distinct(bulk_density_driftless$Field) == 31)

# ---- Mixed models -------------------------------------------------------
model_BD     <- lmer(Bulk_Density      ~ Type + (1 | Field), data = bulk_density_driftless)
model_BD_log <- lmer(log(Bulk_Density) ~ Type + (1 | Field), data = bulk_density_driftless)

shapiro_raw <- shapiro.test(residuals(model_BD))
shapiro_log <- shapiro.test(residuals(model_BD_log))
cat("Shapiro-Wilk, raw-scale residuals:  p =", signif(shapiro_raw$p.value, 3), "\n")
cat("Shapiro-Wilk, log-scale residuals:  p =", signif(shapiro_log$p.value, 3), "\n")

use_log  <- shapiro_log$p.value > shapiro_raw$p.value
final_BD <- if (use_log) model_BD_log else model_BD
cat("Using", if (use_log) "LOG-transformed" else "RAW-scale", "model for inference.\n\n")

summary(final_BD)

## ---- Emmeans and CLD ----------------------------------------------------
emm_BD       <- emmeans(final_BD, "Type")
emm_BD_table <- summary(emm_BD, type = "response")
cld_BD       <- cld(emm_BD, Letters = letters, adjust = "sidak")

print(emm_BD_table)
print(cld_BD)

mean_col <- if ("response" %in% names(emm_BD_table)) "response" else "emmean"

BD_summary <- as.data.frame(emm_BD_table) %>%
  dplyr::select(Type, Mean = !!mean_col) %>%
  dplyr::left_join(
    as.data.frame(cld_BD) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  )

print(BD_summary)

## ---- Plots --------------------------------------------------------------
label_positions <- bulk_density_driftless %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(y_pos = max(Bulk_Density, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(BD_summary %>% dplyr::select(Type, CLD), by = "Type")

p_overall <- ggplot(bulk_density_driftless, aes(x = Type, y = Bulk_Density, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 6, fontface = "bold", color = "black"
  ) +
  scale_x_discrete(labels = type_labels) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  labs(
    title = "Bulk Density by Management Type",
    x = NULL,
    y = expression(paste("Bulk Density (g cm"^-3, ")"))
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = 50)
  )

print(p_overall)
ggsave("BD_overall_boxplot.png", p_overall, width = 10, height = 7, dpi = 300)

label_positions_county <- bulk_density_driftless %>%
  dplyr::group_by(County, Type) %>%
  dplyr::summarise(y_pos = max(Bulk_Density, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(BD_summary %>% dplyr::select(Type, CLD), by = "Type")

p_county <- ggplot(bulk_density_driftless, aes(x = Type, y = Bulk_Density, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions_county,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ County) +
  scale_fill_manual(values = fill_colors, labels = type_labels, name = "Management") +
  labs(
    title = "Bulk Density by Management Type and County",
    x = NULL,
    y = expression(paste("Bulk Density (g cm"^-3, ")"))
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text   = element_text(size = 13, face = "bold"),
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)
  )

print(p_county)
ggsave("BD_by_county_boxplot.png", p_county, width = 12, height = 7, dpi = 300)

cat("\nSaved BD_overall_boxplot.png and BD_by_county_boxplot.png\n")

# ---- Metadata export ----------------------------------------------------
anova_BD <- anova(final_BD)

BD_meta <- as.data.frame(emm_BD_table) %>%
  dplyr::rename(Mean = if ("response" %in% names(emm_BD_table)) "response" else "emmean") %>%
  dplyr::left_join(
    as.data.frame(cld_BD) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  ) %>%
  mutate(
    Metric         = "Bulk_Density",
    Scale          = if (use_log) "log" else "raw",
    Shapiro_raw_p  = signif(shapiro_raw$p.value, 3),
    Shapiro_log_p  = signif(shapiro_log$p.value, 3),
    LMM_F          = round(anova_BD$`F value`[1], 3),
    LMM_p          = signif(anova_BD$`Pr(>F)`[1], 3)
  ) %>%
  dplyr::select(Metric, Scale, Shapiro_raw_p, Shapiro_log_p,
                LMM_F, LMM_p, Type, Mean, CLD)

write.csv(BD_meta, "BD_metadata.csv", row.names = FALSE)
cat("Saved BD_metadata.csv\n")
