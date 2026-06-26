# DRACE 2024 — POXC (Permanganate Oxidizable Carbon)
# Cleaned analysis pipeline: import -> clean -> mixed model -> plots

library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)

# PART 1 — Import and clean
drpoxc24 <- read_xlsx(
  "C:\\Users\\swils\\Downloads\\01 Projects\\2024 CRP\\Sarah_Driftless_POXC\\Sarah_POXC_all_runs.xlsx",
  sheet = "Runs_and_avg"
)
colnames(drpoxc24) <- make.names(colnames(drpoxc24), unique = TRUE)
colnames(drpoxc24)
# Drop plate-consistency QC rows (Rep == "test")
drpoxc24 <- drpoxc24 %>% filter(Rep != "test")

redo_tubes <- drpoxc24 %>% filter(grepl("redo", Rep)) %>% pull(Tube.Label)
drpoxc24 <- drpoxc24 %>%
  filter(!(Tube.Label %in% redo_tubes & !grepl("redo", Rep)))

n_negative <- sum(drpoxc24$Avg_runs < 0)
cat("Remaining negative Avg_runs after redo fix:", n_negative, "\n")
drpoxc24 <- drpoxc24 %>% filter(Avg_runs >= 0)

#factor levels
drpoxc24$Type <- factor(
  drpoxc24$Type,
  levels = c("Conventional", "Soil_Health", "Rotational_Grazing", "New_CRP", "Old_CRP")
)

drpoxc24 <- drpoxc24 %>%
  mutate(Field = paste(Farmer, Type, sep = "_"))

cat("Rows after cleaning:", nrow(drpoxc24), "\n")
cat("Distinct fields:", n_distinct(drpoxc24$Field), "(expect 31)\n")
cat("Distinct sample locations (UID):", n_distinct(drpoxc24$UID), "(expect 154)\n\n")
stopifnot(n_distinct(drpoxc24$Field) == 31)

# PART 2 — Mixed model
model_POXC     <- lmer(Avg_runs ~ Type + (1 | Field), data = drpoxc24)
model_POXC_log <- lmer(log(Avg_runs) ~ Type + (1 | Field), data = drpoxc24)

shapiro_raw <- shapiro.test(residuals(model_POXC))
shapiro_log <- shapiro.test(residuals(model_POXC_log))
cat("Shapiro-Wilk, raw-scale residuals:  p =", signif(shapiro_raw$p.value, 3), "\n")
cat("Shapiro-Wilk, log-scale residuals:  p =", signif(shapiro_log$p.value, 3), "\n")

use_log <- shapiro_log$p.value > shapiro_raw$p.value
final_model <- if (use_log) model_POXC_log else model_POXC
cat("Using", if (use_log) "LOG-transformed" else "RAW-scale", "model for inference.\n\n")

summary(final_model)
#emmeans
emm_POXC       <- emmeans(final_model, "Type")
emm_POXC_table <- summary(emm_POXC, type = "response")
cld_POXC       <- cld(emm_POXC, Letters = letters, adjust = "sidak")

print(emm_POXC_table)
print(cld_POXC)

mean_col <- if ("response" %in% names(emm_POXC_table)) "response" else "emmean"

POXC_summary <- as.data.frame(emm_POXC_table) %>%
  dplyr::select(Type, Mean = !!mean_col) %>%
  dplyr::left_join(
    as.data.frame(cld_POXC) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  )

print(POXC_summary)

# PART 3 — Plotting
type_labels <- c(
  "Conventional"        = "Conventional Cropping & Tillage",
  "Soil_Health"         = "No Till and/or Cover Crops",
  "Rotational_Grazing"  = "Rotational Grazing",
  "New_CRP"             = "0-10 years CRP",
  "Old_CRP"             = "10+ years CRP"
)

fill_colors <- c(
  "Conventional"        = "#E69F00",
  "Soil_Health"         = "#56B4E9",
  "Rotational_Grazing"  = "#009E73",
  "New_CRP"             = "#F0E442",
  "Old_CRP"             = "#0072B2"
)

label_positions <- drpoxc24 %>%
  group_by(Type) %>%
  summarise(y_pos = max(Avg_runs, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(POXC_summary %>% dplyr::select(Type, CLD), by = "Type")

p_overall <- ggplot(drpoxc24, aes(x = Type, y = Avg_runs, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 6, fontface = "bold", color = "black"
  ) +
  scale_x_discrete(labels = type_labels) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  labs(
    title = "Permanganate Oxidizable Carbon by Management",
    x = NULL,
    y = "moles permanganate reduced per kg soil"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50)
  )
print(p_overall)
ggsave("POXC_overall_boxplot.png", p_overall, width = 10, height = 7, dpi = 300)

label_positions_county <- drpoxc24 %>%
  group_by(County, Type) %>%
  summarise(y_pos = max(Avg_runs, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(POXC_summary %>% dplyr::select(Type, CLD), by = "Type")

p_county <- ggplot(drpoxc24, aes(x = Type, y = Avg_runs, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions_county,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ County) +
  scale_fill_manual(values = fill_colors, labels = type_labels, name = "Management") +
  labs(
    title = "Permanganate Oxidizable Carbon by Management and County",
    x = NULL,
    y = "moles permanganate reduced per kg soil"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  )
print(p_county)
ggsave("POXC_by_county_boxplot.png", p_county, width = 12, height = 7, dpi = 300)

cat("\nSaved POXC_overall_boxplot.png and POXC_by_county_boxplot.png\n")

# ---- Metadata export ----------------------------------------------------
anova_POXC <- anova(final_model)

POXC_meta <- as.data.frame(emm_POXC_table) %>%
  dplyr::rename(Mean = if ("response" %in% names(emm_POXC_table)) "response" else "emmean") %>%
  dplyr::left_join(
    as.data.frame(cld_POXC) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  ) %>%
  mutate(
    Metric         = "POXC",
    Scale          = if (use_log) "log" else "raw",
    Shapiro_raw_p  = signif(shapiro_raw$p.value, 3),
    Shapiro_log_p  = signif(shapiro_log$p.value, 3),
    LMM_F          = round(anova_POXC$`F value`[1], 3),
    LMM_p          = signif(anova_POXC$`Pr(>F)`[1], 3)
  ) %>%
  dplyr::select(Metric, Scale, Shapiro_raw_p, Shapiro_log_p,
                LMM_F, LMM_p, Type, Mean, CLD)

write.csv(POXC_meta, "POXC_metadata.csv", row.names = FALSE)
cat("Saved POXC_metadata.csv\n")