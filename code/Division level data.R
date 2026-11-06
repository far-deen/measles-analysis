


library(readxl); library(tidyverse); library(lubridate)
library(incidence2); library(MASS); library(fitdistrplus)
library(patchwork); library(scales)




library(readxl); library(tidyverse); library(lubridate)
library(scales);  library(patchwork); library(spdep)
library(sf);      library(geodata);   library(ggrepel)
library(cowplot); library(viridis);   library(RColorBrewer)




# Load the Long_Format sheet


dwr = read_excel("Division level data.xlsx", sheet = "Long_Format")


dwr

head(dwr)

colnames(dwr)


dwr$date <- as.Date(dwr$date)

str(dwr$date)



summary_tbl = read_excel("Division level data.xlsx",
                         sheet = "Summary_Table")




epi <- dwr %>%
  group_by(date) %>%
  summarise(new_cases = sum(suspected_new_cases_24h, na.rm=TRUE))




ggplot(epi, aes(date, new_cases)) +
  geom_col(fill="#2E75B6", width=0.8) +
  geom_smooth(se=FALSE, colour="firebrick", linewidth=0.8) +
  scale_x_date(date_labels="%b %d", 
               breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),) +
  labs(title="National Epidemic Curve",
       subtitle="Daily new suspected cases, Apr 15 – May 15 2026",
       x="Date", y="New suspected cases") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("epicurve_national.png", width=10, height=5, dpi=150)





ggplot(dwr, aes(date, suspected_new_cases_24h, fill=division)) +
  geom_col(show.legend=FALSE) +
  facet_wrap(~division, scales="free_y", ncol=4) +
  scale_x_date(date_labels="%d %b", date_breaks="10 days") +
  labs(title="Division-wise Epidemic Curves",
       x="Date", y="New suspected cases/day") +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("epicurve_divisions.png", width=14, height=7, dpi=150)


peak_days <- dwr %>%
  group_by(division) %>%
  slice_max(order_by = suspected_new_cases_24h, n = 1) %>%
  arrange(date)

peak_days

peak_days %>%
  dplyr::select(division, date, suspected_new_cases_24h)




ggplot(peak_days, aes(x = date, y = reorder(division, date))) +
  geom_point(size = 4, color = "#2E75B6") +
  geom_line(aes(group = 1), color = "firebrick", linewidth = 1) +
  scale_x_date(date_labels="%b %d", 
               breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),)+
  labs(title = "Division ‑ Level Peak Timing (Apr 15–May 15, 2026)",
       x = "Peak Date", y = "Division") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



ggplot(summary_tbl,
       aes(reorder(division, attack_rate_per_100k),
           attack_rate_per_100k, fill=attack_rate_per_100k)) +
  geom_col(show.legend=FALSE) +
  geom_text(aes(label=round(attack_rate_per_100k,1)),
            hjust=-0.1, size=3.5) +
  coord_flip() +
  scale_fill_gradient(low="#AED6F1", high="#1A5276") +
  labs(title="Attack Rate by Division",
       x="", y="Suspected cases per 100,000 population") +
  theme_minimal()
ggsave("attack_rates.png", width=8, height=5, dpi=150)


























# Quick descriptive check



dwr$date     <- as.Date(dwr$date)
dwr$division <- as.factor(dwr$division)


x <- as.integer(na.omit(dwr$suspected_new_cases_24h))
cat("\nDescriptive stats on daily new cases:\n")
cat("Mean:", round(mean(x),2), "\n")
cat("Variance:", round(var(x),2), "\n")
cat("Dispersion (var/mean):", round(var(x)/mean(x),2),
    "→ if >1: overdispersed, NB preferred over Poisson\n")



# Barishal is reference (first alphabetically) — confirm:

cat("Reference level:", levels(dwr$division)[1], "\n")






# Negative binomial

neg_binom_model = glm.nb(suspected_new_cases_24h ~ division, data = dwr)


summary(neg_binom_model)



theta_hat <- neg_binom_model$theta
theta_se  <- neg_binom_model$SE.theta
cat("Theta:", round(theta_hat, 4),
    "| SE:", round(theta_se, 4), "\n")
cat("Variance inflation vs Poisson (at mean):",
    round(1 + mean(x)/theta_hat, 3), "x\n")







# Poison model 

poisson_model <- glm(
  suspected_new_cases_24h ~ division,
  family = poisson(link = "log"),
  data   = dwr
)
summary(poisson_model)





# Over dispersion test

pois_disp <- deviance(poisson_model) / df.residual(poisson_model)

cat("\nPoisson residual deviance/df:", round(pois_disp, 2),
    "→ If >2: severe overdispersion\n")












long <- read_excel("Division level data.xlsx", sheet = "Long_Format")
pop  <- read_excel("Division level data.xlsx", sheet = "Population")



dwr <- long %>%
  left_join(pop, by = "division")

dwr <- dwr %>%
  mutate(
    failures = population - suspected_new_cases_24h
  )


binomial_model <- glm(
  cbind(suspected_new_cases_24h, failures) ~ division,
  family = binomial(link = "logit"),
  data = dwr
)

summary(binomial_model)









AIC = AIC(poisson_model,neg_binom_model,binomial_model)

write.csv(AIC, "AIC.csv", row.names = FALSE)






pvals <- list(
  Poisson      = coef(summary(poisson_model))[,4],
  NegBinomial  = coef(summary(neg_binom_model))[,4],
  Binomial     = coef(summary(binomial_model))[,4]
)

pvals


write.csv(pvals, "pvals.csv", row.names = FALSE)











# LIKELIHOOD RATIO TEST — Poisson vs NB
# H0: data follows Poisson (theta = Inf)


lrt_stat <- 2 * (logLik(neg_binom_model) - logLik(poisson_model))
lrt_p    <- pchisq(as.numeric(lrt_stat), df = 1, lower.tail = FALSE)/2
cat("\nLRT statistic:", round(lrt_stat, 2))
cat("\nLRT p-value:", format(lrt_p, scientific = TRUE))
cat("\nConclusion:", ifelse(lrt_p < 0.05,
                            "Reject Poisson — NB is significantly better",
                            "Cannot reject Poisson"), "\n")







# INCIDENCE RATE RATIOS (IRR) TABLE — for your paper


irr_table <- tidy(neg_binom_model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(
    term = gsub("division", "", term),
    IRR = round(estimate, 3),
    CI_low = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p_value = round(p.value, 4),
    signif = case_when(
      p.value < 0.001 ~ "****",
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      TRUE ~ "ns"
    )
  ) %>%
  rename(Division = term) %>%
  dplyr::select(Division, IRR, CI_low, CI_high, p_value, signif)



cat("\n=== INCIDENCE RATE RATIOS (for paper Table) ===\n")
print(irr_table)
write.csv(irr_table, "irr_table_NB.csv", row.names = FALSE)










# COMPARISON TABLE


comparison <- tibble(
  Model = c("Poisson", "Negative Binomial"),
  AIC   = c(AIC(poisson_model), AIC(neg_binom_model)),
  BIC   = c(BIC(poisson_model), BIC(neg_binom_model)),
  Resid_Deviance = c(
    round(deviance(poisson_model), 2),
    round(deviance(neg_binom_model), 2)
  ),
  Resid_df   = c(
    df.residual(poisson_model),
    df.residual(neg_binom_model)
  ),
  Dev_per_df = round(c(
    deviance(poisson_model) / df.residual(poisson_model),
    deviance(neg_binom_model) / df.residual(neg_binom_model)
  ), 3),
  Log_lik = c(
    round(as.numeric(logLik(poisson_model)), 2),
    round(as.numeric(logLik(neg_binom_model)), 2)
  ),
  Notes = c(
    "Overdispersion ignored — WRONG for this data",
    "Best fit — overdispersion correctly modelled"
  )
)


print(comparison)

write.csv(comparison, "model_comparison_table.csv", row.names = FALSE)











































































# Calculate nationwide doubling time

nat_cumul = dwr %>%
  group_by(date) %>%
  summarise(cumul = sum(suspected_cases_cumulative, na.rm=TRUE)/8) %>%
  # (divide by 8 to avoid double-counting: use max per date instead)
  ungroup()

early = nat_cumul %>% slice_head(n=10) %>%
  mutate(t = row_number(), log_cumul = log(cumul))

model = lm(log_cumul ~ t, data=early)


doubling_time = log(2) / coef(model)["t"]


cat("Estimated doubling time:", round(doubling_time, 1), "days")
cat("R-squared:", round(summary(model)$r.squared, 3))





ggplot(doubling_summary, aes(x = reorder(division, doubling_time), 
                             y = doubling_time, fill = division)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Estimated Doubling Time by Division",
    x = "Division",
    y = "Doubling Time (days)"
  ) +
  theme_minimal()









# Calculate division-wise doubling time

doubling_summary <- dwr %>%
  group_by(division, date) %>%
  summarise(cumul = sum(suspected_cases_cumulative, na.rm = TRUE), .groups = "drop_last") %>%
  group_modify(~ {
    early <- .x %>%
      slice_head(n = 10) %>%
      mutate(t = row_number(), log_cumul = log(cumul))
    model <- lm(log_cumul ~ t, data = early)
    doubling_time <- log(2) / coef(model)["t"]
    tibble(
      doubling_time = round(doubling_time, 1),
      r_squared = round(summary(model)$r.squared, 3)
    )
  }) %>%
  ungroup()

# View summary

print(doubling_summary)



# Visualize division-wise doubling time

ggplot(doubling_summary, aes(x = reorder(division, doubling_time), y = doubling_time, fill = division)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = paste0(doubling_time, " days")), hjust = -0.1, size = 3.5) +
  labs(
    title = "Estimated Doubling Time by Division",
    x = "Division",
    y = "Doubling Time (days)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))













# Primary Report (CFR)


# Convert date column to Date class after summarizing

cfr_daily <- dwr %>%
  group_by(date, division) %>%
  summarise(
    cfr_susp = suspected_deaths_cumulative / suspected_cases_cumulative * 100,
    .groups = "drop"
  )

# Ensure date is Date type

cfr_daily$date <- as.Date(cfr_daily$date, format = "%Y-%m-%d")

# Plot again

ggplot(cfr_daily, aes(x = date, y = cfr_susp, colour = division)) +
  geom_line(linewidth = 0.7) +
  scale_x_date(date_labels="%b %d", 
               breaks = seq(min(rdm$date) , max(rdm$date), by = "7 day"),) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "CFR Trend by Division (Suspected)",
    x = "Date",
    y = "Case Fatality Rate (%)",
    colour = "Division"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# Secondary Report (CFR) 


# Convert date column to Date class after summarizing

cfr_daily <- dwr %>%
  group_by(date, division) %>%
  summarise(
    cfr_conf = confirmed_deaths_cumulative / confirmed_cases_cumulative * 100,
    .groups = "drop"
  )

# Ensure date is Date type

cfr_daily$date <- as.Date(cfr_daily$date, format = "%Y-%m-%d")

# Plot again

ggplot(cfr_daily, aes(x = date, y = cfr_conf, colour = division)) +
  geom_line(linewidth = 0.7) +
  scale_x_date(date_labels="%b %d", 
               breaks = seq(min(rdm$date) , max(rdm$date), by = "7 day"),) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "CFR Trend by Division (Confirm)",
    x = "Date",
    y = "Case Fatality Rate (%)",
    colour = "Division"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))





# Summarize cumulative suspected cases by division and date

suspected_division <- dwr %>%
  group_by(date, division) %>%
  summarise(suspected_cases_cumulative = sum(suspected_cases_cumulative, na.rm = TRUE),
            .groups = "drop")

# Plot division-wise cumulative suspected cases

ggplot(suspected_division, aes(x = date, y = suspected_cases_cumulative, fill = division)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  labs(
    title = "Suspected Cases - Cumulative by Division",
    x = "Date",
    y = "Number of Suspected Cases (Cumulative)",
    fill = "Division"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Summarize your dataset by division


division_summary <- dwr %>%
  group_by(division) %>%
  summarise(reported_cases = sum(confirmed_cases_cumulative, na.rm = TRUE)) %>%
  mutate(percentage = round((reported_cases / sum(reported_cases)) * 100, 2))



# Plot division-wise reported cases (Confirmed)


ggplot(division_summary, aes(x = reorder(division, reported_cases),
                             y = reported_cases, fill = division)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = paste0(reported_cases, " (", percentage, "%)")),
            hjust = -0.2, color = "black", size = 3.5) +
  labs(
    title = "Division-wise Reported Cases",
    x = "Division",
    y = "Reported Cases"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))






# Log-linear regression

doubling_summary <- dwr %>%
  group_by(division, date) %>%
  summarise(cumul = sum(suspected_cases_cumulative, na.rm = TRUE), .groups = "drop_last") %>%
  group_modify(~ {
    df <- .x %>%
      slice_head(n = 10) %>%        # first 10 days
      mutate(t = row_number(),
             log_cumul = log(cumul))
    
    model <- lm(log_cumul ~ t, data = df)
    
    tibble(
      growth_rate = coef(model)["t"],
      doubling_time = round(log(2) / coef(model)["t"], 1),
      r_squared = round(summary(model)$r.squared, 3)
    )
  }) %>%
  ungroup()



# Visualizng the growth rate

ggplot(doubling_summary, aes(x = reorder(division, growth_rate), y = growth_rate, fill = r_squared)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = paste0(round(growth_rate, 3))), hjust = -0.1, size = 3.5) +
  scale_fill_gradient(low = "orange", high = "darkgreen") +
  labs(
    title = "Division-wise Growth Rate (Colored by R²)",
    x = "Division",
    y = "Growth Rate",
    fill = "R²"
  ) +
  theme_minimal()




# Calculate total cases per division

division_cases <- dwr %>%
  group_by(division) %>%
  summarise(total_cases = sum(suspected_cases_cumulative, na.rm = TRUE)) %>%
  mutate(share = total_cases / sum(total_cases),
         hhi_component = share^2)

# Compute HHI

hhi_value <- sum(division_cases$hhi_component)

# View results

division_cases
cat("Herfindahl-Hirschman Index (HHI):", round(hhi_value, 3))


ggplot(division_cases, aes(x = reorder(division, share), y = share, fill = division)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = paste0(round(share * 100, 1), "%")), hjust = -0.1, size = 3.5) +
  labs(
    title = paste0("Division-wise Case Share (HHI = ", round(hhi_value, 3), ")"),
    x = "Division",
    y = "Share of Total Cases"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0.05, 0.25)))











# HHI trend over time


hhi_trend <- dwr %>%
  group_by(date, division) %>%
  summarise(cases = sum(suspected_cases_cumulative, na.rm = TRUE), .groups = "drop_last") %>%
  group_by(date) %>%
  mutate(share = cases / sum(cases),
         hhi_component = share^2) %>%
  summarise(HHI = sum(hhi_component))

ggplot(hhi_trend, aes(x = date, y = HHI)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(title = "HHI Trend Over Time", x = "Date", y = "Herfindahl-Hirschman Index") +
  theme_minimal()+
  theme(
    plot.title = element_text(
      hjust = 0.5,          # 0 = left, 0.5 = center, 1 = right
      # optional: make it bold
      size = 14             # optional: adjust font size
    ))









division_totals <- Long_Format %>%
  filter(date >= ymd("2026-04-15") & date <= ymd("2026-05-15")) %>%
  group_by(division) %>%
  summarise(
    total_confirmed_cases  = sum(confirmed_new_cases_24h, na.rm = TRUE),
    total_confirmed_deaths = sum(confirmed_new_deaths_24h, na.rm = TRUE)
  ) %>%
  arrange(desc(total_confirmed_cases))

division_totals




long <- read_excel("measles_outbreak_organised.xlsx",
                   sheet = "Long_Format")
dwr$date <- as.Date(dwr$date)

summary_tbl <- read_excel("measles_outbreak_organised.xlsx",
                          sheet = "Summary_Table")

# Population lookup (2022 estimates)
pop_lookup <- c(
  Barishal   = 8325666,  Chattogram = 28423019,
  Dhaka      = 36054418, Mymensingh = 11370000,
  Khulna     = 15563000, Rajshahi   = 18484858,
  Rangpur    = 15787000, Sylhet     = 9910219
)











# CHOROPLETH MAPS 



bgd_shp = geodata::gadm(country = "BGD", level = 1,
                        path = tempdir(), version = "4.1")
bgd_sf  = sf::st_as_sf(bgd_shp)



# Name harmonisation lookup

name_map = c(
  "Barisal"    = "Barishal",
  "Chittagong" = "Chattogram",
  "Dhaka"      = "Dhaka",
  "Khulna"     = "Khulna",
  "Mymensingh" = "Mymensingh",
  "Rajshahi"   = "Rajshahi",
  "Rangpur"    = "Rangpur",
  "Sylhet"     = "Sylhet"
)

bgd_sf = bgd_sf %>%
  mutate(division = recode(NAME_1, !!!name_map))




start_date <- as.Date("2026-04-15")
end_date   <- as.Date("2026-05-15")




map_data <- dwr %>%
  filter(date >= start_date & date <= end_date) %>%
  group_by(division) %>%
  summarise(
    cases  = sum(confirmed_new_cases_24h, na.rm = TRUE),
    deaths = sum(confirmed_new_deaths_24h, na.rm = TRUE)
  )


# Merge with shapefile

bgd_map <- bgd_sf %>%
  left_join(map_data, by = "division")


centroids = bgd_map %>%
  mutate(
    cx = map_dbl(geometry, ~st_centroid(.) %>% st_coordinates() %>% .[1]),
    cy = map_dbl(geometry, ~st_centroid(.) %>% st_coordinates() %>% .[2])
  ) %>%
  st_drop_geometry()


# Map theme 

map_theme = theme_void(base_family = "sans") +
  theme(
    plot.title       = element_text(size = 11, face = "bold",
                                    hjust = 0.5, margin = margin(b = 8)),
    legend.position  = "right",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height= unit(1.2, "cm"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 5, 10, 5)
  )


# MAP 1 — Cases

p_cases = ggplot(bgd_map) +
  geom_sf(aes(fill = cases), colour = "black", linewidth = 0.4) +
  scale_fill_gradient(
    low      = "#FDDBC7",   # very light salmon/pink
    high     = "#67001F",   # very dark red (matches image)
    name     = "Measles\ncases",
    labels   = comma,
    guide    = guide_colorbar(barwidth = 0.8, barheight = 8,
                              title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(data = centroids,
            aes(x = cx, y = cy,
                label = paste0(division, "\n", comma(cases))),
            size = 2.5, fontface = "bold",
            colour = ifelse(centroids$cases > 4000, "white", "black"),
            lineheight = 0.85) +
  labs(title = "Division-wise Measles Cases in Bangladesh") +
  map_theme


# MAP 2 — Deaths

p_deaths = ggplot(bgd_map) +
  geom_sf(aes(fill = deaths), colour = "black", linewidth = 0.4) +
  scale_fill_gradient(
    low      = "#C6DBEF",   # very light blue
    high     = "#084594",   # dark navy blue (matches image)
    name     = "Measles\ndeaths",
    labels   = comma,
    guide    = guide_colorbar(barwidth = 0.8, barheight = 8,
                              title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(data = centroids,
            aes(x = cx, y = cy,
                label = paste0(division, "\n", deaths)),
            size = 2.5, fontface = "bold",
            colour = ifelse(centroids$deaths > 40, "white", "black"),
            lineheight = 0.85) +
  labs(title = "Division-wise Measles Deaths in Bangladesh") +
  map_theme

# Combine side by side 

maps_combined = plot_grid(p_cases, p_deaths,
                          ncol = 2, align = "hv",
                          rel_widths = c(1, 1))

















# ── 2A. Daily active hospital occupancy per division ─────────────────────────
# Active beds = cumulative admitted - cumulative discharged

dwr <- dwr %>%
  mutate(active_beds = suspected_admitted_cumulative -
           suspected_discharged_cumulative)

# ── 2B. Stacked area chart — occupancy by division over time ─────────────────

occ_palette <- c(
  Dhaka      = "#1A5276", Rajshahi   = "#C0392B",
  Chattogram = "#27AE60", Barishal   = "#E67E22",
  Khulna     = "#8E44AD", Sylhet     = "#2980B9",
  Mymensingh = "#D35400", Rangpur    = "#7F8C8D"
)



p_hosp_stack <- dwr %>%
  ggplot(aes(date, active_beds, fill = division)) +
  geom_area(position = "stack", alpha = 0.85, colour = "white",
            linewidth = 0.2) +
  scale_fill_manual(values = occ_palette, name = "Division") +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Hospital Occupancy Over Time — All Divisions",
    subtitle = "Active beds = Cumulative admissions − Cumulative discharges",
    x        = "Date", y = "Estimated active hospital beds"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))










# ── 2C. Daily new admissions vs discharges — national level ──────────────────


nat_hosp <- dwr %>%
  group_by(date) %>%
  summarise(
    new_admissions  = sum(suspected_admitted_24h,    na.rm = TRUE),
    new_discharges  = sum(suspected_discharged_24h,  na.rm = TRUE),
    net_change      = new_admissions - new_discharges
  )



p_hosp_flow <- nat_hosp %>%
  pivot_longer(cols = c(new_admissions, new_discharges),
               names_to = "type", values_to = "count") %>%
  mutate(type = recode(type,
                       new_admissions  = "New admissions",
                       new_discharges  = "New discharges")) %>%
  ggplot(aes(date, count, colour = type)) +
  geom_line(linewidth = 1.2) +
  geom_area(aes(fill = type), alpha = 0.12, position = "identity") +
  scale_colour_manual(values = c("New admissions"  = "#C0392B",
                                 "New discharges"  = "#27AE60"),
                      name = "") +
  scale_fill_manual(values  = c("New admissions"  = "#C0392B",
                                "New discharges" = "#27AE60"),
                    guide = "none") +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  labs(
    title    = "Daily Hospital Admissions vs Discharges — National",
    subtitle = "Gap between lines = net increase in hospital burden",
    x = "Date", y = "Patients"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))








# ── 2D. Division-level peak occupancy bar ────────────────────────────────────


peak_occ <- dwr %>%
  group_by(division) %>%
  summarise(peak_occupancy = max(active_beds, na.rm = TRUE),
            final_admitted  = last(suspected_admitted_cumulative),
            final_discharge = last(suspected_discharged_cumulative),
            final_occupancy = last(active_beds))

p_peak_occ <- peak_occ %>%
  mutate(division = reorder(division, peak_occupancy)) %>%
  ggplot(aes(division, peak_occupancy,
             fill = reorder(division, peak_occupancy))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(peak_occupancy)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = occ_palette) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Peak Hospital Occupancy by Division",
    x     = "", y = "Maximum active beds (simultaneously)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))










# ── 2E. Discharge rate (health system efficiency) ────────────────────────────


disch_rate <- dwr %>%
  group_by(date) %>%
  summarise(
    total_admitted  = sum(suspected_admitted_cumulative),
    total_discharged= sum(suspected_discharged_cumulative)
  ) %>%
  mutate(discharge_rate_pct = total_discharged / total_admitted * 100)

p_discharge_rate <- ggplot(disch_rate, aes(date, discharge_rate_pct)) +
  geom_line(colour = "#1A5276", linewidth = 1.2) +
  geom_hline(yintercept = 80, linetype = "dashed",
             colour = "orange", alpha = 0.8) +
  annotate("text", x = min(disch_rate$date) + 3,
           y = 81, label = "80% discharge threshold",
           colour = "orange", size = 3) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  scale_y_continuous(labels = percent_format(scale = 1),
                     limits = c(50, 100)) +
  labs(
    title    = "Cumulative Discharge Rate Over Time",
    subtitle = "Discharged / Admitted × 100 — proxy for system clearance capacity",
    x = "Date", y = "Discharge rate (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Combine hospital panels

hosp_combined <- (p_hosp_flow / p_hosp_stack) |
  (p_peak_occ / p_discharge_rate)
ggsave("hospital_burden_analysis.png", hosp_combined,
       width = 16, height = 12, dpi = 150, bg = "white")
cat("✓ Saved: hospital_burden_analysis.png\n")












# ── 3A. Daily confirmation rate per division ──────────────────────────────────

conf_rate <- dwr %>%
  mutate(
    conf_rate_pct = confirmed_cases_cumulative /
      suspected_cases_cumulative * 100
  )

# ── 3B. Line chart: confirmation rate over time, faceted by division ──────────
p_conf_facet <- conf_rate %>%
  ggplot(aes(date, conf_rate_pct, colour = division)) +
  geom_line(linewidth = 0.9, show.legend = FALSE) +
  geom_point(size = 0.8, show.legend = FALSE) +
  facet_wrap(~division, ncol = 4, scales = "free_y") +
  scale_x_date(date_labels = "%d %b", date_breaks = "10 days") +
  scale_y_continuous(labels = function(x) paste0(round(x, 0), "%")) +
  labs(
    title    = "Confirmation Rate Trend by Division",
    subtitle = "Confirmed cumulative / Suspected cumulative × 100",
    x = "Date", y = "Confirmation rate (%)"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text   = element_text(face = "bold"))




# ── 3C. National confirmation rate vs daily new suspected ─────────────────────


nat_conf <- dwr %>%
  group_by(date) %>%
  summarise(
    nat_suspected  = sum(suspected_cases_cumulative),
    nat_confirmed  = sum(confirmed_cases_cumulative),
    nat_new_susp   = sum(suspected_new_cases_24h, na.rm = TRUE),
    nat_new_conf   = sum(confirmed_new_cases_24h, na.rm = TRUE)
  ) %>%
  mutate(
    cumul_conf_rate = nat_confirmed / nat_suspected * 100,
    daily_conf_rate = ifelse(nat_new_susp > 0,
                             nat_new_conf / nat_new_susp * 100, NA)
  )

# Dual-axis: new suspected (bars) + confirmation rate (line)


p_conf_national <- ggplot(nat_conf, aes(x = date)) +
  geom_col(aes(y = nat_new_susp),
           fill = "#AED6F1", alpha = 0.8, width = 0.8) +
  geom_line(aes(y = cumul_conf_rate * 80),  # scaled for dual axis
            colour = "#C0392B", linewidth = 1.2) +
  geom_point(aes(y = cumul_conf_rate * 80),
             colour = "#C0392B", size = 1.5) +
  scale_y_continuous(
    name     = "Daily new suspected cases",
    labels   = comma,
    sec.axis = sec_axis(~./80,
                        name   = "Cumulative confirmation rate (%)",
                        labels = function(x) paste0(round(x, 1), "%"))
  ) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  labs(
    title    = "Confirmation Rate vs Epidemic Activity — National",
    subtitle = "Bars = daily new suspected | Red line = cumulative confirmation rate",
    x = "Date"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(colour = "#C0392B"),
    axis.text.y.right  = element_text(colour = "#C0392B"),
    plot.title    = element_text(hjust = 0.5,face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )









# ── Final confirmation rate summary table ─────────────────────────────────


conf_summary <- dwr %>%
  filter(date == snapshot_date) %>%
  mutate(
    total_suspected=division_totals$total_suspected_cases[match(division, division_totals$division)],
    total_confirmed  = division_totals$total_confirmed_cases[match(division, division_totals$division)],
    confirmation_rate_pct = round(total_confirmed /
                                    total_suspected * 100, 1)
  ) %>%
  select(division, total_suspected, total_confirmed,confirmation_rate_pct)


cat("\n=== Confirmation Rate by Division (May 15) ===\n")

print(conf_summary)




# Confirmation Rate by Division (April 15 - May 15)

p_conf_bar <- conf_summary %>%
  mutate(division = reorder(division, confirmation_rate_pct)) %>%
  ggplot(aes(division, confirmation_rate_pct,
             fill = confirmation_rate_pct)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(confirmation_rate_pct, "%")),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Confirmation Rate by Division (April 15 - May 15)",
    x = "", y = "Confirmation rate (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme( plot.title    = element_text(hjust = 0.5,face = "bold"))








#  Confirmed vs Suspected case trajectory


dwr$date <- as.Date(dwr$date)

nat_traj <- dwr %>%
  group_by(date) %>%
  summarise(
    suspected = sum(suspected_cases_cumulative),
    confirmed = sum(confirmed_cases_cumulative)
  ) %>%
  pivot_longer(-date, names_to = "type", values_to = "count")

p_traj <- ggplot(nat_traj, aes(date, count, colour = type)) +
  geom_line(linewidth = 1.3) +
  scale_colour_manual(values = c(suspected = "#AED6F1",
                                 confirmed  = "#1A5276"),
                      labels = c("Confirmed (lab)", "Suspected (clinical)"),
                      name   = "Case definition") +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Confirmed vs Suspected Cases — Cumulative Trajectory",
    x = "Date", y = "Cumulative cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))








# Confirmed CFR vs Suspected CFR trend


cfr_trend <- dwr %>%
  group_by(date) %>%
  summarise(
    cfr_suspected = sum(suspected_deaths_cumulative) /
      sum(suspected_cases_cumulative)  * 100,
    cfr_confirmed = sum(confirmed_deaths_cumulative) /
      sum(confirmed_cases_cumulative)  * 100
  ) %>%
  pivot_longer(-date, names_to = "type", values_to = "cfr")

p_cfr_trend <- ggplot(cfr_trend, aes(date, cfr, colour = type)) +
  geom_line(linewidth = 1.3) +
  scale_colour_manual(
    values = c(cfr_suspected = "#E67E22", cfr_confirmed = "#C0392B"),
    labels = c("CFR (confirmed)", "CFR (suspected)"),
    name   = ""
  ) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), "%")) +
  labs(
    title    = "Case Fatality Rate Trend — Confirmed vs Suspected",
    x = "Date", y = "CFR (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title    = element_text(hjust = 0.5,face = "bold"))









# Confirmed Cases: Actual vs Population-Expected


total_pop    <- sum(pop_lookup)
total_conf   = division_totals$total_confirmed_cases
names(total_conf) = division_totals$division

conf_vs_expected <- tibble(
  division      = names(pop_lookup),
  population    = pop_lookup,
  pop_share_pct = pop_lookup / total_pop * 100,
  confirmed     = as.numeric(total_conf[names(pop_lookup)])
) %>%
  mutate(
    total_confirmed    = sum(confirmed, na.rm = TRUE),
    expected_confirmed = pop_share_pct / 100 * total_confirmed,
    excess_pct         = (confirmed - expected_confirmed) /
      expected_confirmed * 100
  )

p_excess <- conf_vs_expected %>%
  mutate(
    division = reorder(division, excess_pct),
    direction = ifelse(excess_pct >= 0, "Above expected", "Below expected")
  ) %>%
  ggplot(aes(division, excess_pct, fill = direction)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_text(aes(label = paste0(ifelse(excess_pct > 0, "+", ""),
                               round(excess_pct, 0), "%")),
            hjust = ifelse(conf_vs_expected$excess_pct[
              order(conf_vs_expected$excess_pct)] > 0, -0.1, 1.1),
            size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("Above expected" = "#C0392B",
                               "Below expected" = "#2980B9"),
                    name = "") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Confirmed Cases: Actual vs Population-Expected",
    x = "", y = "Deviation from expected (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title    = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


































# Moran's I




bgd_shp <- gadm(country="BGD", level=1, path=tempdir(), version="4.1")
bgd_sf  <- st_as_sf(bgd_shp)



name_map <- c("Barisal"="Barishal","Chittagong"="Chattogram",
              "Dhaka"="Dhaka","Khulna"="Khulna",
              "Mymensingh"="Mymensingh","Rajshahi"="Rajshahi",
              "Rangpur"="Rangpur","Sylhet"="Sylhet")
bgd_sf <- bgd_sf %>%
  mutate(division = recode(NAME_1, !!!name_map))


summary_tbl <- read_excel("Division level data.xlsx", sheet="Summary_Table")


bgd_sf <- bgd_sf %>%
  left_join(
    summary_tbl %>%
      dplyr::select(
        division,
        attack_rate = attack_rate_per_100k,
        total_susp = total_suspected_cases,
        total_deaths = total_suspected_deaths,
        cfr = cfr_suspected_pct,
        hosp_occ = hospital_occupancy_May15
      ),
    by = "division"
  )


cents <- bgd_sf %>%
  mutate(cx=map_dbl(geometry,~st_centroid(.)%>%st_coordinates()%>%.[1]),
         cy=map_dbl(geometry,~st_centroid(.)%>%st_coordinates()%>%.[2])) %>%
  st_drop_geometry()


nb <- poly2nb(bgd_sf, queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=TRUE)


moran_g <- moran.test(bgd_sf$attack_rate, lw,
                      randomisation=TRUE, zero.policy=TRUE)
cat("\n=== Global Moran's I — Attack Rate ===\n"); print(moran_g)

# Extract I and p for annotation
moran_I <- round(moran_g$statistic["Moran I statistic"], 4)
moran_p <- round(moran_g$p.value, 4)
moran_interp <- ifelse(moran_I > 0 & moran_p < 0.05,
                       "Significant positive clustering",
                       ifelse(moran_I > 0, "Positive but not significant (n=8)", "No clustering"))
cat("Interpretation:", moran_interp, "\n")


moran_df <- data.frame(
  division   = bgd_sf$division,
  x          = scale(bgd_sf$attack_rate)[,1],
  lag_x      = lag.listw(lw, scale(bgd_sf$attack_rate)[,1]),
  attack_rate= bgd_sf$attack_rate
)
moran_df$quadrant <- case_when(
  moran_df$x > 0 & moran_df$lag_x > 0 ~ "HH (hot spot)",
  moran_df$x < 0 & moran_df$lag_x < 0 ~ "LL (cold spot)",
  moran_df$x > 0 & moran_df$lag_x < 0 ~ "HL (outlier)",
  moran_df$x < 0 & moran_df$lag_x > 0 ~ "LH (outlier)",
  TRUE ~ "Neutral"
)

p_moran_scatter <- ggplot(moran_df, aes(x, lag_x, colour=quadrant)) +
  geom_hline(yintercept=0, linetype="dashed", colour="gray60") +
  geom_vline(xintercept=0, linetype="dashed", colour="gray60") +
  geom_smooth(aes(group=1), method="lm", se=TRUE,
              colour="#1A5276", fill="#AED6F1", linewidth=1,
              linetype="solid") +
  geom_point(size=4, alpha=0.9) +
  geom_text_repel(aes(label=paste0(division,"\n",attack_rate,"/100k")),
                  size=3, fontface="bold", show.legend=FALSE) +
  scale_colour_manual(values=c("HH (hot spot)"="#C0392B",
                               "LL (cold spot)"="#2980B9",
                               "HL (outlier)"  ="#E67E22",
                               "LH (outlier)"  ="#F39C12"),
                      name="LISA quadrant") +
  annotate("text", x=-1.6, y=max(moran_df$lag_x)*0.95,
           label=paste0("Global Moran's I = ", moran_I,
                        "\np-value = ", moran_p,
                        "\n", moran_interp),
           size=3, hjust=0, colour="#1A5276",
           family="mono") +
  labs(title="Moran Scatter Plot — Measles Attack Rate by Division",
       subtitle="Each point = one division. Positive slope = spatial clustering.",
       x="Standardised attack rate (z-score)",
       y="Spatially lagged attack rate") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        legend.position="right")




lisa        <- localmoran(bgd_sf$attack_rate, lw, zero.policy=TRUE)
bgd_sf$local_I  <- lisa[,1]
bgd_sf$local_p  <- lisa[,5]
bgd_sf$z_score  <- scale(bgd_sf$attack_rate)[,1]
bgd_sf$lag_z    <- lag.listw(lw, bgd_sf$z_score)

bgd_sf$cluster <- case_when(
  bgd_sf$z_score >  0 & bgd_sf$lag_z >  0 & bgd_sf$local_p < 0.05 ~ "HH — Hot spot",
  bgd_sf$z_score <  0 & bgd_sf$lag_z <  0 & bgd_sf$local_p < 0.05 ~ "LL — Cold spot",
  bgd_sf$z_score >  0 & bgd_sf$lag_z <  0 & bgd_sf$local_p < 0.05 ~ "HL — Spatial outlier",
  bgd_sf$z_score <  0 & bgd_sf$lag_z >  0 & bgd_sf$local_p < 0.05 ~ "LH — Spatial outlier",
  TRUE ~ "Not significant"
)

map_theme <- theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=11),
        plot.subtitle=element_text(hjust=0.5,size=9,colour="gray40"),
        legend.position="right",
        legend.title=element_text(size=9,face="bold"),
        legend.text=element_text(size=8),
        plot.background=element_rect(fill="white",colour=NA))

p_lisa <- ggplot(bgd_sf) +
  geom_sf(aes(fill=cluster), colour="black", linewidth=0.5) +
  scale_fill_manual(
    values=c("HH — Hot spot"       = "#C0392B",
             "LL — Cold spot"      = "#2980B9",
             "HL — Spatial outlier"= "#E67E22",
             "LH — Spatial outlier"= "#F1C40F",
             "Not significant"     = "#ECECEC"),
    name="LISA cluster") +
  geom_text(data=cents,
            aes(x=cx, y=cy,
                label=paste0(division,"\n",attack_rate,"/100k")),
            size=2.6, fontface="bold",
            colour=ifelse(cents$attack_rate>45,"white","black"),
            lineheight=0.85) +
  labs(title="LISA Cluster Map — Measles Attack Rate",
       subtitle=paste0("Global Moran's I = ",moran_I," | p = ",moran_p)) +
  map_theme

# Attack rate choropleth

p_attack_map <- ggplot(bgd_sf) +
  geom_sf(aes(fill=attack_rate), colour="black", linewidth=0.4) +
  scale_fill_gradient2(low="#FDFEFE", mid="#F1948A", high="#7B241C",
                       midpoint=median(bgd_sf$attack_rate, na.rm=TRUE),
                       name="Cases\nper 100k",
                       guide=guide_colorbar(barwidth=0.8, barheight=7)) +
  geom_text(data=cents,
            aes(x=cx, y=cy,
                label=paste0(division,"\n",round(attack_rate,1))),
            size=2.6, fontface="bold",
            colour=ifelse(cents$attack_rate>50,"white","black"),
            lineheight=0.85) +
  labs(title="Attack Rate per 100,000 Population",
       subtitle="Suspected cumulative cases | May 15, 2026") +
  map_theme

moran_combined <- (p_attack_map | p_lisa) / p_moran_scatter +
  plot_layout(heights=c(2,1.5)) +
  plot_annotation(
    title="Spatial Analysis — Measles Outbreak Concentration",
    theme=theme(plot.title=element_text(face="bold",size=14))
  )


print(bgd_sf %>% st_drop_geometry() %>%
        dplyr:: select(division, attack_rate, local_I, local_p, cluster))



Local_Mor_I = bgd_sf %>% st_drop_geometry() %>%
  dplyr:: select(division, attack_rate, local_I, local_p, cluster)



write.csv(Local_Mor_I, "Local_Mor_I.csv", row.names = FALSE)





library(readxl); library(tidyverse); library(tseries)
library(forecast); library(patchwork)



dwr$date <- as.Date(dwr$date)

dwr <- dwr %>% filter(division != "Total") %>%
  mutate(division = droplevels(as.factor(division)))

divisions <- levels(dwr$division)




# ACF plots for all 8 divisions 

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1))
for (div in divisions) {
  x <- dwr %>%
    filter(division == div) %>%
    arrange(date) %>%
    pull(suspected_new_cases_24h) %>%
    na.omit()
  
  acf(
    x,
    main = paste(div),   
    col = "steelblue",
    lwd = 1.5,
    lag.max = 15
  )
}



dev.copy(
  png,
  "plot1a_ACF_all_divisions.png",
  width = 1400, height = 1100, res = 130
)
dev.off()






# PACF plots


par(mfrow = c(3, 3), mar = c(3, 3, 3, 1))
for (div in divisions) {
  x <- dwr %>%
    filter(division == div) %>%
    arrange(date) %>%
    pull(suspected_new_cases_24h) %>%
    na.omit()
  
  pacf(
    x,
    main = paste(div),   # show division name
    col = "firebrick",
    lwd = 1.5,
    lag.max = 15
  )
  
  
}


dev.copy(
  png,
  "plot1b_PACF_all_divisions.png",
  width = 1400, height = 1100, res = 130
)
dev.off()



# ADF Unit Root Test — Table 1 equivalent



adf_results <- map_dfr(divisions, function(div) {
  x <- dwr %>%
    filter(division == div) %>%
    arrange(date) %>%
    pull(suspected_new_cases_24h) %>%
    na.omit()
  
  # Level
  adf_lev <- adf.test(x, alternative = "stationary")
  # First difference
  adf_dif <- adf.test(diff(x), alternative = "stationary")
  
  tibble(
    Division       = div,
    ADF_level      = round(adf_lev$statistic, 4),
    p_level        = round(adf_lev$p.value, 4),
    Decision_level = ifelse(adf_lev$p.value < 0.05,
                            "Stationary","Unit root"),
    ADF_diff       = round(adf_dif$statistic, 4),
    p_diff         = round(adf_dif$p.value, 4),
    Decision_diff  = ifelse(adf_dif$p.value < 0.05,
                            "Stationary","Unit root")
  )
})


print(adf_results)

write.csv(adf_results, "adf_unit_root_results.csv", row.names = FALSE)








# BINOMIAL GLM

library(MASS); library(broom); library(scales)

# Population lookup

pop_lookup <- c(Barishal=8325666, Chattogram=28423019,
                Dhaka=36054418,   Khulna=15563000,
                Mymensingh=11370000, Rajshahi=18484858,
                Rangpur=15787000, Sylhet=9910219)

# Build binomial dataset: one row per division × date
# Response = (daily new cases, population_at_risk - daily new cases)
# Note: population_at_risk = total division population
# (Using cumulative cases as "successes" and pop-cumul as "failures"
# is the correct binomial framing matching the paper)



binom_data <- dwr %>%
  filter(division != "Total") %>%
  left_join(
    tibble(division = names(pop_lookup), pop_value = pop_lookup),
    by = "division"
  ) %>%
  mutate(
    pop_value = as.numeric(pop_value),
    suspected_cases_cumulative = as.numeric(suspected_cases_cumulative),
    suspected_new_cases_24h = as.numeric(suspected_new_cases_24h),
    cases = suspected_cases_cumulative,
    non_cases = pop_value - suspected_cases_cumulative,
    d_cases = suspected_new_cases_24h,
    d_noncases = pop_value - suspected_new_cases_24h
  ) %>%
  filter(!is.na(non_cases), non_cases > 0, cases >= 0)





# ── Binomial GLM (cumulative proportion) ─────────────────────

binom_model_cumul <- glm(
  cbind(cases, non_cases) ~ division,
  family  = binomial(link = "logit"),
  data    = binom_data
)

cat("=== BINOMIAL MODEL (cumulative) ===\n")
print(summary(binom_model_cumul))

# Goodness of fit

binom_disp <- deviance(binom_model_cumul) /
  df.residual(binom_model_cumul)
cat("Residual deviance/df:", round(binom_disp, 3), "\n")
cat("AIC:", round(AIC(binom_model_cumul), 2), "\n")

# ── Binomial GLM (daily proportion) ──────────────────────────

binom_model_daily <- glm(
  cbind(d_cases, d_noncases) ~ division,
  family  = binomial(link = "logit"),
  data    = binom_data
)

cat("\n=== BINOMIAL MODEL (daily) ===\n")
print(summary(binom_model_daily))





# FITTED vs OBSERVED per division — all 3 models


obs_fit <- binom_data %>%
  mutate(
    obs_rate  = suspected_new_cases_24h / 1000,
    fit_pois  = fitted(poisson_model) / 1000,
    fit_nb    = fitted(neg_binom_model) / 1000
  )

p_fit <- obs_fit %>%
  pivot_longer(cols = c(obs_rate, fit_pois, fit_nb),
               names_to = "type", values_to = "value") %>%
  mutate(type = recode(type,
                       obs_rate = "Observed",
                       fit_pois = "Poisson fitted",
                       fit_nb   = "NB fitted")) %>%
  ggplot(aes(date, value, colour = type, linetype = type)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  facet_wrap(~division, ncol = 4, scales = "free_y") +
  scale_colour_manual(
    values = c(Observed       = "black",
               "Poisson fitted" = "#C0392B",
               "NB fitted"    = "#1A5276"),
    name = "") +
  scale_linetype_manual(
    values = c(Observed = "solid",
               "Poisson fitted" = "dashed",
               "NB fitted" = "dotted"),
    name = "") +
  scale_x_date(date_labels = "%d %b",
               date_breaks = "10 days") +
  labs(title    = "Observed vs Fitted Values by Division",
       subtitle = "Poisson (dashed) and NB (dotted) fitted over observed (solid)",
       x = "Date", y = "Daily cases (thousands)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        strip.text    = element_text(face = "bold"),
        plot.title    = element_text(face = "bold", size = 12),
        legend.position = "top")

ggsave("plot6_fitted_vs_observed.png", p_fit,
       width = 14, height = 8, dpi = 180, bg = "white")














# Final models comparison 






library(readxl); library(MASS); library(tidyverse); library(broom)

dwr <- read_excel("Division_wise_data.xlsx", sheet="Long_Format")
dwr$date <- as.Date(dwr$date)
dwr <- dwr %>% filter(division != "Total") %>%
  mutate(division = droplevels(as.factor(division)))

pop <- c(Barishal=8325666,Chattogram=28423019,Dhaka=36054418,
         Khulna=15563000,Mymensingh=11370000,Rajshahi=18484858,
         Rangpur=15787000,Sylhet=9910219)

dwr <- dwr %>%
  left_join(tibble(division=names(pop), pop_count=pop), by="division") %>%
  mutate(cases     = as.integer(suspected_new_cases_24h),
         non_cases = as.integer(pop_count - suspected_new_cases_24h)) %>%
  filter(!is.na(cases), non_cases > 0, cases >= 0)

# Binomial
binom_model <- glm(cbind(cases, non_cases) ~ division,
                   family = binomial(link="logit"), data = dwr)

# Poisson (already done — re-run clean)
pois_model <- glm(cases ~ division,
                  family = poisson(link="log"), data = dwr)

# NB (already done)
nb_model <- glm.nb(cases ~ division, data = dwr)

# Comparison table
Comparison_table = tibble(
  Model          = c("Poisson","Neg. Binomial","Binomial"),
  AIC            = round(c(AIC(pois_model),AIC(nb_model),AIC(binom_model)),1),
  Resid_Dev_df   = round(c(deviance(pois_model)/df.residual(pois_model),
                           deviance(nb_model)/df.residual(nb_model),
                           deviance(binom_model)/df.residual(binom_model)),3),
  Delta_AIC      = round(c(AIC(pois_model),AIC(nb_model),AIC(binom_model)) -
                           min(c(AIC(pois_model),AIC(nb_model),AIC(binom_model))),1),
  Verdict        = c("Reject — overdispersed",
                     "Best count model",
                     "Best proportion model")
) %>% print()

write.csv(Comparison_table, "Comparison_table.csv", row.names = FALSE)

summary(binom_model)


sum_pois = summary(pois_model)

write.csv(sum_pois, "summmary_pois.csv", row.names = FALSE)




# IRR for all three models


make_irr_table <- function(model, model_name, exp_coef) {
  tidy(model, conf.int = TRUE,
       exponentiate = exp_coef) %>%
    mutate(
      model  = model_name,
      term   = gsub("division", "", term),
      IRR    = round(estimate, 4),
      CI_low = round(conf.low, 4),
      CI_high= round(conf.high, 4),
      p_val  = round(p.value, 4),
      sig    = case_when(p.value<0.001~"***",
                         p.value<0.01 ~"**",
                         p.value<0.05 ~"*",
                         TRUE~"ns")
    ) %>%
    rename(
      Model = model,
      Term = term,
      Estimate = IRR,
      CI_2.5 = CI_low,
      CI_97.5 = CI_high,
      p = p_val,
      Sig = sig
    )
  
  
}

full_table <- bind_rows(
  make_irr_table(pois_model,    "Poisson",           TRUE),
  make_irr_table(nb_model,  "Negative Binomial", TRUE),
  make_irr_table(binom_model_cumul,"Binomial",          TRUE)
)

cat("\n=== FULL IRR TABLE (all 3 models) ===\n")
print(full_table, n = 30)
write.csv(full_table, "irr_all_models.csv", row.names = FALSE)



















