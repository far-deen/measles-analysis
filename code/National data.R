
rdm= read.csv("National data.xlsx")

rdm

head(rdm)


pacman::p_load(
  rio,          # file import/export
  here,         # relative filepaths 
  lubridate,    # working with dates/epiweeks
  aweek,        # alternative package for working with dates/epiweeks
  incidence2,   # epicurves of linelist data
  i2extras,     # supplement to incidence2
  stringr,      # search and manipulate character strings
  forcats,      # working with factors
  RColorBrewer, # Color palettes from colorbrewer2.org
  tidyverse     # data management + ggplot2 graphics
) 


library(incidence2)



# Convert to class date

rdm = rdm %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

colnames(df)




# Suspected Cases - Daily

ggplot(rdm, aes(x = date, y = suspected_new_cases_24h)) +
  geom_col(fill = "steelblue") +
  scale_x_date(
    breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Suspected Cases - Daily",
    x = "Date",
    y = "Number of Suspected Cases (24h)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))





# Confirmed Cases - Daily

ggplot(rdm, aes(x = date, y = confirmed_new_cases_24h)) +
  geom_col(fill = "steelblue") +
  scale_x_date(
    breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Confirmed Cases - Daily",
    x = "Date",
    y = "Number of Confirmed Cases (24h)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Confirmed Deaths - Daily

ggplot(rdm, aes(x = date, y = confirmed_new_deaths_24h)) +
  geom_col(fill = "steelblue") +
  scale_x_date(
    breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Confirmed Deaths - Daily",
    x = "Date",
    y = "Number of Deaths (24h)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# Suspected Cases - Cumulative

ggplot(rdm, aes(x = date, y = suspected_cases_cumulative)) +
  geom_col(fill = "steelblue") +
  scale_x_date(
    breaks = seq(min(rdm$date) , max(rdm$date), by = "1 day"),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Suspected Cases - Cumulative",
    x = "Date",
    y = "Number of Suspected Cases (Cumulative)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))











# Average Daily Reported Measles Cases by Week


weekly = read_excel("Extracted_Data.xlsx", sheet = "Weekly_Trend_Data")



ggplot(weekly, aes(reorder(period_label, week_number), avg_daily_cases)) +
  geom_col(fill="#2E75B6", width=0.65) +
  geom_text(aes(label=comma(avg_daily_cases)), vjust=-0.4, size=4.5, fontface="bold") +
  labs(title="Average Daily Reported Measles Cases by Week",
       x="", y="Average cases per day") +
  theme_minimal(base_size=12)+
  theme(plot.title    = element_text(hjust = 0.5))




























































































