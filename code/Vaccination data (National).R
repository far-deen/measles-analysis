

vd = read_excel("Vaccination data (National).xlsx")


names(vd)



# Age Distribution of Reported Cases 


age_data <- data.frame(
  AgeGroup = c("Under 2 years", "Under 5 years", "Under 9 months"),
  Percentage = c(66, 79, 33)
)

ggplot(age_data, aes(x = AgeGroup, y = Percentage, fill = AgeGroup)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.3, size = 5) +
  labs(
    title = "Age Distribution of Reported Cases",
    x = "Age Group",
    y = "Percentage of Reported Cases"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(10, 10, 30, 10),  
    axis.title.x = element_text(vjust = -1.5),
    plot.title    = element_text(hjust = 0.5)
  ) +
  ylim(0, 100)  # extends y-axis so the 79% label fits





































