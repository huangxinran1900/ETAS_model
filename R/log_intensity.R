### log-intensity plot of 2018 earthquake
# parameters from the stationary ETAS model

Taiwan <- read.csv("/Users/huangxinran/Desktop/Code/data/complete_dataset/taiwan_complete.csv")

Taiwan <- subset(Taiwan, select = -c(No., Place_2, Location))
Taiwan <- Taiwan %>%
  rename(time=Orgin.date, longitude=Longitude.E.,latitude=Latitude.N., mag=Magnitude, depth=Depth, place=Place)
taiwan_data <- subset(Taiwan, select = -c(place))
taiwan_data$time <- as.POSIXct(taiwan_data$time, format="%Y-%m-%d %H:%M:%S")
# Assuming taiwan_data has been loaded with the necessary columns and cleaned

# Extract year, month, and day from the timestamp
taiwan_data$year <- as.numeric(format(taiwan_data$time, "%Y"))
taiwan_data$month <- as.numeric(format(taiwan_data$time, "%m"))
taiwan_data$day <- as.numeric(format(taiwan_data$time, "%d"))

# Order data by time
taiwan_data <- taiwan_data[order(taiwan_data$time), ]
taiwan_data$seq_num <- 1:nrow(taiwan_data)

# Calculate days from the first event (t0)
t0 <- min(as.POSIXct(taiwan_data$time))
taiwan_data$days_from_t0 <- as.numeric(difftime(taiwan_data$time, t0, units = "days"))

# Initialize the fixed parameters
mu <- 0.164
K0 <- 2.293
c <- 0.003
alpha <- 1.257
p <- 1.051

intensity <- numeric(nrow(taiwan_data))

# Calculate the conditional intensity function for each event
for (i in 1:nrow(taiwan_data)) {
  t <- taiwan_data$days_from_t0[i]
  
  if (i == 1) {
    # For the first event, only the background rate applies
    intensity[i] <- mu
  } else {
    # Calculate intensity based on all previous events (t_i < t)
    contributions <- sapply(1:(i-1), function(j) {
      if (t > taiwan_data$days_from_t0[j]) {
        K0 * exp(alpha * (taiwan_data$mag[j] - 3.6)) * (t - taiwan_data$days_from_t0[j] + c)^(-p)
      } else {
        0
      }
    })
    
    # Calculate total intensity as background rate plus contributions from previous events
    intensity[i] <- mu + sum(contributions)
  }
}

# Add the intensity values back to the data frame
taiwan_data$intensity <- log(intensity)

# Filter the data to include only the events between "2018-01-20" and "2018-02-20"
start_date <- as.Date("2018-01-20")
end_date <- as.Date("2018-03-04")
filtered_data <- taiwan_data[taiwan_data$time >= start_date & taiwan_data$time <= end_date, ]

line1_position <- 1325.11416666667
line2_position <- 1327.9612037037
png(filename = "/Users/huangxinran/Desktop/intensity.png", width = 4, height = 4, units = "in", res = 300)

ggplot(filtered_data, aes(x = days_from_t0, y = intensity)) +
  geom_line(color = "black") +
  geom_vline(xintercept = line1_position, linetype = "dashed", color = "red") +
  geom_vline(xintercept = line2_position, linetype = "dashed", color = "red") +
  geom_text(aes(x = line1_position, y = max(filtered_data$intensity), label = "June 4th"),
            vjust = -0.9, hjust = 1.1, angle = 90, color = "black") +
  geom_text(aes(x = line2_position, y = max(filtered_data$intensity), label = "June 6th"),
            vjust = -0.9, hjust = 1.1, angle = 90, color = "black") +
  labs(x = "Time (days from t0)",
       y = expression(log(lambda * "(" * t ~ "|" ~ H[t] * ")")))+
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", size = 0.4), # Add a black border around the plot
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_text(size = 8)  )

dev.off()