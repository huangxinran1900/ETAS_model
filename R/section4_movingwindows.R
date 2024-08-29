#section4: moving time windows
#mu_p: time variation in mu
png(filename = "/Users/huangxinran/Desktop/mu_p_mu.png", width = 14, height = 6, units = "in", res = 300)

par(pin = c(9, 3))
wind_1 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu_p/output_file/mu_p_values_1day.csv")
wind_5 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu_p/output_file/mu_p_values_5days.csv")
wind_10 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu_p/output_file/mu_p_values_10days.csv")

plot(wind_1$Time..days., wind_1$Mu, type = "l", lty = 1, col = "black", lwd = 2,
     xlab = "Time", ylab = expression(mu(t)))
lines(wind_5$Time..days., wind_5$Mu, lty = 3, col = "red", lwd = 2)
lines(wind_10$Time..days., wind_10$Mu, lty = 2, col = "blue", lwd = 2)

# Horizontal reference line
abline(h = 1.139694, col = "black", lty = 2)

# Vertical reference line at the desired point
abline(v = wind_5$Time..days.[25], col = "gray", lty = 2)
arrows(x0 = wind_5$Time..days.[25] - 2, y0 = max(wind_1$Mu) * 0.8, 
       x1 = wind_5$Time..days.[25] - 0.5, y1 = max(wind_1$Mu) * 0.8, 
       col = "black", length = 0.1)

# Adjust the text position to align with the shorter arrow
text(x = wind_5$Time..days.[25] - 2.5, y = max(wind_1$Mu) * 0.8, 
     labels = expression("Aug 17"^th), adj = 1)
legend("topleft", legend = c("1 Day", "5 Days", "10 Days", "Stationary"),
       col = c("black", "red", "blue", "black"), lty = c(1, 3, 2, 2), lwd = 2)

# Finish the plot
dev.off()

#mu_p: time variation in p
png(filename = "/Users/huangxinran/Desktop/mu_p_p.png", width = 14, height = 6, units = "in", res = 300)

par(pin = c(9, 3))

plot(wind_1$Time..days., wind_1$P, type = "l", lty = 1, col = "black", lwd = 2,
     xlab = "Time", ylab = "p(t)")
lines(wind_5$Time..days., wind_5$P, lty = 3, col = "red", lwd = 2)
lines(wind_10$Time..days., wind_10$P, lty = 2, col = "blue", lwd = 2)


# Add a horizontal reference line
abline(h = 1.090183, col = "black", lty = 2)

legend("topleft", legend = c("1 Day", "5 Days", "10 Days", "stationary"),
       col = c("black", "red", "blue", "black"), lty = c(1, 3, 2, 2), lwd = 2)
dev.off()



#mu
png(filename = "/Users/huangxinran/Desktop/mu.png", width = 14, height = 6, units = "in", res = 300)

par(pin = c(9, 3))
wind_1 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu/output_file/mu_values_1day.csv")
wind_5 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu/output_file/mu_values_5days.csv")
wind_10 <- read.csv("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section4/vary_mu/output_file/mu_values_10days.csv")

plot(wind_1$Time..days.,wind_1$Mu, type = "l", lty = 1, col = "black", lwd = 2,
     xlab = "Time", ylab = expression(mu(t)))
lines(wind_5$Time..days., wind_5$Mu, lty = 3, col = "red", lwd = 2)
lines(wind_10$Time..days., wind_10$Mu, lty = 2, col = "blue", lwd = 2)

abline(h=1.139694, col="black", lty=2)
legend("topright", legend = c("1 Day", "5 Days","10 Days", "Stationary"),
       col = c("black", "red", "blue", "black"), lty = c(1, 3, 2, 2), lwd = 2)

dev.off()