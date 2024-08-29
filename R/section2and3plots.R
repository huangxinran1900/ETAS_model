# Necessary R packages
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(ggspatial)
library(patchwork)
library(grid)
library(gtable)

# Preprocess the raw dataset to run in the ZMAP software to analyze the completeness of magnitude
preprocess_for_zmap <- function(input_csv, output_csv) {
  # Read the input CSV file
  Taiwan <- read.csv(input_csv)
  
  # Remove unnecessary columns and rename the required ones
  Taiwan <- subset(Taiwan, select = -c(No., Place_2, Location))
  Taiwan <- Taiwan %>%
    rename(time=Orgin.date, longitude=Longitude.E., latitude=Latitude.N., mag=Magnitude, depth=Depth, place=Place)
  
  # Drop the 'place' column
  taiwan_data <- subset(Taiwan, select = -c(place))
  
  # Convert the 'time' column to POSIXct format
  taiwan_data$time <- as.POSIXct(taiwan_data$time, format = "%Y-%m-%d %H:%M:%S")
  
  # Extract year, month, day, hour, minute, and second
  taiwan_data$year <- as.numeric(format(taiwan_data$time, "%Y"))
  taiwan_data$month <- as.numeric(format(taiwan_data$time, "%m"))
  taiwan_data$day <- as.numeric(format(taiwan_data$time, "%d"))
  taiwan_data$hour <- as.numeric(format(taiwan_data$time, "%H"))
  taiwan_data$minute <- as.numeric(format(taiwan_data$time, "%M"))
  taiwan_data$second <- as.numeric(format(taiwan_data$time, "%S"))
  
  # Remove the 'time' column
  taiwan_data <- subset(taiwan_data, select = -c(time))
  
  # Reorder the columns
  taiwan_data <- taiwan_data[, c("longitude", "latitude", "year", "month", "day", "mag", "depth", "hour", "minute", "second")]
  
  # Write the table without column names
  write.table(taiwan_data, file = output_csv, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Full raw dataset to detect the completeness of magnitude in ZMAP
preprocess_for_zmap("/Users/huangxinran/Desktop/Code/data/complete_dataset/Taiwan_updated.csv", "/Users/huangxinran/Desktop/taiwan_Mc.csv")



# Process the dataset above Mc into work.etas file to run in the Fortran program
preprocess_for_etas <- function(input_csv) {
  # Read the input CSV file
  Taiwan <- read.csv(input_csv)
  
  # Remove unnecessary columns and rename the required ones
  Taiwan <- subset(Taiwan, select = -c(No., Place_2, Location))
  Taiwan <- Taiwan %>%
    rename(time = Orgin.date, longitude = Longitude.E., latitude = Latitude.N., mag = Magnitude, depth = Depth, place = Place)
  
  # Drop the 'place' column
  taiwan_data <- subset(Taiwan, select = -c(place))
  
  # Convert the 'time' column to POSIXct format
  taiwan_data$time <- as.POSIXct(taiwan_data$time, format = "%Y-%m-%d %H:%M:%S")
  
  # Extract year, month, and day
  taiwan_data$year <- as.numeric(format(taiwan_data$time, "%Y"))
  taiwan_data$month <- as.numeric(format(taiwan_data$time, "%m"))
  taiwan_data$day <- as.numeric(format(taiwan_data$time, "%d"))
  
  # Order the data by time and add a sequential number
  taiwan_data <- taiwan_data[order(taiwan_data$time), ]
  taiwan_data$seq_num <- 1:nrow(taiwan_data)
  
  # Calculate days from t0
  t0 <- min(as.POSIXct(taiwan_data$time))
  taiwan_data$days_from_t0 <- as.numeric(difftime(taiwan_data$time, t0, units = "days"))
  
  # Remove the time column
  taiwan_data_0 <- subset(taiwan_data, select = -c(time))
  
  # Create the final dataset in the required format
  taiwan_data <- taiwan_data_0[, c("seq_num", "longitude", "latitude", "mag", "days_from_t0", "depth", "year", "month", "day")]
  
  # Rename columns to ensure they are correctly named
  colnames(taiwan_data) <- c("seq_num", "longitude", "latitude", "mag", "time", "depth", "year", "month", "day")
  
  # Convert the depth column to negative values
  taiwan_data[, 6] <- -abs(taiwan_data[, 6])
  return(taiwan_data)
}

# Preprocess the dataset beyond Mc to run into the Fortran program
taiwan_data<-preprocess_for_etas("/Users/huangxinran/Desktop/Code/data/complete_dataset/taiwan_complete.csv")

# Write the table without column names
write.table(taiwan_data, file =  "/Users/huangxinran/Desktop/work.etas", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# Create initial input and setting necessary to run the Fortran program
# Define a function to create the etas.open file
create_etas_open <- function(output_file_path, 
                             nfunct, 
                             iappr, 
                             zts, 
                             zte, 
                             tstart, 
                             amx1, 
                             xmag0, 
                             xini1, 
                             neg_log_likelihood, 
                             xini2) {
  
  # Create the etas.open content
  etas_open_content <- c(
    paste(nfunct, iappr),
    paste(zts, zte, tstart),
    paste(amx1, xmag0),
    paste(xini1, collapse = " "),
    paste(zts, zte, tstart),
    paste(amx1, xmag0),
    neg_log_likelihood,
    paste(xini2, collapse = " "),
    paste(nfunct, iappr),
    paste(zts, zte, tstart),
    paste(amx1, xmag0),
    paste(xini2, collapse = " "),
    paste(zts, zte, tstart),
    paste(amx1, xmag0),
    neg_log_likelihood
  )
  
  # Write to file
  writeLines(etas_open_content, output_file_path)
}

create_etas_open(
  output_file_path = "/Users/huangxinran/Desktop/etas.open",
  nfunct = 9,   ##subrountine func9 (fast log likelihood)
  iappr = 2,    ## level 2 (lower level has lower computation efficiency but higher accuracy)
  zts = 0.0,    ## start time
  zte = 3653.85103009259,   ## end time
  tstart = 5.01,   ## start of study time
  amx1 = 3.6,     ## Mc (completeness of mag)
  xmag0 = 7.2,    ## Mz (largest magnitude among dataset)
  xini1 = c(0.34378E+00, 0.24709E+00, 0.49796E-02, 0.63649E+00, 0.11189E+01),   ## initial value of mu, K0, c, alpha, p
  neg_log_likelihood = -0.4753056E+04, 
  xini2 = c(0.34378E+00, 0.24709E+00, 0.49796E-02, 0.63649E+00, 0.11189E+01)
)


# Process the dataset of magnitude
# complete dataset of whole period (Mag>1.5)
Taiwan <- read.csv("/Users/huangxinran/Desktop/Code/data/complete_dataset/Taiwan_updated.csv")
Taiwan <- subset(Taiwan, select = -c(No., Place, Place_2, Location))
Taiwan <- Taiwan %>%
  rename(time=Orgin.date, longitude=Longitude.E.,latitude=Latitude.N., mag=Magnitude, depth=Depth)
taiwan_data_plot <- Taiwan %>%
  mutate(time=ymd_hms(time)) %>%
  dplyr::select(time, latitude, longitude, mag, depth)
taiwan_data_plot$time <- as.POSIXct(taiwan_data_plot$time, format="%Y-%m-%d %H:%M:%S")
taiwan_data_plot <- taiwan_data_plot[order(taiwan_data_plot$time), ]
taiwan_data_plot <- subset(taiwan_data_plot, latitude >= 23 & latitude <= 25)
taiwan_data_plot <- subset(taiwan_data_plot, longitude >= 120.5 & latitude <= 123.1)
high_magnitude_data <- subset(taiwan_data_plot, mag > 6)
mag_1 <- high_magnitude_data[5,]
mag_2 <- high_magnitude_data[22,]
high_magnitude_data <- rbind(mag_1, mag_2)

df_filtered <- taiwan_data_plot[taiwan_data_plot$time != high_magnitude_data$time, ]
df_sf_filtered <- st_as_sf(df_filtered, coords = c("longitude", "latitude"), crs = 4326)
high_mag_coords <- high_magnitude_data[, c("longitude", "latitude")]

# Convert the extracted data to a Simple Features object
high_mag_sf <- st_as_sf(high_mag_coords, coords = c("longitude", "latitude"), crs = 4326)
world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
# Region_analysis plot
main_plot <- ggplot() +
  geom_sf(data = world, fill = "gray80", color = "white") +
  geom_sf(data = df_sf_filtered, aes(color = mag)) +
  geom_point(data = mag_1, aes(x = longitude, y = latitude), 
             shape = 17, size = 3, color = "red", stroke = 2) +
  geom_point(data = mag_2, aes(x = longitude, y = latitude), 
             shape = 24, size = 3, color = "red", stroke = 2) +
  scale_color_gradient(low = "yellow", high = "red", name = "Magnitude") +
  coord_sf(xlim = c(min(taiwan_data_plot$longitude), max(taiwan_data_plot$longitude)), ylim = c(min(taiwan_data_plot$latitude), max(taiwan_data_plot$latitude))) +
  labs(title = "Earthquake Magnitudes",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()+
  theme_void()+
  theme(
    panel.grid = element_blank())+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    aspect.ratio = 1.1/1 
  )
plot(main_plot)

taiwan_map <- st_as_sf(maps::map("world", regions = "Taiwan", plot = FALSE, fill = TRUE))

inset_plot <- ggplot() +
  geom_sf(data = taiwan_map, fill = "gray80", color = "white") +
  geom_rect(data = taiwan_data_plot, aes(xmin = min(taiwan_data_plot$longitude), xmax = max(taiwan_data_plot$longitude), ymin = min(taiwan_data$latitude), ymax = max(taiwan_data$latitude)),
            fill = NA, color = "black", size = 1) +
  coord_sf(xlim = c(119.9, 122.05), ylim = c(21.93, 25.3)) +  
  theme_void()  +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
print(inset_plot)
vp <- viewport(width = 0.2, height =0.29, x = 0.1589,y = 0.835,just=c("left","top"))
print(main_plot)
print(inset_plot,vp=vp)
save_combined_plot <- function(filename, main_plot, inset_plot, vp, width, height, dpi) {
  png(filename, width = width, height = height, units = "in", res = dpi)
  grid.newpage()
  print(main_plot)
  # Draw the inset plot in the specified viewport
  print(inset_plot, vp = vp)
  dev.off()
}
save_combined_plot("/Users/huangxinran/Desktop/region_analysis.png", main_plot, inset_plot, vp, width = 6.5, height = 7.9, dpi = 300)

# Produce plots of Fortran output (Ogata: SASeis2006) of stationary ETAS model
# Plot1: Cumulative number of events versus transformed time
png(filename = "/Users/huangxinran/Desktop/stationary_residual_plot.png", width = 14, height = 6, units = "in", res = 300)
para <- scan("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/stationary_ETAS/output_file/etas.open")
mag1 <- para[6]
fmt <- scan("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/stationary_ETAS/output_file/work.res")
hypo <- matrix(fmt, ncol = 7, byrow = TRUE)
hypo <- hypo[hypo[,4]>0,]
level <- min(0, fmt[1])
cn <- 1:length(hypo[, 1])+hypo[1,1]
ti <- hypo[, 7]
xrange <- c(min(ti), max(ti))
mgrange <- max(cn)/4
bottom <- min(cn) - mgrange
yrange <- c(bottom, max(max(cn),max(ti)))
par(pty="s",xaxs='r')
plot(xrange, yrange, type = "n", main = "ETAS Residual", xlab
     = "TRANSFORMED TIME", ylab = "CUMULATIVE NUMBER OF EVENTS",
     lwd = 1, xaxt = "n", yaxt = "n")
axis(1, at = seq(from = 0, to = 6000, length.out = 8), labels = seq(from = 0, to = 3500, length.out = 8))
axis(2, at = seq(from = 0, to = 6000, length.out = 7), labels = seq(from = 0, to = 3000, length.out = 7))

points(ti, 1:length(ti)+hypo[1,1]+1, type='s')
mag <- hypo[, 4]
mgmax <- max(mag - mag1 + 1)
mag <- mag - mag1 + 0.5
segments(ti, bottom, ti, mag/mgmax * mgrange + bottom)	
abline(h=bottom)
abline(h=0)
abline(v=0,lty=2)
abline(0,1,lty=1,col='red')
s <- para[5]
t <- para[14]
timax <- max(hypo[hypo[,5]<=t, 7])
abline(v=timax,lty=2)
te <- para[4]
text(max(ti)*0.4,max(cn)*0.9,paste('S=',s,' T=',t))
text(max(ti)*0.4,max(cn)*.83,paste('M>=',mag1,' Tend=',te))
mu <- para[8]
k0 <- para[9]
c <- para[10]
a <- para[11]
p <- para[12]

conf_level <- 0.99  # 99% confidence level

n <- length(ti)

# Confidence bands for k < n (using beta distribution)
conf_upper <- sapply(1:n, function(k) qbeta(1 - (1 - conf_level) / 2, k+1, n - k + 1))
conf_lower <- sapply(1:n, function(k) qbeta((1 - conf_level) / 2, k+1, n - k + 1))
conf_upper <- conf_upper * n
conf_lower <- conf_lower * n


# Add the confidence bands
lines(1:n, conf_upper + hypo[1,1] + 1, lty = 2, col = "blue")
lines(1:n, conf_lower + hypo[1,1] + 1, lty = 2, col = "blue")

# Add boxes and labels
anomalies <- list(
  list(start=1610, end=2350, bottom= bottom+2770, top = max(cn)-3380, label="2018\nM6.2", line_pos = 1611.286),
  list(start=3130, end=4200, bottom= bottom+4250, top = max(cn)-1900, label="2021\nApr", line_pos = 3120),
  list(start=4200, end=5100, bottom= bottom+5500, top = max(cn)-840, label="2024\nM7.2", line_pos = 4819)
)
max_y_for_labels <- max(cn) + 500
for (anomaly in anomalies) {
  rect(anomaly$start, anomaly$bottom, anomaly$end, anomaly$top, border="grey", lwd=2, lty=2)
  mtext(anomaly$label, side = 3, line = 0, at = anomaly$line_pos, cex = 0.8)
}

# 2018-02-06: 1327.9612
hypo[hypo[,5]==1327.9612,]
abline(v = 1611.286, col = "gray", lwd = 2, lty = 2)
# 2021: swarms
abline(v = 3120, col = "gray", lwd = 2, lty = 2)
# 2024: Hualien
abline(v = 4819.6182, col = "grey", lwd = 2, lty = 2)
dev.off()

##########
# Plot2: Cumulative number of events versus ordinary time
png(filename = "/Users/huangxinran/Desktop/combined_plots.png", width = 14, height = 6, units = "in", res = 300)

para <- scan("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/stationary_ETAS/output_file/etas.open")
mag0 <- para[6]; ts <- para[3]; te <- para[4]
data <- scan("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/stationary_ETAS/work.etas",skip=1)

hypo <- matrix(data, ncol = 9, byrow = T)
hypo <- hypo[hypo[, 4] >= mag0,  ]
mag1 <- min(hypo[, 4])
mag <- hypo[, 4]
ti <- hypo[, 5]
zero <- hypo[, 1] * 0
cn <- 1:length(hypo[, 1])
cna <- append(cn, 0, after = 0)
cn1 <- cna[1:length(cn)]
tia <- append(ti, 0, after = 0)
ti1 <- tia[1:length(ti)]
#       par(xaxs='i',yaxs='i')
par(pty="s",xaxs='r',yaxs='r')
xrange <- c(min(ti), max(ti))
mgrange <- max(cn)/4
bottom <- min(cn) - mgrange
inc <- matrix(scan("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/stationary_ETAS/output_file/work.res"),ncol=7,byrow=T)
yrange <- c(bottom, max( max(cn),max(inc[,7]-inc[1,1]) ))
plot(xrange,yrange, type = "n", 
     main = "ETAS Fit and Prediction",
     xlab = "ORDINARY TIME (DAYS)", ylab = "CUMULATIVE NUMBER OF EVENTS", lwd = 1, xlim=c(xrange))
#ETAS increasing function:
lines(inc[,5],inc[,7]-inc[1,1],type='l',col='red')
segments(ti1, cn1, ti, cn1)
segments(ti, cn1, ti, cn)
mgmax <- max(mag - mag1 + 1)
mag <- mag - mag1 + 0.5
segments(ti, bottom, ti, mag/mgmax * mgrange + bottom)	
#	text(ti, mag/mgmax * mgrange + bottom, "-")
abline(h = 0)
abline(h = bottom)
mark0 <- para[ 5]; abline(v = mark0,lty=2)
mark1 <- para[18]; abline(v = mark1,lty=2)
mark2 <- para[ 4]; abline(v = mark2,lty=2)
s <- para[5]
t <- para[14]
t0 <- para[5]
abline(v=t0,lty=2)
abline(v=t,lty=2)
te <- para[4]
text(max(ti)*0.4,max(cn)*0.9,paste('S=',s,' T=',t))
text(max(ti)*0.4,max(cn)*.83,paste('M>=',mag1,' Tend=',te))
mu <- para[8]
k0 <- para[9]
c  <- para[10]
a  <- para[11]
p  <- para[12]
dev.off()


############
# Single change point analysis
## subperiod 1 transformed time plot [5.01, 3570.901]
# Convert to cumulative counts
png(filename = "/Users/huangxinran/Desktop/sub_1.png", width = 7, height = 6, units = "in", res = 300)

generate_etas_plot <- function(etas_open_path, work_res_path, confidence_level = 0.99) {
  # Read parameters and data
  para <- scan(etas_open_path)
  mag1 <- para[6]
  fmt <- scan(work_res_path)
  
  # Process the data
  hypo <- matrix(fmt, ncol = 7, byrow = TRUE)
  hypo <- hypo[hypo[,4] > 0, ]
  cn <- 1:length(hypo[, 1]) + hypo[1,1]
  ti <- hypo[, 7]
  
  # Set up plot ranges
  xrange <- c(min(ti), max(ti))
  mgrange <- max(cn) / 4
  bottom <- min(cn) - mgrange
  yrange <- c(bottom, max(max(cn), max(ti)))
  
  # Generate the plot
  par(pty = "s", xaxs = 'r')
  plot(xrange, yrange, type = "n", main = "ETAS Residual", 
       xlab = "TRANSFORMED TIME", ylab = "CUMULATIVE NUMBER OF EVENTS", lwd = 1)
  points(ti, 1:length(ti) + hypo[1,1] + 1, type = 's')
  
  mag <- hypo[, 4]
  mgmax <- max(mag - mag1 + 1)
  mag <- mag - mag1 + 0.5
  segments(ti, bottom, ti, mag / mgmax * mgrange + bottom)
  
  abline(h = bottom)
  abline(h = 0)
  abline(v = 0, lty = 2)
  abline(0, 1, lty = 1, col = 'red')
  
  # Annotate the plot
  s <- para[5]
  t <- para[14]
  timax <- max(hypo[hypo[,5] <= t, 7])
  abline(v = timax, lty = 2)
  te <- para[4]
  text(max(ti) * 0.4, max(cn) * 0.9, paste('S=', s, ' T=', t))
  text(max(ti) * 0.4, max(cn) * 0.81, paste('M>=', mag1, ' Tend=', te))

}
generate_etas_plot("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/single_changepoint/subperiod_1/output_file/etas.open", "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/single_changepoint/subperiod_1/output_file/work.res")
conf_level <- 0.99  # 99% confidence level
n <- length(ti[ti <= timax])

# Confidence bands for k < n (using beta distribution)
conf_upper <- sapply(1:n, function(k) qbeta(1 - (1 - conf_level) / 2, k+1, n - k + 1))
conf_lower <- sapply(1:n, function(k) qbeta((1 - conf_level) / 2, k+1, n - k + 1))
conf_upper <- conf_upper * n
conf_lower <- conf_lower * n

# Add the confidence bands
lines(1:n, conf_upper + hypo[1,1] + 1, lty = 2, col = "blue")
lines(1:n, conf_lower + hypo[1,1] + 1, lty = 2, col = "blue")
dev.off()

# subperiod 2 transformed time plot [3570.901, 3653.85103009259]
png(filename = "/Users/huangxinran/Desktop/sub_2.png", width = 7, height = 6, units = "in", res = 300)

generate_etas_plot("/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/single_changepoint/subperiod_2/output_file/etas.open", "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/single_changepoint/subperiod_2/output_file/work.res")
# Confidence bands for k < n (using beta distribution)
conf_upper <- sapply(1:n, function(k) qbeta(1 - (1 - conf_level) / 2, k+1, n - k + 1))
conf_lower <- sapply(1:n, function(k) qbeta((1 - conf_level) / 2, k+1, n - k + 1))
conf_upper <- conf_upper * n
conf_lower <- conf_lower * n

# Add the confidence bands
lines(1:n, conf_upper + hypo[1,1] + 1, lty = 2, col = "blue")
lines(1:n, conf_lower + hypo[1,1] + 1, lty = 2, col = "blue")
dev.off()


###############
# multistage ETAS model plots
### subperiod_11: [5.01,1353.034]
png(filename = "/Users/huangxinran/Desktop/sub_11_section.png", width = 5, height = 6, units = "in", res = 300)
generate_etas_section_plot <- function(etas_open_path, work_res_path) {
  para <- scan(etas_open_path)
  mag1 <- para[6]
  
  # Read data from the work.res file
  fmt <- scan(work_res_path)
  hypo <- matrix(fmt, ncol = 7, byrow = TRUE)
  hypo <- hypo[hypo[,4] > 0, ]
  
  # Calculate values for plotting
  cn <- 1:length(hypo[, 1]) + hypo[1, 1]
  ti <- hypo[, 7]
  xrange <- c(min(ti), max(ti))
  mgrange <- max(cn) / 4
  bottom <- min(cn) - mgrange
  yrange <- c(bottom, max(max(cn), max(ti)))
  
  # Determine timax for the plotting range
  timax <- max(hypo[hypo[, 5] <= para[14], 7])
  
  # Define a layout matrix for larger upper plot
  layout_matrix <- matrix(c(1, 1, 2), nrow = 3, byrow = TRUE) # 2 rows for the first plot, 1 row for the second plot
  layout(layout_matrix)
  
  # Plot the main ETAS residuals plot (larger plot)
  par(mar = c(4, 3, 2, 1) + 0.1)
  plot(xrange, yrange, type = "n", main = "ETAS Residual", 
       xlab = "TRANSFORMED TIME", ylab = "CUMULATIVE NUMBER OF EVENTS", lwd = 1)
  points(ti, 1:length(ti) + hypo[1, 1] + 1, type = 's')
  abline(h = bottom)
  abline(h = 0)
  abline(v = 0, lty = 2)
  abline(0, 1, lty = 1, col = 'red')
  abline(v = timax, lty = 2)
  
  # Add text annotations for model parameters
  s <- para[5]
  t <- para[14]
  te <- para[4]
  text(max(ti) * 0.4, max(cn) * 0.9, paste('S =', s, ' T =', t))
  text(max(ti) * 0.4, max(cn) * 0.81, paste('M >=', mag1, ' Tend =', te))
  
  # Plot the magnitude segments (smaller plot)
  # Filter to only include data up to timax
  filtered_indices <- ti <= timax & ti >= 0
  ti_filtered <- ti[filtered_indices]
  mag_filtered <- hypo[filtered_indices, 4]
  
  # Plot for magnitudes
  par(mar = c(4, 4, 1, 1) + 0.1) # Adjust margins for better appearance
  plot(ti_filtered, mag_filtered, type = "n", main = "Magnitude Plot", 
       xlab = "TRANSFORMED TIME", ylab = "MAGNITUDE", ylim = c(min(mag_filtered), max(mag_filtered) + 1))
  segments(ti_filtered, 0, ti_filtered, mag_filtered, col = "black")
}

# Example usage:
generate_etas_section_plot(
  etas_open_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_1/output_file/etas.open",
  work_res_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_1/output_file/work.res")
dev.off()

### subperiod_21: [1353.034,2537.57]
png(filename = "/Users/huangxinran/Desktop/sub_21_section.png", width = 5, height = 6, units = "in", res = 300)

generate_etas_section_plot(
  etas_open_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_2/output_file/etas.open",
  work_res_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_2/output_file/work.res")
dev.off()

### subperiod_22: [2537.57,3570.901]
png(filename = "/Users/huangxinran/Desktop/sub_22_section.png", width = 5, height = 6, units = "in", res = 300)

generate_etas_section_plot(
  etas_open_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_3/output_file/etas.open",
  work_res_path = "/Users/huangxinran/Desktop/Code/Fortran_prgram/Section3/multistage_changepoints/subperiod_3/output_file/work.res")
dev.off()


# Plot only the desired section for [2537.57,3570.901] and the confidence band
# Filter indices based on the desired time range
filtered_indices <- ti >= 0 & ti <= timax

# Filter the time and cumulative number vectors
ti_filtered <- ti[filtered_indices]
cn_filtered <- cn[filtered_indices]
png(filename = "/Users/huangxinran/Desktop/sub_12_section.png", width = 7, height = 6, units = "in", res = 300)
par(pty = "s", xaxs = 'r')
# Update xrange and yrange for the filtered data
xrange <- range(ti_filtered)  # Range of the filtered time values
yrange <- range(cn_filtered)  # Range of the filtered cumulative number of events

plot(xrange, yrange, type = "n", 
     xlab = "TRANSFORMED TIME", ylab = "CUMULATIVE NUMBER OF EVENTS", lwd = 1)

# Plot the filtered cumulative number of events
points(ti_filtered, cn_filtered, type = 's', col = 'black')
abline(0,1,lty=1,col='red')

# Plot magnitude bars for the filtered data
mag_filtered <- hypo[filtered_indices, 4]
mgmax <- max(mag_filtered - mag1 + 1)
mag_scaled <- mag_filtered - mag1 + 0.5

# Plot reference lines for clarity
abline(h = 0)
abline(v = 0, lty = 2)
abline(v = timax, lty = 2)

# Add text annotations
s <- para[5]
t <- para[14]
te <- para[4]

conf_level <- 0.99  # 99% confidence level

n <- sum(filtered_indices)
# Confidence bands for k < n (using beta distribution)
conf_upper <- sapply(1:n, function(k) qbeta(1 - (1 - conf_level) / 2, k+1, n - k + 1))
conf_lower <- sapply(1:n, function(k) qbeta((1 - conf_level) / 2, k+1, n - k + 1))
conf_upper <- conf_upper * n
conf_lower <- conf_lower * n
# Add the confidence bands
lines(1:n, conf_upper, lty = 2, col = "blue")
lines(1:n, conf_lower, lty = 2, col = "blue")
dev.off()

