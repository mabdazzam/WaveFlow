# Load required libraries
# The `fftw` library is used for performing fast Fourier transform operations.
library(fftw)

# Section 1: Download and Load Flow Data

# Define the URL for downloading daily flow data from the USGS NWIS website.
url <- "https://nwis.waterdata.usgs.gov/nwis/dv?cb_00010=on&cb_00060=on&cb_00065=on&cb_00095=on&cb_00300=on&format=rdb&site_no=02334430&legacy=&referred_module=sw&period=&begin_date=1950-01-01&end_date=2023-12-31"

# Define file path to save the downloaded data.
file_path <- "usgs_02334430_daily_data.txt"

# Download the file and save it locally.
download.file(url, destfile = file_path, method = "auto")

# Read the tab-separated data, skip the first metadata row, and clean column names.
data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-1, ]
rownames(data) <- 1:nrow(data)  # Reset row names for clarity.

# Convert the date column to `Date` type and rename relevant columns.
data$date <- as.Date(data[, 3], format = "%Y-%m-%d")
colnames(data)[c(22)] <- "q"  # Rename column 22 to `q` (discharge/streamflow).
data$q <- as.numeric(data[, "q"])  # Ensure the `q` column is numeric.

# Extract streamflow and date as a time series object.
q <- data[, "q"]
q.timeseries <- data[, c("date", "q")]

# Section 2: Apply Fast Fourier Transform (FFT) and Analyze Seasonal Signals

# Remove any missing values from the streamflow data.
q <- na.omit(q)

# Perform FFT to analyze frequency components.
fft_result <- fft(q, inverse = FALSE)

# Extract the magnitude of the FFT result for plotting.
fft_magnitude <- abs(fft_result)

# Calculate the frequencies corresponding to the FFT result.
frequencies <- seq(0, length(q) - 1) / length(q)

# Plot the original streamflow data.
plot(q, type = "l", xlab = "Time", ylab = "Streamflow", main = "Original Streamflow Data")

# Plot the FFT magnitude to visualize dominant frequencies.
plot(frequencies, fft_magnitude, type = "l", xlab = "Frequency", ylab = "Magnitude", main = "FFT Magnitude Spectrum")

# Section 3: Reconstruct and Validate Original Time Series

# Extract magnitude and phase information for inverse FFT.
fft_magnitude <- Mod(fft_result)
fft_phase <- Arg(fft_result)

# Reconstruct the signal using inverse FFT.
reconstructed_signal <- fft(fft_result, inverse = TRUE)
reconstructed_signal <- Re(reconstructed_signal) / length(q)

# Plot original and reconstructed streamflow data for comparison.
plot(q, type = "l", col = "blue", xlab = "Time", ylab = "Streamflow", main = "Original vs. Reconstructed Streamflow")
lines(reconstructed_signal, col = "red")
legend("topright", legend = c("Original", "Reconstructed"), col = c("blue", "red"), lty = 1)

# Calculate and print correlation between the original and reconstructed data.
if (length(q) == length(reconstructed_signal)) {
  correlation_value <- cor(q, reconstructed_signal)
  print(paste("Correlation between original and reconstructed data:", correlation_value))
} else {
  print("The lengths of the original and reconstructed data do not match. Cannot calculate correlation.")
}

# Section 4: Denoising with Moving Average and FFT

# Define window size for smoothing (e.g., 15-day moving average).
window_size <- 15
filter_weights <- rep(1 / window_size, window_size)  # Uniform filter weights.

# Apply moving average filter for smoothing.
q_smoothed <- stats::filter(q, filter_weights, sides = 2)

# Remove NA values introduced by the filter.
q_smoothed <- na.omit(q_smoothed)

# Perform FFT on the smoothed data.
fft_result <- fft(q_smoothed, inverse = FALSE)
fft_magnitude <- abs(fft_result)  # Extract magnitude.

# Reconstruct the smoothed signal using inverse FFT.
reconstructed_signal <- fft(fft_result, inverse = TRUE)
reconstructed_signal <- Re(reconstructed_signal) / length(q_smoothed)

# Plot smoothed vs. reconstructed data.
plot(q_smoothed, type = "l", col = "blue", xlab = "Time", ylab = "Streamflow", main = "Smoothed vs. Reconstructed Streamflow")
lines(reconstructed_signal, col = "red")
legend("topright", legend = c("Smoothed", "Reconstructed"), col = c("blue", "red"), lty = 1)

# Section 5: Effect of Variable Smoothing Window Sizes

# Plot original data and smoothed versions with different window sizes.
plot(q.timeseries$date, q.timeseries$q, type = "l", xlab = "Date", ylab = "Streamflow", main = "Effect of Different Window Sizes")

# Define an array of window sizes to explore smoothing effects.
window_sizes <- c(30, 90, 180)  # Example: 30, 90, and 180 days.
colors <- c("red", "blue", "green")

# Loop through different window sizes and plot.
for (i in seq_along(window_sizes)) {
  filter_weights <- rep(1 / window_sizes[i], window_sizes[i])
  q_smoothed <- stats::filter(q.timeseries$q, filter_weights, sides = 2)
  q_smoothed <- na.omit(q_smoothed)  # Remove NA values.
  
  # Plot only if lengths match.
  if (length(q_smoothed) == length(q.timeseries$date)) {
    lines(q.timeseries$date, q_smoothed, col = colors[i])
  } else {# Add a legend to the plot.
legend("topright", legend = paste(window_sizes, "days"), col = colors, lty = 1)

# Section 6: Extract and Analyze Low-Frequency Components

# Identify and retain low-frequency components (e.g., below a threshold).
low_freq_threshold <- 0.01  # Define a low-frequency threshold.
fft_filtered <- ifelse(frequencies < low_freq_threshold, fft_result, 0)

# Perform inverse FFT on filtered data to reconstruct low-frequency signal.
filtered_signal <- fft(fft_filtered, inverse = TRUE)
filtered_signal <- Re(filtered_signal) / length(q)

# Plot original and low-frequency signals for comparison.
plot(q, type = "l", col = "blue", main = "Original vs. Low-Frequency Signal", ylab = "Streamflow", xlab = "Time")
lines(filtered_signal, col = "red")
legend("topright", legend = c("Original", "Low-Frequency"), col = c("blue", "red"), lty = 1)

# Plot the low-frequency FFT components.
    warning("Length of smoothed data and date data do not match for window size", window_sizes[i])
  }
}

# Add a legend to the plot.
legend("topright", legend = paste(window_sizes, "days"), col = colors, lty = 1)

# Section 6: Extract and Analyze Low-Frequency Components

# Identify and retain low-frequency components (e.g., below a threshold).
low_freq_threshold <- 0.01  # Define a low-frequency threshold.
fft_filtered <- ifelse(frequencies < low_freq_threshold, fft_result, 0)

# Perform inverse FFT on filtered data to reconstruct low-frequency signal.
filtered_signal <- fft(fft_filtered, inverse = TRUE)
filtered_signal <- Re(filtered_signal) / length(q)

# Plot original and low-frequency signals for comparison.
plot(q, type = "l", col = "blue", main = "Original vs. Low-Frequency Signal", ylab = "Streamflow", xlab = "Time")
lines(filtered_signal, col = "red")
legend("topright", legend = c("Original", "Low-Frequency"), col = c("blue", "red"), lty = 1)

# Plot the low-frequency FFT components.
plot(frequencies[frequencies < low_freq_threshold], fft_magnitude[frequencies < low_freq_threshold], type = "l", main = "Low-Frequency FFT Components", xlab = "Frequency", ylab = "Magnitude")

# Section 7: Statistical Analysis

# Compute basic statistics for original and filtered data.
mean_original <- mean(q)
mean_filtered <- mean(filtered_signal)
var_original <- var(q)
var_filtered <- var(filtered_signal)

# Print the results.
cat("Mean (Original):", mean_original, "\n")
cat("Mean (Filtered):", mean_filtered, "\n")
cat("Variance (Original):", var_original, "\n")
cat("Variance (Filtered):", var_filtered, "\n")

# Compute and plot autocorrelation functions.
acf(q, main = "ACF of Original Streamflow", plot = TRUE)
acf(filtered_signal, main = "ACF of Filtered Streamflow", plot = TRUE)
