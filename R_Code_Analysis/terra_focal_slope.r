library(terra)
library(Rcpp)

# Function to create weight matrix for any window size - edges and center only
create_edges_center_matrix <- function(window_size) {
  if (window_size %% 2 == 0) {
    stop("Window size must be odd")
  }
  
  # Create matrix of zeros
  weights <- matrix(0, nrow = window_size, ncol = window_size)
  
  # Set edges to 1
  weights[1, ] <- 1                    # Top edge
  weights[window_size, ] <- 1          # Bottom edge
  weights[, 1] <- 1                    # Left edge
  weights[, window_size] <- 1          # Right edge
  
  # Set center to 1
  center <- ceiling(window_size / 2)
  weights[center, center] <- 1
  
  return(weights)
}

# Function to create weight matrix for corners and center only
create_corners_center_matrix <- function(window_size) {
  if (window_size %% 2 == 0) {
    stop("Window size must be odd")
  }
  
  # Create matrix of zeros
  weights <- matrix(0, nrow = window_size, ncol = window_size)
  
  # Set corners to 1
  weights[1, 1] <- 1                           # Top-left
  weights[1, window_size] <- 1                 # Top-right
  weights[window_size, 1] <- 1                 # Bottom-left
  weights[window_size, window_size] <- 1       # Bottom-right
  
  # Set center to 1
  center <- ceiling(window_size / 2)
  weights[center, center] <- 1
  
  return(weights)
}

# Advanced slope calculation function for any window size
calculate_slope_variable_window <- function(dem_raster, window_size = 3, 
                                          method = "edges_center", 
                                          algorithm = "horn") {
  
  if (window_size %% 2 == 0) {
    stop("Window size must be odd")
  }
  
  # Create appropriate weight matrix
  if (method == "edges_center") {
    weights <- create_edges_center_matrix(window_size)
  } else if (method == "corners_center") {
    weights <- create_corners_center_matrix(window_size)
  } else {
    stop("Method must be 'edges_center' or 'corners_center'")
  }
  
  cell_size <- res(dem_raster)[1]
  center_idx <- ceiling(window_size^2 / 2)
  
  # Define the focal function based on algorithm
  if (algorithm == "horn") {
    focal_fun <- function(x) {
      if (sum(!is.na(x)) < 5) return(NA)
      
      # Reshape to matrix for easier indexing
      x_matrix <- matrix(x, nrow = window_size, ncol = window_size)
      
      # Calculate gradients using Horn's algorithm adapted for edges
      # Use weighted differences across the window
      dzdx <- 0
      dzdy <- 0
      weight_sum <- 0
      
      # For each edge pixel, calculate its contribution
      for (i in 1:window_size) {
        for (j in 1:window_size) {
          if (weights[i, j] == 1 && !is.na(x_matrix[i, j])) {
            # Distance from center
            di <- i - ceiling(window_size / 2)
            dj <- j - ceiling(window_size / 2)
            
            if (di != 0 || dj != 0) {  # Skip center pixel for gradient calc
              # Weight by inverse distance
              weight <- 1 / sqrt(di^2 + dj^2)
              dzdx <- dzdx + (x_matrix[i, j] * dj * weight)
              dzdy <- dzdy + (x_matrix[i, j] * di * weight)
              weight_sum <- weight_sum + weight
            }
          }
        }
      }
      
      if (weight_sum > 0) {
        dzdx <- dzdx / (weight_sum * cell_size)
        dzdy <- dzdy / (weight_sum * cell_size)
        
        slope_rad <- atan(sqrt(dzdx^2 + dzdy^2))
        return(slope_rad * 180 / pi)
      } else {
        return(NA)
      }
    }
    
  } else if (algorithm == "simple_gradient") {
    focal_fun <- function(x) {
      if (sum(!is.na(x)) < 5) return(NA)
      
      # Reshape to matrix
      x_matrix <- matrix(x, nrow = window_size, ncol = window_size)
      center <- ceiling(window_size / 2)
      
      # Simple gradient using opposite edges
      # East-West gradient
      east_vals <- x_matrix[, window_size][!is.na(x_matrix[, window_size])]
      west_vals <- x_matrix[, 1][!is.na(x_matrix[, 1])]
      
      # North-South gradient  
      north_vals <- x_matrix[1, ][!is.na(x_matrix[1, ])]
      south_vals <- x_matrix[window_size, ][!is.na(x_matrix[window_size, ])]
      
      if (length(east_vals) > 0 && length(west_vals) > 0 && 
          length(north_vals) > 0 && length(south_vals) > 0) {
        
        dzdx <- (mean(east_vals) - mean(west_vals)) / (2 * (window_size - 1) * cell_size)
        dzdy <- (mean(south_vals) - mean(north_vals)) / (2 * (window_size - 1) * cell_size)
        
        slope_rad <- atan(sqrt(dzdx^2 + dzdy^2))
        return(slope_rad * 180 / pi)
      } else {
        return(NA)
      }
    }
    
  } else {
    stop("Algorithm must be 'horn' or 'simple_gradient'")
  }
  
  # Apply focal function
  result <- focal(dem_raster, w = weights, fun = focal_fun)
  return(result)
}

# Function to visualize weight patterns
visualize_weight_pattern <- function(window_size, method = "edges_center") {
  if (method == "edges_center") {
    weights <- create_edges_center_matrix(window_size)
    title <- paste0("Edges + Center (", window_size, "x", window_size, ")")
  } else {
    weights <- create_corners_center_matrix(window_size)
    title <- paste0("Corners + Center (", window_size, "x", window_size, ")")
  }
  
  # Convert to raster for plotting
  weight_raster <- rast(weights)
  plot(weight_raster, main = title, axes = FALSE, legend = FALSE)
  
  # Add grid lines
  for (i in 0.5:(window_size + 0.5)) {
    abline(h = i, col = "gray", lwd = 0.5)
    abline(v = i, col = "gray", lwd = 0.5)
  }
  
  return(weights)
}

# Create example DEM data
set.seed(123)
dem_matrix <- matrix(runif(400, 100, 200), nrow = 20, ncol = 20)
# Add realistic terrain structure
for(i in 1:20) {
  for(j in 1:20) {
    dem_matrix[i,j] <- dem_matrix[i,j] + 
      10 * sin(i/3) * cos(j/3) +  # Rolling hills
      5 * (i-10)^2/100 +          # Slope
      3 * (j-10)^2/100            # Cross slope
  }
}

# Create SpatRaster
dem <- rast(dem_matrix)
# ext(dem) <- c(0, 2000, 0, 2000)  # 2000m x 2000m extent
# res(dem) <- 100  # 100m resolution

print("Original DEM:")
plot(dem, main = "Digital Elevation Model")

# Test different window sizes and methods
window_sizes <- c(3, 5, 7, 9, 15)
methods <- c("edges_center", "corners_center")

# Visualize weight patterns
par(mfrow = c(2, 3))
for (size in c(3, 5, 7)) {
  visualize_weight_pattern(size, "edges_center")
}
for (size in c(3, 5, 7)) {
  visualize_weight_pattern(size, "corners_center")
}

# Calculate slopes with different window sizes
slope_results <- list()

for (size in window_sizes) {
  for (method in methods) {
    key <- paste0("w", size, "_", method)
    
    print(paste("Calculating slope with window size", size, "and method", method))
    
    slope_results[[key]] <- calculate_slope_variable_window(
      dem, 
      window_size = size, 
      method = method,
      algorithm = "simple_gradient"
    )
  }
}

# Compare results visually
par(mfrow = c(2, 3))
plot(dem, main = "Original DEM")
plot(slope_results[["w3_edges_center"]], main = "3x3 Edges+Center")
plot(slope_results[["w7_edges_center"]], main = "7x7 Edges+Center")
plot(slope_results[["w15_edges_center"]], main = "15x15 Edges+Center")
plot(slope_results[["w7_corners_center"]], main = "7x7 Corners+Center")
plot(terrain(dem, "slope", neighbors = 8), main = "Built-in Slope")

# Advanced: Create custom patterns
create_custom_pattern <- function(window_size, pattern_type = "cross") {
  weights <- matrix(0, nrow = window_size, ncol = window_size)
  center <- ceiling(window_size / 2)
  
  if (pattern_type == "cross") {
    # Cross pattern: edges in cardinal directions only
    weights[center, ] <- 1  # Horizontal line
    weights[, center] <- 1  # Vertical line
    
  } else if (pattern_type == "diagonals") {
    # Diagonal pattern: corners and diagonals
    for (i in 1:window_size) {
      weights[i, i] <- 1                    # Main diagonal
      weights[i, window_size + 1 - i] <- 1  # Anti-diagonal
    }
    
  } else if (pattern_type == "ring") {
    # Ring pattern: only the outermost ring
    weights[1, ] <- 1                    # Top edge
    weights[window_size, ] <- 1          # Bottom edge
    weights[, 1] <- 1                    # Left edge
    weights[, window_size] <- 1          # Right edge
    
  } else if (pattern_type == "sparse_ring") {
    # Sparse ring: every other pixel on the edge
    for (i in seq(1, window_size, by = 2)) {
      weights[1, i] <- 1                 # Top edge
      weights[window_size, i] <- 1       # Bottom edge
      weights[i, 1] <- 1                 # Left edge
      weights[i, window_size] <- 1       # Right edge
    }
  }
  
  # Always include center
  weights[center, center] <- 1
  
  return(weights)
}

# Test custom patterns
custom_patterns <- c("cross", "diagonals", "ring", "sparse_ring")
par(mfrow = c(2, 2))

for (pattern in custom_patterns) {
  weights <- create_custom_pattern(7, pattern)
  weight_raster <- rast(weights)
  plot(weight_raster, main = paste("7x7", pattern, "pattern"), 
       axes = FALSE, legend = FALSE)
}

# Function to analyze pattern effectiveness
analyze_pattern_effectiveness <- function(dem_raster, window_sizes, methods) {
  results <- data.frame(
    window_size = integer(),
    method = character(),
    mean_slope = numeric(),
    sd_slope = numeric(),
    non_na_pixels = integer()
  )
  
  for (size in window_sizes) {
    for (method in methods) {
      slope_raster <- calculate_slope_variable_window(
        dem_raster, 
        window_size = size, 
        method = method,
        algorithm = "simple_gradient"
      )
      
      slope_values <- values(slope_raster)
      slope_values <- slope_values[!is.na(slope_values)]
      
      results <- rbind(results, data.frame(
        window_size = size,
        method = method,
        mean_slope = mean(slope_values),
        sd_slope = sd(slope_values),
        non_na_pixels = length(slope_values)
      ))
    }
  }
  
  return(results)
}

# Analyze effectiveness
effectiveness <- analyze_pattern_effectiveness(dem, c(3, 5, 7, 9, 15), methods)
print("Pattern effectiveness analysis:")
print(effectiveness)

# Performance comparison
print("\nPerformance comparison:")
system.time({
  slope_3x3 <- calculate_slope_variable_window(dem, 3, "edges_center")
})

system.time({
  slope_15x15 <- calculate_slope_variable_window(dem, 15, "edges_center")
})

system.time({
  slope_builtin <- terrain(dem, "slope")
})
