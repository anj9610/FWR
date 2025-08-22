FWR code 
library(sp)
library(spdep)
library(MASS)

apartment_data <- read.csv("C:/Users/USER/Desktop/GIS/개인연구/부동산/아파트데이터/아파트 데이터 기본(7).csv")

migration_data <- w2  


region_names <- rownames(migration_data)


apartment_data <- apartment_data[match(region_names, apartment_data$ADM_DR_NM), ]


unique_names <- region_names


W_f <- matrix(NA, nrow = length(unique_names), ncol = length(unique_names), 
              dimnames = list(unique_names, unique_names))

for (i in 1:ncol(migration_data)) {  
  total_inflow <- sum(migration_data[, i])  
  
  for (j in 1:nrow(migration_data)) {  
    D_ij <- migration_data[j, i]  
    W_f[i, j] <- ifelse(total_inflow > 0, D_ij / total_inflow, 0)
  }
}



ranked_regions <- apply(W_f, 1, function(x) order(x, decreasing = TRUE))


local_coefficients <- list()


local_coefficients <- list()

for (i in 1:length(region_names)) {
  region_index <- i
  
  
  selected_regions <- unique(c(region_index, ranked_regions[1:33, region_index])) 
  
  
  selected_migration <- W_f[region_index, selected_regions]
  
  
  selected_migration[selected_regions == region_index] <- max(selected_migration)
  
  
  W_i <- diag(selected_migration)
  
  
  X <- as.matrix(apartment_data[selected_regions, c("meanYR", "elite", "univ_rate", "high_dist", "low_eld", "subway")])
  y <- apartment_data[selected_regions, "Price"]
  
  
  if (nrow(X) >= ncol(X) && all(!is.na(X)) && all(!is.na(y))) {
    XtW <- t(X) %*% W_i
    XtWX <- XtW %*% X + lambda * diag(ncol(X))  # Ridge 정규화
    beta_i <- tryCatch(solve(XtWX) %*% XtW %*% y, error = function(e) rep(NA, ncol(X)))
    
    
    local_coefficients[[i]] <- beta_i
  } else {
    local_coefficients[[i]] <- rep(NA, ncol(X))
  }
}