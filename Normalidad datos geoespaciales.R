library(geoR)
library(nortest)
library(MASS)

# Crear grillas 
grid_25 <- expand.grid(seq(0, 9, length.out = 5), seq(0, 9, length.out = 5))
grid_49 <- expand.grid(seq(0, 9, length.out = 7), seq(0, 9, length.out = 7))
grid_100 <- expand.grid(seq(0, 9, length.out = 10), seq(0, 9, length.out = 10))

set.seed(2024)
x_cluster <- rnorm(100, mean = 5, sd = 1.5)
y_cluster <- rnorm(100, mean = 5, sd = 1.5)
grid_cluster_100 <- data.frame(x = x_cluster, y = y_cluster)

# Funciones necesarias
campo_aleatorio <- function(grid, phi, sigma, mu) { #Campo aleatorio
  coords <- as.matrix(grid)
  n <- nrow(coords)
  grf_result <- grf(n = n, grid = coords, cov.pars = c(sigma, phi), mean = mu, messages = FALSE)
  return(grf_result$data)
}

Q_stat<- function(z, mu, cov_matrix) { #Estadística Q
  diff_z_mu <- z - mu
  term <- solve(cov_matrix) %*% diff_z_mu
  q_statistic <- t(diff_z_mu) %*% term
  n <- length(z)
  p_value <- 1 - pchisq(q_statistic, df = n)
  return(list(q_statistic = q_statistic, p_value = p_value))
}

covar <- function(grid, sigma, phi) { #Matriz var-cov
  distancia <- dist(grid, diag = TRUE, upper = TRUE)
  distancia <- as.matrix(distancia)
  covar <- cov.spatial(distancia, cov.model = "exponential", cov.pars = c(sigma, phi))
  return(covar)
}

normalidad <- function(data, mu, y_bar, covarianza) { # Pruebas de normalidad
  shapiro <- shapiro.test(data)
  ad <- ad.test(data)
  sf <- sf.test(data)
  Q1 <- Q_stat(data, mu, covarianza)
  Q2 <- Q_stat(data, y_bar, covarianza)
  return(list(shapiro = shapiro, 
              ad = ad, 
              sf = sf,
              Q1 = Q1,
              Q2 = Q2))
}

simulacion <- function(num_simul, grid, phi, sigma, mu) {
  rechazo <- matrix(0, nrow = num_simul, ncol = 5)  # Almacenar rechazos 
  
  for (i in 1:num_simul) {
    z <- campo_aleatorio(grid, phi, sigma, mu) # Generar campo aleatorio normal
    cov_matrix <- covar(grid, sigma, phi)  # Calcular la matriz de covarianza
    y_bar <- rep(mean(z), nrow(grid))
    mu_1<-rep(mu, nrow(grid))
    tests <- normalidad(z, mu_1, y_bar, cov_matrix)
    
    # Evaluar rechazos
    rechazo[i, 1] <- as.numeric(tests$shapiro$p.value < 0.05)
    rechazo[i, 2] <- as.numeric(tests$ad$p.value < 0.05)
    rechazo[i, 3] <- as.numeric(tests$sf$p.value < 0.05)
    rechazo[i, 4] <- as.numeric(tests$Q1$p_value < 0.05)
    rechazo[i, 5] <- as.numeric(tests$Q2$p_value < 0.05)
  }
  
  # Tasa de rechazo promedio por cada prueba
  tasa_rechazo <- colMeans(rechazo)
  return(tasa_rechazo)
}

# Configuración de parámetros
phi <- 3
sigma <- 8
mu <- 100
num_simul <- 1000

# Realizar simulación de Monte Carlo para cada configuración de grilla
tasa_rechazos_25 <- simulacion(num_simul, grid_25, phi, sigma, mu)
tasa_rechazos_49 <- simulacion(num_simul, grid_49, phi, sigma, mu)
tasa_rechazos_100 <- simulacion(num_simul, grid_100, phi, sigma, mu)
tasa_rechazos_cluster_100 <- simulacion(num_simul, grid_cluster_100, phi, sigma, mu)

# Mostrar resultados
tasa_rechazos_25
tasa_rechazos_49
tasa_rechazos_100
tasa_rechazos_cluster_100
