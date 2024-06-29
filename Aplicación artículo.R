library(ggplot2)
library(gridExtra)

# Datos de la tabla en R
diseno_pb <- data.frame(
  Aceite = c(1, 1, 1, 1, -1, -1, 1, -1, +1, -1, -1, -1),
  Glucosa = c(-1, -1, -1, 1, -1, 1, 1, -1, 1, 1, -1, +1),
  Mezcla1 = c(-1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1),
  Mezcla2 = c(1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, 1),
  Mezcla3 = c(-1, 1, 1, -1, -1, 1, -1, 1, 1, 1, -1, -1),
  Extrac_lev = c(1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1),
  Mezcla4 = c(1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1),
  Mezcla5 = c(1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, 1),
  Tam_inoculo = c(-1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1),
  pH = c(-1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1),
  Vol_cultivo = c(-1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1),
  Consum_aceite = c(9.8,2.98,16.87,12.78,3.41,2.56,19.18,0.51,12.78,1.7,3.58,1.02)
)

lm1 <- lm(Consum_aceite ~ ., data = diseno_pb)
lm <- aov(Consum_aceite ~ ., data = diseno_pb)
summary(lm);summary(lm1)

# Efectos y efectos principales
promedios_positivos <- sapply(diseno_pb[, -ncol(diseno_pb)], function(x) mean(diseno_pb$Consum_aceite[x == 1]))
promedios_negativos <- sapply(diseno_pb[, -ncol(diseno_pb)], function(x) mean(diseno_pb$Consum_aceite[x == -1]))

promedios <- data.frame(
  Factor = names(promedios_positivos),
  Promedio_positivo = promedios_positivos,
  Promedio_negativo = promedios_negativos,
  Efecto = (promedios_positivos-mean(diseno_pb$Consum_aceite))-(promedios_negativos-mean(diseno_pb$Consum_aceite))
)

# Crear una función para generar el gráfico de un factor con título
generar_grafico <- function(factor, promedios) {
  df <- promedios[promedios$Factor == factor, ]
  
  df_plot <- data.frame(
    Nivel = c("Positivo", "Negativo"),
    Promedio = c(df$Promedio_positivo, df$Promedio_negativo)
  )
  
  p <- ggplot(df_plot, aes(x = Nivel, y = Promedio)) +
    geom_point(size = 3) +
    geom_line(aes(group = 1)) +
    labs(title = factor,
         x = "Nivel",
         y = "Prom_consumo") +
    theme_minimal() +
    ylim(2, 13)  # Ajustar los límites del eje y
  
  return(p)
}

ggplot(promedios, aes(x = reorder(Factor, Efecto), y = Efecto)) +
  geom_point(size = 4) +
  geom_segment(aes(x = Factor, xend = Factor, y = 0, yend = Efecto), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Efectos de los Factores", x = "Factores", y = "Efecto del Factor") +
  theme_minimal() +
  coord_flip()

# Generar los gráficos para todos los factores
graficos <- lapply(rownames(promedios), function(factor) generar_grafico(factor, promedios))

# Mostrar los gráficos en una disposición de varios paneles
do.call("grid.arrange", c(graficos, ncol = 3))


## Diseño Box Behnken

library(rsm)

diseño_bb<- data.frame(
  Glucosa =c(1,1,-1,-1,1,1,-1,-1,0,0,0,0,0),
  Aceite =c(1,-1,1,-1,0,0,0,0,1,1,-1,-1,0),
  ph =c(0,0,0,0,1,-1,1,-1,1,-1,1,-1,0),
  Consum_aceite =c(54.9,41.8,18.8,12.5,0.7,9.9,25.0,18.0,22.6,40.0,71.0,25.0,22.0)
)

model<- rsm(Consum_aceite~SO(Glucosa, Aceite, ph), data = diseño_bb)
summary(model)


### Gráficos


par(mfrow = c(1, 3))

contour(model, ~ Glucosa + Aceite, image = TRUE, at = list(ph = 0))
contour(model, ~ Glucosa + ph, image = TRUE, at = list(Aceite = 0))
contour(model, ~ Aceite + ph, image = TRUE, at = list(Glucosa = 0))


# Crear una cuadrícula de valores para las variables independientes
Glucosa_seq <- seq(-1, 1, length.out = 30)
Aceite_seq <- seq(-1, 1, length.out = 30)
ph_seq <- seq(-1, 1, length.out = 30)

# Crear una cuadrícula para Glucosa y Aceite a un valor fijo de ph
grid_glucosa_aceite <- expand.grid(Glucosa = Glucosa_seq, Aceite = Aceite_seq, ph = 0)
pred_glucosa_aceite <- predict(model, newdata = grid_glucosa_aceite)

# Convertir las predicciones a una matriz
z_glucosa_aceite <- matrix(pred_glucosa_aceite, nrow = length(Glucosa_seq), ncol = length(Aceite_seq))

# Crear una cuadrícula para Glucosa y ph a un valor fijo de Aceite
grid_glucosa_ph <- expand.grid(Glucosa = Glucosa_seq, Aceite = 0, ph = ph_seq)
pred_glucosa_ph <- predict(model, newdata = grid_glucosa_ph)

# Convertir las predicciones a una matriz
z_glucosa_ph <- matrix(pred_glucosa_ph, nrow = length(Glucosa_seq), ncol = length(ph_seq))

# Crear una cuadrícula para Aceite y ph a un valor fijo de Glucosa
grid_aceite_ph <- expand.grid(Glucosa = 0, Aceite = Aceite_seq, ph = ph_seq)
pred_aceite_ph <- predict(model, newdata = grid_aceite_ph)

# Convertir las predicciones a una matriz
z_aceite_ph <- matrix(pred_aceite_ph, nrow = length(Aceite_seq), ncol = length(ph_seq))


# Graficar en 3D utilizando persp
par(mfrow = c(1, 3))

# Gráfico para Glucosa y Aceite
persp(Glucosa_seq, Aceite_seq, z_glucosa_aceite, theta = 30, phi = 30, col = "lightblue",
      xlab = "Glucosa", ylab = "Aceite", zlab = "Consumo de Aceite",
      main = "Glucosa y Aceite a ph = 0")

# Gráfico para Glucosa y ph
persp(Glucosa_seq, ph_seq, z_glucosa_ph, theta = 30, phi = 30, col = "lightgreen",
      xlab = "Glucosa", ylab = "ph", zlab = "Consumo de Aceite",
      main = "Glucosa y ph a Aceite = 0")

# Gráfico para Aceite y ph
persp(Aceite_seq, ph_seq, z_aceite_ph, theta = 30, phi = 30, col = "lightpink",
      xlab = "Aceite", ylab = "ph", zlab = "Consumo de Aceite",
      main = "Aceite y ph a Glucosa = 0")
