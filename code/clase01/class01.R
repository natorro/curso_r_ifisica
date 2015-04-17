getwd()
setwd()
ls()

# Evaluation
x <- 2
print(x)
x
msg <- "hello"

# Grammar
x <- ;

# Evaluation
x <- 5 
x
print(x)

# Printing
x <- 1:20
x

# Creating vectors
x <- c(0.5, 0.6)
x <- c(TRUE, FALSE)
x <- c(T, F)
x <- c("a", "b", "c", "d")
x <- 9:29
x <- c(1+0i, 2+4i)
x <- vector("numeric", length = 10)

# Mixing objects
x <- c(1.7, "a")
y <- c(TRUE, 2)
y <- c("a", TRUE)

# Explicit Coercing 

x <- 0:6
class(x)
as.numeric(x)
as.logical(x)
as.character(x)
as.complex(x)

x <- c("a", "b", "c")
as.numeric(x)
as.logical(x)
as.complex(x)

# Matrices: Vectors with dimension attributes

m <- matrix(nrow = 2, ncol = 3)
m
dim(m)
attributes(m)

m <- matrix(1:6, nrow = 2, ncol = 3)
m

# cbind, rbind

x <- 1:3
y <- 10:12

cbind(x, y)
rbind(x, y)

# lists
x <- list(1, "a", TRUE, 1+4i)
x

# factor, special for lm() glm()
# ordered or not ordered

x <- factor(c("yes", "yes", "no", "yes", "no"))
x
table(x)
unclass(x)
x <- factor(c("yes", "yes", "no", "yes", "no"), levels = c("yes", "no"))
x

# missing values
is.na()
is.nan()

x <- c(1, 2, NA, 10, 3)
is.na(x)
is.nan(x)
x <- c(1, 2, NaN, NA, 4)
is.na(x)
is.nan(x)

# Dataframes

read.table()
read.csv()
data.matrix()

x <- data.frame(coluno = 1:4, col2 = c(T, T, F, F))
x
nrow(x)
ncol(x)

# Names

x <- 1:3
names(x)
names(x) <- c("uno", "dos", "tres")
x
names(x)

x <- list(a = 1, b = 2, c = 3)
x

m <- matrix(1:4, nrow = 2, ncol = 2)
m
dimnames(m) <- list(c("a", "b"), c("c", "d"))
m

# Subsetting

x <- c("a", "b", "c", "c", "d", "a")
x[1]
x[2]
x[1:4]
u <- x > "a"
u
x[u]

# Subsetting a matrix

x <- matrix(1:6, 2, 3)
x
x[1, 2]
x[2, 1]

x[1, ]
x[, 2]
x[1, 2, drop=FALSE]

x <- matrix(1:6, 2, 3)
x[1, ]
x[1, , drop=FALSE]

# Subsetting lists
x <- list(uno = 1:4, dos = 0.6)
x
x[1]
x[[1]]
x$uno
x$dos
x[["uno"]]
x["uno"]
x
x <- list(uno = 1:4, dos = 0.6, msg = "hello")
x
x[c(1, 3)]

name <- "uno"
x[[name]]
x$name
x$uno

x <- list(a = list(10, 12, 14), b= c(2.4, 3.8))
x
x[[c(1, 3)]]
x[[1]][[3]]
x[[c(2, 1)]]

# Partial Matching
x <- list(aardvardkss = 1:5)
x
x$a
x[["a"]]
x[["a", exact=FALSE]]

# Removing NA values
x <- c(1, 2, NA, 4, NA, 5)
bad <- is.na(x)
x[!bad]

x <- c(1, 2, NA, 4, NA, 5)
y <- c("a", "b", NA, "d", NA, "f")

good <- complete.cases(x, y)
good
x[good]
y[good]

airquality
airquality[1:6, ]
good <- complete.cases(airquality)
airquality[good, ][1:6, ]
# Read and write data
# read.table(), read.csv(), readLines(), source(), dget(), load(), unserialize()
# write.table(), writeLines(), dump(), dput(), save(), serialize()
data <- read.table("myfile.txt")

# Using dput to save an r object
x <- data.frame(a = 1, b = "a")
x
dput(x)
dput(x, "dputdata.R")
recoveredData <- dget("dputdata.R")
recoveredData

# Using dump to save an r object
x <- "somechar"
y <- data.frame(a = 1, b = "a")
dump(c("x", "y"), file = "datadump.R")

rm(x, y)
rm("recoveredData")
ls()
x
y
source("datadump.R")
x
y

# Connections
conection <- file("hola.txt", "r")
data <- read.csv(conection)
close(con)
data <- read.csv("hola.txt")
# Conection to a gz file
conection2 <- gzfile("words.gz")
x <- readLines(conection2, 10)
x
# Conection to a URL
conection3 <- url("http://www.jhsph.edu", "r")
conection3
x <- readLines(conection3)
head(x)
summary(x)
x[24]

# Control structures
# if, else
# for 
# while
# repeat
# break
# next
# return
#if(condition) {
#  do something
#} else {
#  do something else
#}

if(x > 3) {
  y <- 10
} else {
  y <- 0
}

y <- if(x > 3){
  10
} else {
  0
}

for(i in 1:10){
  print(i)
}

x <- c("a", "b", "c", "d")

for(i in x){
  print(i)
}

for(i in seq_along(x)){
  print(x[i])
}

for(letter in x){
  print(letter)
}

for(i in 1:4) print(x[i])

# for loops can be nested

x <- matrix(1:6, 2, 3)
x
for(i in seq_len(nrow(x))){
  for(j in seq_len(ncol(x))){
    print(x[i, j])
  }
}

#while loop
count <- 0
while(count < 10){
  print(count)
  count <- count + 1
}

# when will this stop? 

z <- 5 
while(z >= 3 && z <= 10){
  print(z)
  coin <- rbinom(1, 1, 0.5)
  
  if(coin == 1) {
    z <- z+1
  } else {
    z <- z - 1
  }
}

#repeat
x <- 1
tol <- 1e-8

repeat{
  x1 <- computeEstimate()
  
  if(abs(x1 - x) < tol){
    break
  } else {
    x <- x1
  }
}

# next - return
for (i in 1:100) {
  if (i <= 20) {
    next
  }
  #Do something
  return()
}

#functions, 
# 1. can be passed as arguments to other functions
# 2. can be nested
# 3. named arguments

nameOfFunction <- function(<arguments>){
  #do something here
}

alreves <- function(x){
  if(TRUE == is.matrix(x)){
    print("Error: x is not a matrix")
  } else{
    size <- length(x)
    y <- vector()
    for(i in 1:size){
      j <- i - 1
      y[i] <- x[size-j]
    }
    y
  }
}

# argument matching
misdatos <- rnorm(100)
var(data)
var(x = data)
var(x = data, na.rm = FALSE)
var(na.rm = FALSE, x = data )
var(na.rm = FALSE, data)

# Definition
f <- function(a, b = 1, c = 2, d = NULL){
  
}

# Lazy evaluation
g <- function(a, b){
  a ^ 2
}
g(2)
g(2, 4)

g <- function(a, b){
  print(a)
  print(b)
}

g(5)

# ... argument

mygraphic <- function(x, y, type="l", ...) {
  plot(x, y, type = type, ...)
}

paste("hol", "a", sep=" ")
paste("hol", "a", se=" ")

# Create a function that takes two vectors and sum them together and 
# takes the dot product and then returns both results in a list.


# Create a function that takes as an argument a matrix and gives back the
# mean o every row in the form of a vector

# Crear una función que tome dos vectores y los sume y los multiplique 
# y regrese ambos vectores en forma de lista

# Search lists

# Functions defined inside other functions:

toThePower <- function(n){
  pow <- function(x){
    x^n
  }
  pow
}

cube <- toThePower(3)
cube(2)
squared <- toThePower(2)
squared(2)

ls(environment(cube))
ls(environment(squared))

get("n", environment(cube))
get("n", environment(squared))

y <- 10

f <- function(x) {
  y <- 2
  y^2 + g(x)
}

g <- function(x){
  x * y
}

f(3)

# vectorized operations
x <- 1:4
y <- 10:13
x > 2
x >= 2
y == 10
x * y
x / y

# matrix multiplication
A <- matrix(1:16, 4, 4)
B <- matrix(21:36, 4, 4)
B
A %*% B

# dates and times
# Date Class

as.Date("2015-01-01")-as.Date("1977-12-25")
unclass(as.Date("2015-01-01"))-unclass(as.Date("1977-12-25"))
weekdays()
months()

as.Date("2015-01-01"")
Sys.time()
x <- as.Date(Sys.time(),format = "d%M%Y%")
names(unclass(as.POSIXlt(x)))
as.POSIXlt(x)$sec

datestring <- c("January 10, 2012 10:40", "December 9, 2011 9:12")
datestring
class(strptime(datestring, "%B %d, %Y %H:%M"))

fix(function)

# Graficación basica

x <- c(12, 15, 13, 20, 14, 16, 10, 10, 8, 15)
hist(x)
x <- seq(1, 10)
y <- x^2 - 10*x
plot(x, y)
curve(expr = sin, from = 0, to = 10)
curve(x^2 -10 * x, from = 1, to = 10)
#Funciones basicas estadísticas
median(x)
var(x)
summary(x)


# Graficos
# Usaremos el dataset VADeaths
VADeaths
        
        # Haremos nuestro primer grafica de barras
        
        barplot(VADeaths, beside=TRUE, legend=TRUE,
        ylim = c(0, 90), ylab="Muertes c/1000",
        main = "Tasa de muertes en Virginia")
        
        dotchart(VADeaths, xlim=c(0, 75),
        xlab = "Muertes c/1000",
        main = "Tasa de muertes en Virginia")
        
        # Graficas de pay
        groupsizes <- c(18, 30, 32, 10, 10)
        labels <- c("A", "B", "C", "D", "F")
        
        pie(groupsizes, labels, col=c("grey40", "white",
        "red", "grey",
        "grey90"))
        
        # Ejemplo que hicimos en clase para generar una
        # grafica de pastel usando los colores del sistema
        # (los primeros veinte)
        
        colores <- colors()
        mis_valores <- rep(1, 20)
        pie(mis_valores, col=colores[1:20])
        
        # Histogramas
        # Los intervalos son log_2(n) + 1 esto se conoce
        # Como la regla de Sturges
        
        x <- rnorm(100)
        hist(x)
hist(x, breaks = "Sturges")
hist(x, breaks = "Scott")
hist(x, breaks = "Freedman-Diaconis", main="HOLA MUNDO")
        
        # Boxplots
        
        boxplot(Sepal.Length ~ Species, data = iris,
        ylab = "Longitud de sepalo (cm)",
        main = "Medidas de Iris",
        boxwex = 0.5)
        
        # Graficas de dispersion
        x <- rnorm(100)
        y <- rpois(100, 30)
        mean(y)
        
        plot(x, y, main = "Poisson vs. Normal")
        plot(x, y, type="l")
        plot(sort(x), sort(y), type = "l")
        
        # qqplots
        
        X <- rnorm(1000)
        
        A <- rnorm(1000)
        qqplot(X, A, main = "A y X son iguales")
        
        B <- rnorm(1000, mean = 3, sd = 2)
        qqplot(X, B, main = "B es X re-escalada")
        
        C <- rt(1000, df=2)
        qqplot(X, C, main = "C tiene colas pesadas")
        
        D <- exp(rnorm(1000))
        qqplot(X, D, main = "D esta sesgada a la derecha")
# Funciones de graficacion de bajo nivel

# points(x, y, ...)
        # lines(x, y, ...)
        # text(x, y, labels, ...)
        # abline(a, b, ...) linea y = a + bx
        # abline(h=y, ...) linea horizontal
        # abline(v=x, ...) linea vertical
        # polygon(x, y, ...) añade un poligono cerrado
        # segments(x0, y0, x1, y1, ..) dibuja segmentos de linea
        # arrows(x0, y0, x1, y1, ...) dibuja flechas
        # symbols(x, y, ... ) dibuja circulos, cuadrados, termometros, etc
        # legend(x, y, legend, ...) dibuja una leyenda
        
        # Ejemplo
        
        sex <- c("M", "F", "M", "F", "F", "M", "F", "M")
        length.index <- c(7.9, 6.5, 8.4, 5.5, 6.5, 8.0, 7.0, 7.5)
        width <- c(2.3, 1.7, 2.6, 1.7, 1.9, 2.1, 1.8, 1.9)
        indexfinger <- data.frame(sex, length.index, width)
        indexfinger
        
        #plot(width ~ length.index, data=indexfinger)
        
        with(indexfinger[c(3, 7), ], points(length.index, width, pch=17))
        
        plot(width ~ length.index, pch=as.character(sex), data=indexfinger)
        
        abline(lm(width ~ length.index, data=indexfinger, subset=sex=="M"), lty=1)
        abline(lm(width ~ length.index, data=indexfinger, subset=sex=="F"), lty=2)
        legend("topleft", legend=c("Male", "Female"), lty=1:2)
        
        # Funciones para anotar
        
        # title(main, sub, xlab, ylab, ...) agrega un titulo
        
# mtext(text, side, line, ...) dibuja texto en los margenes
# axis(side, at, labels, ...) agrega un eje a la grafica
# box(...) agrega una caja

# Usando locator para localizar exactamente:

        mi_vector <- rnorm(10000)
        hist(mi_vector)
        text(c(-3, 1866),"nv=75")
        
        par(bg="yellow")
        
        hist(mi_vector)
        
        par(bg="white")
        # Podemos intentar usar lowess para ajustar una recta:
        x <- rnorm(100)
        y <- rnorm(100)
        plot(data.frame(x, y))
        lines(lowess(data.frame(x, y)))
        
        # Funciones explí?citas con curve
        
        g <- function(t){
        return (t^2+1)^0.5
        }
        
        x <- seq(0,5,length=10000)
        
        y <- g(x)
        
        plot(x,y,type="l")
        
        curve((x^2+1)^0.5,0,5)
        
        curve(x^2, 0, 10)
        curve((x^2+1)^0.5,0,5, add=TRUE)
        curve((x^2+1)^0.5,0,5,add=TRUE)
# optimization of code

x <- rnorm(10000)
y <- rnorm(10000)
z <- c()
system.time(
  for (i in 1:10000){
    z <- c(z, x[i] + y[i])
  }
  )

z <- rep(NA, 10000)

system.time(
  for (i in 1:10000){
    z[i] <- x[i] + y[i]
  }
  )
system.time(z <- x + y)
f <- function(x) return((x^2+1)^0.5)

plot(f,0,5)
        
        # Guardando a archivos:
        
        #pdf("c://archivo.pdf")
        
        #dev.list()
        
        #dev.set()
        
        #dev.off()
        
        # gráficas en tres dimensiones
        
        library(lattice)
        a <- 1:10
        b <- 1:15
        eg <- expand.grid(x=a,y=b)
        eg$z <- eg$x^2 + eg$x * eg$y
        wireframe(z ~ x+y, eg)
        
        #cloud:
        mis_datos <- data.frame(rnorm(100),
        rexp(100),
        rnorm(100))
        names(mis_datos) <- c("x", "y", "z")
        
        cloud(z ~ x + y, data=mis_datos)
        cloud(z ~ x + y, eg, add=TRUE)
        # Ya sabemos hacer histogramas
        datos <- rexp(10000)
        hist(datos)
        
        hist(datos, breaks=15)
hist(datos, breaks=15, freq=F)

density(datos)
        
        #lines()
        
        
        holas <- rnorm(1000)
        hist(holas)
        density(holas)
        lines(density(holas))
        
        hist(holas, freq=F)
        lines(density(holas, na.rm=T))
        
        media <- mean(holas)
        est_dev <- sd(holas)
        hist(holas, freq=F)
        curve(dnorm(x, media, est_dev), add=T)
        
        # Y podemos agregar más opciones
        
        hist(holas, breaks = 15, freq = F,
        main = "Nuestro histograma editado",
        xlab = "Satisfaccion de vida",
        #  xlim = c(5,40),
        #  ylim = c(0,.1),
        col = "grey", las = 1)
        
        curve(dnorm(x, media, est_dev), add = T,
        lty = 2, lwd = 2, col = "blue")
        
        abline(h = 0)
        # en vez de coordenadas podemos usar locator(1) para escoger
# exactamente donde queremos la gráfica.

symbols(100, 100, circles = 2, add = T,
        inches = F, fg = "red", lwd = 2)
        # Exactamente lo mismo, podemos usar locator(1) en vez de 100, 100
        
        text(100, 110, "Punto raro", col = "red")
        
        # Comparar dos gráficas
        
        ls.1 <- rnorm(300)
        
        m <- mean(ls.1, na.rm = T)
        
        sd <- sd(ls.1, na.rm = T)
        
        curve(dnorm(x, m, sd), m-3*sd, m+3*sd,lty = 2)
        
        lines(density(ls.1, na.rm = T), lwd = 2)
        
        legend("topleft", c("Función de densidad normal",
        "Densidad de Kernel"),
        lty = c(2,1), lwd = c(1,2), box.lty = 0)
        
        # QQ Norm y QQPlot
        
        qqnorm(ls.1)
        library(car)
        qqPlot(ls.1)
        
        # Una más
        
        ls.1 <- rnorm(300)
        valores <- 1:300
        
        lin <- lm(ls.1 ~ valores)
        
        plot(ls.1 ~ valores)

abline(lin, lwd = 2)
        
        x = predict(lin, interval = "confidence")
        
        x[1:5, ]
        lines(valores, x[,2], lty = 2, col = 2)
        lines(valores, x[,3], lty = 2, col = 2)
        
        # podemos crear bandas
        plot(ls.1 ~ valores)
        polygon(c(valores, rev(valores)), c(x[,2], rev(x[,3])), col = "lightgrey", border = NA)
        abline(lin, lwd = 2)
        axis(2, labels = F)
        axis(4, labels = F, lwd.ticks = 0)
        
        #Podemos poner más de una gráfica
        
        par(mfrow=c(2, 2))
        plot(lin)
        
        
# Una más

ls.1 <- rnorm(300)
        valores <- 1:300
        
        lin <- lm(ls.1 ~ valores)
        
        plot(ls.1 ~ valores)
        
        abline(lin, lwd = 2)
        
        x = predict(lin, interval = "confidence")
        
        x[1:5, ]
        lines(valores, x[,2], lty = 2, col = 2)
        lines(valores, x[,3], lty = 2, col = 2)
        
        # podemos crear bandas
        plot(ls.1 ~ valores)
        polygon(c(valores, rev(valores)), c(x[,2], rev(x[,3])), col = "lightgrey", border = NA)
        abline(lin, lwd = 2)
        axis(2, labels = F)
        axis(4, labels = F, lwd.ticks = 0)
        
        # Podemos poner más de una gráfica
        
        par(mfrow=c(2, 2))
        plot(lin)
        
# Algebra lineal con matrices
# Determinante:
det(X)
        
        # Diagonal
        diag(X)
        
        # Podemos encontrar la traza:
        
        traza <- function(datos) sum(diag(datos))
        
        # Podemos crear una matriz diagonal con diag:
        
        diag(c(1:3))
        
        # La transpuesta:
        
        t(X)
        
        # Matrices triangulares
        # Usamos las funciones lower.tri() y uppper.tri()
        # para obtener matrices
        X
        upper.tri(X)
        lower.tri(X)
        X3 <- X
        X3[upper.tri(X, diag=TRUE)] <- 0
        
        # Aritmetica de matrices
        
        # Multiplicacion por un escalar
        Y <- 2 * X
        Y

# Suma elemento a elemento

X
        Y
        X + Y
        
        # Multiplicacion de matrices e
        # inversion de una matriz
        
        # Recordemos que para que podamos multiplicar dos matrices
        # estas deben ser conformes entre ellas.
        # El operador es %*%
        
        X %*% Y
        
        # Creemos matrices que no son conformes
        
        # Una manera mas eficiente de crear t(Y) %*% X
        
        #crossprod()
        
        # Inversion de matrices
        
        A <- matrix(c(3, 1, -4, 2), nrow=2)
        
        # Usamos las funciones solve() y qr.solve()
        
        solve(A)
        solve(A) %*% A
        
        qr.solve(A) %*% A
        
        # Como resolver un sistema de ecuaciones lineales:
        
        b <- c(1, 2, 3)
        
        b
X <-1/cbind(seq(1, 3), seq(2, 4), seq(3, 5))
X
        x <- solve(X, b)
        x
        
        # Eigenvalues y eigenvectors
        
        eigen(X)
        
        # Descomposicion de valor singular
        # Tres matrices cuadradas U D V
        # D es diagonal
        # U y V son ortogonales
        # La relacion entre estas tres matrices es:
        # A = U D t(V)
        
        # Podemos encontrar estos valores de la siguiente manera:
        X.svd <- svd(X)
        
        # Podemos ver que esta matriz es correcta
        X.svd$u %*% diag(X.svd$d) %*% t(X.svd$v)
        # Podemos encontrar la inversa con esto
        X.svd$v %*% diag(1/X.svd$d) %*% t(X.svd$u)
        
        # La decomposicion Cholesky
        # Si A es una matriz positiva definida, tiene raiz cuadrada
        # De hecho hay varias matrices tal que Bcuadrada = A
        # La decomposicion Cholesky toma esa idea de tal manera que
        # t(U)U = A
        
        X.chol <- chol(X)
        X.chol
        crossprod(X.chol, X.chol)
        chol2inv(X.chol)
        
        # Para encontrar la solucion a un problema lineal
b <- seq(1, 3)
y <- forwardsolve(t(X.chol), b)
        backsolve(X.chol, y)
        
        # Descomposicion QR
        
        X.qr <- qr(X)
        X.qr
        Q <- qr.Q(X.qr)
        Q
        R <- qr.R(X.qr)
        R
        Q %*% R
        
        qr.solve(R) %*% t(Q)
        
        # Numero de condicion
        kappa(X)
        
        kappa(diag(c(1, 1, 1)))
        
        # Funciones outer()
        x1 <- seq(1, 5)
        x1
        outer(x1, x1, "/")
        outer(x1, x1, "-")
        y <- seq(5, 10)
        outer(x1, y, "+")
        
        # Funcion apply() mucho mas eficiente computacionalmente
        apply(X, 1, sum)
        X
        sum(X[3, ])
        
        
        # Graficando datos con base
library(ggplot2)
#library(fImport)
        library(tseries)
        telmex <- get.hist.quote(instrument="TELMEXL.MX", start="2010-04-15", end="2011-09-27", quote=c("O","H","L","C","A","V"), provider="yahoo", retclass="zoo")
        
        america_movil <- get.hist.quote(instrument="AMX", start="2010-04-15", end="2011-09-27", quote=c("O","H","L","C","A","V"), provider="yahoo", retclass="zoo")
        
        
        nrow(telmex)
        nrow(america_movil)
        
        telmex[1:5, ]
        america_movil[1:5, ]
        
        
        plot(telmex$Open)
        
        plot(telmex)
        
        plot(telmex$Open - telmex$Close)
        
        plot(as.numeric(telmex$Open- telmex$Close), type="line")
        
        hola <- as.numeric(telmex$Open-telmex$Close)
        x<- 1:length(hola)
        qplot(x, hola, geom="line", ylab="mi titulo")
        
        datos <- NULL
        datos[1] <- telmex$Volume[1]
        for(i in 2:length(telmex$Volume)){
        datos[i] <- as.numeric(telmex$Volume[i]) - as.numeric(telmex$Volume[i-1])
        }
        
        plot(datos, type="line", ylab="cambios en Volumen día a día")
        
        hola2 <- NULL
        for(i in 2:length(telmex$Volume)){
        hola2[i] <- as.numeric(telmex$Open[i]) - as.numeric(telmex$Close[i-1])
        }
        
        
        #par(mfrow=c(3, 1))
        plot(hola, type="line", ylab="cambios en Volumen dia x dia")
        plot(telmex$Open-telmex$Close, ylab="Open(i)-Close(i)", xlab="")
        plot(hola2, type="line", ylab="cambios cierre y apertura")
        
        plot(density(hola))
        plot(density(telmex$Open-telmex$Close))
        hola2 <- hola2[!(is.na(hola2))]
        plot(density(hola2))
        