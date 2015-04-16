getwd()
setwd()
ls()

# Evaluation
x <- 1
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

for(i in 1:4){
  print(x[i])
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


# Crear una función que tome dos vectores y los sume y los multiplique 
# y regrese