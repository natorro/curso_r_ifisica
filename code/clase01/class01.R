getwd()

setwd()

ls()

source()

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

