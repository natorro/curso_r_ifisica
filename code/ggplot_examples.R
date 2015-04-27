library(ggplot2)

qplot(displ, hwy, data = mpg, geom=c("point", "smooth"))
ggplot(data=mpg) + aes(x=displ, y=hwy) + geom_point() + geom_smooth()

qplot(hwy, data=mpg, fill = drv)
ggplot(data = mpg) + aes(x=hwy, fill = drv) + geom_histogram()

qplot(hwy, displ, data=mpg, facets = . ~ drv)
ggplot(data=mpg) + aes(x=hwy, y=displ, facets= . ~ drv) + geom_point()

qplot(hwy, data=mpg, facets = drv ~ ., binwidth=2)
ggplot(data=mpg) + aes(x=hwy, y=displ, facets= . ~ drv) + geom_histogram()

ggplot(data=mpg) + aes(x=hwy, fill=drv) + geom_histogram(alpha=0.5)



