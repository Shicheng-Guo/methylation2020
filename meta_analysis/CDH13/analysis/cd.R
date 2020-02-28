

install.packages("metamisc")
install.packages("rjags")
library("metamisc")

data("Kertai")
fit2 <- riley(Kertai,type="test.accuracy")
fit2
summary(fit2)
