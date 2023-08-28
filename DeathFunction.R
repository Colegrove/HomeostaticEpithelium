# Title     : Shows how death function scales death with EGF concentration
# Objective : plots that display probabilies of death with different exponents
# Created by: schencro
# Created on: 12/15/17

t <- seq(0,0.008,0.00001)

eq <- function(x){
    return((1-x/0.005939094)^2)
}

par(mfrow=c(1,3))
out <- sapply(t, eq)
ex <- paste("f([E]) = (1-[E]/α)","^2",sep="")
plot(t,out,xlab = "[E]",ylab="f([E])", type="l",main=ex)
abline(v=0.005939094, col="red", lty=2)
abline(v=0, col="green", lty=2)
abline(h=0, col="green", lty=2)

eq <- function(x){
    return((1-x/0.005939094)^3)
}

out <- sapply(t, eq)
ex <- paste("f([E]) = (1-[E]/α)","^3",sep="")
plot(t,out,xlab = "[E]",ylab="f([E])", type="l",main=ex)
abline(v=0.005939094, col="red", lty=2)
abline(v=0, col="green", lty=2)
abline(h=0, col="green", lty=2)

eq <- function(x){
    return((1-x/0.005939094)^5)
}

out <- sapply(t, eq)
ex <- paste("f([E]) = (1-[E]/α)","^5",sep="")
plot(t,out,xlab = "[E]",ylab="f([E])", type="l",main=ex)
abline(v=0.005939094, col="red", lty=2)
abline(v=0, col="green", lty=2)
abline(h=0, col="green", lty=2)

