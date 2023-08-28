library(ggplot2)
library(gridExtra)
library(reshape2)

GetLine <- function(df){

  return(fit)
}

#BenchMarking
twoD <- c(67,0.780,100,1.061,200,2.100,300,3.148,1000,10.997,2000,49.485,3000,99.625,5000,271.884,8000,644.676,10000,980.436)
threeD <- c(67,NA,100,3.581,200,16.034,300,9.581,1000,41.921,2000,97.037,3000,169.512,5000,417.149,8000,795.338,10000,1136.026)
DF <- as.data.frame(matrix(NA,nrow=0,ncol=4))
colnames(DF) <- c("Size","Time","Model")
for(i in seq(1,length(twoD),by=2)){
    print(i)
    tmp <- data.frame(Size=twoD[i], Time=twoD[i+1]/(365*2), Model="2D")
    tmp3D <- data.frame(Size=threeD[i], Time=threeD[i+1]/(365*2), Model="3D")
    DF <- rbind(DF,tmp, tmp3D)
}

fit2D <- lm(Time~Size,subset(DF,DF$Model=="2D"))
fit3D <- lm(Time~Size,subset(DF,DF$Model=="3D"))
line2D <- function(x){
  t <- tF*365
  y <- fit2D[[1]][2]*x+fit2D[[1]][1]
  return(y*t / 60 / 60)
}
line3D <- function(x){
  t <- tF*365
  y <- fit3D[[1]][2]*x+fit3D[[1]][1]
  return(y*t / 60 / 60)
}

twoDF <- subset(DF,DF$Model=="2D")
threeDF <- subset(DF,DF$Model=="3D")
tF <- 65
dfRuns <- data.frame(Size=twoDF$Size, Time=sapply(twoDF$Size, FUN="line2D"), Age="65 year sim")
dfRuns5 <- data.frame(Size=threeDF$Size, Time=sapply(threeDF$Size, FUN="line3D"), Age="65 year sim")
tF <- 55
dfRuns2 <- data.frame(Size=twoDF$Size, Time=sapply(twoDF$Size, FUN="line2D"), Age="55 year sim")
dfRuns6 <- data.frame(Size=threeDF$Size, Time=sapply(threeDF$Size, FUN="line3D"), Age="55 year sim")
tF <- 58
dfRuns3 <- data.frame(Size=twoDF$Size, Time=sapply(twoDF$Size, FUN="line2D"), Age="58 year sim")
dfRuns7 <- data.frame(Size=threeDF$Size, Time=sapply(threeDF$Size, FUN="line3D"), Age="58 year sim")
tF <- 73
dfRuns4 <- data.frame(Size=twoDF$Size, Time=sapply(twoDF$Size, FUN="line2D"), Age="73 year sim")
dfRuns8 <- data.frame(Size=threeDF$Size, Time=sapply(threeDF$Size, FUN="line3D"), Age="73 year sim")
dfRuns2D <- rbind(dfRuns,dfRuns2,dfRuns3,dfRuns4)
dfRuns3D <- rbind(dfRuns5,dfRuns6,dfRuns7,dfRuns8)

p1 <- ggplot(DF, aes(x=Size, y=Time, colour=Model)) + geom_point() + scale_color_manual(values=c("red","blue")) +
  geom_abline(intercept=fit2D[[1]][1], slope=fit2D[[1]][2], color="red") + geom_abline(intercept=fit3D[[1]][1], slope=fit3D[[1]][2], color="blue") +
  ggtitle("Benchmarking") + xlab("Basal Cells") + ylab("Time per Tick (s)") + ylim(0.0,2.0) +
  geom_vline(xintercept = 10000)
p1

p2 <- ggplot(dfRuns2D, aes(x=Size, y=Time, colour=Age)) + geom_line() + ggtitle("2D estimated hours per simulation") + xlab("xDim") + ylab("Hours") +
  geom_vline(xintercept = 10000)
p2

p3 <- ggplot(dfRuns3D, aes(x=Size, y=Time, colour=Age)) + geom_line() + ggtitle("3D estimated hours per simulation") + xlab("xDim") + ylab("Hours") +
  geom_vline(xintercept = 10000)
p3

lay <- rbind(c(1,1),
             c(2,3))

0grid.arrange(p1,p2,p3, layout_matrix = lay)


