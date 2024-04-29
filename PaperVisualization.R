# Title     : Visualize 3D Simulations
# Objective : Create Paper worthy Cell figures
# Created by: schencro
# Created on: 12/6/17

## Edits made 29Aug23HLC

library(rgl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)

# RGL XYZ correspond to vertices...such that x1,y1,z1 are split by x=c(),y=(),z=()
# in this way it's a very R fashion. Can do complex shapes with this method. x y or z
# is not an individual vertices. For triangles this means a length(vector)/3 must equal 0.

SQtoX <- function(i){ return( floor(i/(yDim*zDim))+1 )  }
SQtoY <- function(i){ return( floor((i/zDim)%%yDim)+1 )  }
SQtoZ <- function(i){ return( (i%%zDim)+1 ) }

SubsectHelper <- function(i, section=1){
    # Section corresponds to right chunk out (1), left chunk out (2), or half out (3).
    if(section==1){
        if( (SQtoX(i)<xDim/2 || SQtoZ(i)<xDim/2)){
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else if(section==2){
        if( (SQtoX(i)>xDim/2 || SQtoZ(i)<xDim/2)){
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else if(section==3){
        if( SQtoZ(i)<zDim/2 ){
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
}
SubsectHelperXYZ <- function(x,z,y, section=1){
  # Section corresponds to right chunk out (1), left chunk out (2), or half out (3).
  if(section==1){
    if( (x<xDim/2 || z<xDim/2)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if(section==2){
    if( (x>xDim/2 || z<xDim/2)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if(section==3){
    if( z<zDim/2 ){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

Subsect <- function(df, section=1){
    dfList <- sapply(df$i, FUN = SubsectHelper, section=section)
    return(df[dfList,])
}
SubsectXYZ <- function(df, section=1){
  dfList <- sapply(df$x,df$y,df$z, FUN = SubsectHelperXYZ, section=section)
  return(df[dfList,])
}

GetData <- function(fileName){
    df <- read.csv(fileName, header=F, sep="\t")
    #colnames(df)<-c("i","h","s","v","alpha")
    colnames(df)<-c("x","z","y","h","s","v","alpha") ## output data is "x,z,y" not "i" 29aug23HLC
    whiteCells <- subset(df,df$h==0.0)
    otherCells <- subset(df,df$h!=0.0)
    whiteCells$alpha = rep(0.1,length(whiteCells$alpha))
    otherCells$alpha = otherCells$alpha-0.2
    df <- rbind(whiteCells,otherCells)
    return(df)
}

UpdateColorscheme <- function(df){ ## adjusts reference cells to purple (#AA4499) friendly 29Aug23HLC
  NewH <- 250/360
  NewS <- 0.75
  NewV <- 0.53
  df$h[df$h == 1 & df$s == 1 & df$v == 1]
  df$replacement_h[df$h == 1 & df$s == 1 & df$v == 1] <- NewH
  df$replacement_s[df$h == 1 & df$s == 1 & df$v == 1] <- NewS
  df$replacement_v[df$h == 1 & df$s == 1 & df$v == 1] <- NewV
  df$replacement_h[df$h == 0 & df$s == 0 & df$v == 0] <- 0
  df$replacement_s[df$h == 0 & df$s == 0 & df$v == 0] <- 0
  df$replacement_v[df$h == 0 & df$s == 0 & df$v == 0] <- 0
  return(df)
}

initializeRGL <- function(theta=0,phi=20,myColor="white"){
    par3d(windowRect = c(0, 0, 1600, 800))
    view3d(theta,phi)
    bg3d(color=myColor)
}

placeCells <- function(theData=data.frame(), adjustx=0, adjustz=0, adjusty = 0, myAlpha=1.0, radius=0.8){
    #x=SQtoX(theData$i) + adjustx
    #z=SQtoZ(theData$i) + adjustz
    #y=SQtoY(theData$i) + adjusty
    x=theData$x + adjustx
    z=theData$z + adjustz
    y=theData$y + adjusty
    spheres3d(x=x,y=y,z=z,radius=radius,color=hsv(theData$h,theData$s,theData$v),alpha=myAlpha)
}

BigBoxPlot <- function(){
    l = (xDim*3+xDim/4*2+5)
    x = c( 0 , 0 , 0 , 0 , 0 , 0, 0 , (xDim*3+xDim/4*2+5) , 0 , (xDim*3+xDim/4*2+5) , 0 ,
    (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5) , 0 , 0 , (xDim*3+xDim/4*2+5)
    , (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5) , (xDim*3+xDim/4*2+5))
    y = c( -50 , yDim+5 , -50 , -50 , yDim+5 , yDim+5 , -50 , -50 , yDim+5 , yDim+5 , -50 , -50 , -50 , -50 , yDim+5 , yDim+5 , -50 , yDim+5 , -50 , yDim+5 , -50 , yDim+5)
    z = c( 0 , 0 , 0 , zDim+5 , 0 , zDim+5 , 0 , 0 , 0 , 0 , zDim+5, zDim+5 , 0 , zDim+5 , 0 , zDim+5 , zDim+5, zDim+5 , zDim+5 , zDim+5 , 0 , 0)
    rgl.lines(x=x,y=y,z=z, color="black")
    #box
    x = c(0,l,l, 0,0,l, 0,0,l, 0,l,l, 0,0,0, 0,0,0, l,l,l, l,l,l)
    y = c(-50,-50,yDim+5, -50,yDim+5,yDim+5, -50,-50,-50, -50,-50,-50, -50,-50,yDim+5, -50,yDim+5,yDim+5, -50,-50,yDim+5, -50,yDim+5,yDim+5)
    z = c(0,0,0, 0,0,0, 0,zDim+5,zDim+5, 0,0,zDim+5, 0,zDim+5,zDim+5, 0,0,zDim+5, 0,zDim+5,zDim+5, 0,0,zDim+5)
    triangles3d(x=x,y=y,z=z,color="black")

    x = c(l/3,l/3,l/3, l/3,l/3,l/3, l/3*2,l/3*2,l/3*2, l/3*2,l/3*2,l/3*2)
    y = c(-50,yDim+5,yDim+5, -50,-50,yDim+5, -50,yDim+5,yDim+5, -50,-50,yDim+5)
    z = c(0,0,zDim/2, 0,zDim/2,zDim/2, 0,0,zDim/2, 0,zDim/2,zDim/2)
    triangles3d(x=x,y=y,z=z,color="red", alpha=1)
}

BasalLayer <- function(df){
  df <- df %>% filter(y == 0) ## basal layer only
  return(df)
}

crossSection <- function(df,x_coord){
  df <- df %>% filter(x == x_coord)
  return(df)
}

RemoveRefDead <- function(df){
  # Filter rows where HSV values are not (0, 0, 0) or (1, 1, 1)
  df_filtered <- df[!(df$h == 0 & df$s == 0 & df$v == 0) & !(df$h == 1 & df$s == 1 & df$v == 1),] 
  return(df_filtered)
}

updateMergedColors <- function(df){
  merged_df = UpdateColorscheme(df)
  print(merged_df)
  merged_df <- merged_df[,!names(merged_df) %in% c("h", "s", "v")]
  colnames(merged_df)[colnames(merged_df) == "replacement_h"] ="h"
  colnames(merged_df)[colnames(merged_df) == "replacement_s"] ="s"
  colnames(merged_df)[colnames(merged_df) == "replacement_v"] ="v"
  return(merged_df)
}
#====Main====#
xDim=20
yDim=20
zDim=xDim

#t25 <- GetData("~/Desktop/100xDim.25yrs.txt")
t365 <- GetData("/Volumes/gs-vol1/home/huntc10/dat/HomeostaticEpithelium/2023-09-21-BlockingProb_long/results/VisFile_block_1.0_1.1460.txt")

t365_basal = BasalLayer(t365) ## basal layer only
basal_mutant_only = RemoveRefDead(t365_basal) ## mutants in basal layer
mutant_only_all = RemoveRefDead(t365) ## mutants in entire tissue
unique_hsv <- unique(basal_mutant_only[, c("h", "s", "v")]) ## generates unique clones in basal only
unique_hsv_all <- unique(mutant_only_all[, c("h", "s", "v")]) ## generates unique clones in entire tissue
num_unique <- nrow(unique_hsv) ## number of mutant clones basal only
num_unique_all <- nrow(unique_hsv_all) ## number of mutant clones entire tissue
# Create a mapping table for corrected clones to change colors
replacement_values <- data.frame(
  #hsv(350,50,80) # rose
  #hsv(250,75,53) # indigo
  #hsv(0.861,0.6,0.67)
  #hsv(50,46,87) # sand
  #hsv(140,86,47) # green
  #hsv(200,43,93) # cyan
  #hsv(330,75,53) # wine
  #hsv(170,60,67) # teal
  #hsv(60,67,60) # olive
  #hsv(0,0,87) # grey
  replacement_h = rep(c(350/360, 0.861, 50/360, 140/360, 200/360, 330/360, 170/360, 60/360, 0), 
                      length.out = num_unique),
  replacement_s = rep(c(50/100, .6, 46/100, 86/100, 43/100, 75/100, 60/100, 67/100, 0/100), 
                      length.out = num_unique),
  replacement_v = rep(c(80/100, .67, 87/100, 47/100, 93/100, 53/100, 67/100, 60/100, 87/100), 
                      length.out = num_unique)
)
# Create a mapping table for corrected clones to change colors
replacement_values_all <- data.frame(
  #hsv(350,50,80) # rose
  #hsv(250,75,53) # indigo
  #hsv(0.861,0.6,0.67)
  #hsv(50,46,87) # sand
  #hsv(140,86,47) # green
  #hsv(200,43,93) # cyan
  #hsv(330,75,53) # wine
  #hsv(170,60,67) # teal
  #hsv(60,67,60) # olive
  #hsv(0,0,87) # grey
  replacement_h = rep(c(350/360, 0.861, 50/360, 140/360, 200/360, 330/360, 170/360, 60/360, 0), 
                      length.out = num_unique_all),
  replacement_s = rep(c(50/100, .6, 46/100, 86/100, 43/100, 75/100, 60/100, 67/100, 0/100), 
                      length.out = num_unique_all),
  replacement_v = rep(c(80/100, .67, 87/100, 47/100, 93/100, 53/100, 67/100, 60/100, 87/100), 
                      length.out = num_unique_all)
)
color_mapping <- cbind(unique_hsv, replacement_values)
color_mapping_all <- cbind(unique_hsv_all, replacement_values_all)
print(color_mapping)

# Update the original dataframe with replacement HSV values
merged_df <- merge(t365_basal, color_mapping, by = c("h", "s", "v"), all.x = TRUE)
merged_df_all <- merge(t365, color_mapping_all, by = c("h", "s", "v"), all.x = TRUE)
# Print the updated dataframe
merged_df = updateMergedColors(merged_df)
merged_df_all = updateMergedColors(merged_df_all)


initializeRGL()
clear3d()
placeCells(merged_df, myAlpha=merged_df$alpha) ## basal layer only
view3d(zoom=0.7, phi=90) # top down

clear3d()
merged_df_all_Xsection = crossSection(merged_df_all, 16)
#placeCells(merged_df_all, myAlpha=merged_df_all$alpha) ## entire tissue
placeCells(merged_df_all_Xsection, myAlpha=merged_df_all_Xsection$alpha) ## cross section only
view3d(zoom=0.7, phi=0, theta=90) # side

clear3d()



t25_y0 <- UpdateColorscheme(t25_y0)
#t25 <- GetData("VisFile.txt.2.txt")
tibble(t25_y0) %>%  ggplot() + geom_point(aes(x = x, y = z, col = rgb(h, s, v, alpha= alpha))) + facet_wrap(~y)
tibble(t25) %>%  ggplot() + geom_point(aes(x = x, y = z, col = rgb(h, s, v, alpha= alpha))) + facet_wrap(~y)
#t25 %>% filter(between(x, 4, 6), between(y, 4, 6), z== 0)
#t25 <- UpdateColorscheme(t25)
#t25cut <- Subsect(t25,section=1)
t25cut <- SubsectXYZ(t25, section=1)
#t50 <- GetData("~/Desktop/100xDim.50yrs.txt")
#t50cut <- Subsect(t50,section=3)
#t75 <- GetData("~/Desktop/100xDim.75yrs.txt")
#t75cut <- Subsect(t75,section=2)

#drawCells
initializeRGL()

placeCells(t25, myAlpha=t25$alpha)
placeCells(t25_y0, myAlpha=t25_y0$alpha) ## basal layer only
clear3d()
placeCells(merged_df, myAlpha=merged_df$alpha) ## basal layer only
placeCells(t25cut, myAlpha=t25cut$alpha, adjusty=-50)

placeCells(t50, myAlpha=t50$alpha, adjustx=xDim+xDim/4)
placeCells(t50cut, myAlpha=t50cut$alpha, adjusty=-50, adjustx=xDim+xDim/4)

placeCells(t75, myAlpha=t75$alpha, adjustx=xDim*2+(xDim/4)*2)
placeCells(t75cut, myAlpha=t75cut$alpha, adjusty=-50, adjustx=xDim*2+(xDim/4)*2)

view3d(zoom=0.7, phi=90, theta=0)

BigBoxPlot()
#Fig 2A
rgl.snapshot( "3D.VisFile_block_0.0_3.365.basal.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_0.6_40.365.basal.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_1.0_1.365.basal.png", fmt = "png", top = TRUE )
#Fig 2B
rgl.snapshot( "3D.VisFile_block_0.5_50.1460.basal.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_0.5_50.2920.basal.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_0.5_50.4380.basal.png", fmt = "png", top = TRUE )
#Fig 3B
rgl.snapshot( "3D.VisFile_block_1.0_1.365.cross_section_16.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_1.0_1.1460.cross_section_16.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.VisFile_block_1.0_1.2920.cross_section_16.png", fmt = "png", top = TRUE )



rgl.snapshot( "3D.progression.png", fmt = "png", top = TRUE )
rgl.snapshot( "3D.progression_basal.png", fmt = "png", top = TRUE )

clear3d()
rgl.close()

#---End Paper Visualization---#




#---Movie Picture Time---#
GetData2 <- function(fileName){
    df <- read.csv(fileName, header=F, sep="\t")
    colnames(df)<-c("i","h","s","v","alpha","time")
    df$time <- round(df$time,2)
    whiteCells <- subset(df,df$h==0.0)
    otherCells <- subset(df,df$h!=0.0)
    whiteCells$alpha = rep(0.1,length(whiteCells$alpha))
    otherCells$alpha = otherCells$alpha-0.2
    df <- rbind(whiteCells,otherCells)
    return(df)
}

initializeRGL <- function(theta=40,phi=20,myColor="black"){
    par3d(windowRect = c(0, 0, 720, 720))
    view3d(theta,phi)
    rgl.bg(color=myColor)
}

xDim=50
yDim=20
zDim=xDim

df <- GetData("~/Desktop/EpidermisMovie/EpiVisWounding.parsed.txt")
out <- split( df , f = df$time )
#prep window
initializeRGL()

for(i in 1:length(out)){
    placeCells(theData=out[[i]], myAlpha=out[[i]]$alpha)
    #Sys.sleep(0.1)
    f=paste("~/Desktop/EpidermisMovie/imgsWounding/3D.video.wounded.",unique(out[[i]]$time),".png",sep="")
    rgl.snapshot( filename = f, fmt = "png", top = TRUE )
    rgl.clear()
}
rgl.close()

initializeRGL()
startFrame <- data.frame(i=seq(1,(xDim*zDim*yDim)))
startFrame$h <- 0.0
startFrame$s <- 0.0
startFrame$v <- 1.0
startFrame$alpha <- 0.1
placeCells(theData=startFrame, myAlpha=startFrame$alpha)
rgl.snapshot( filename = "~/Desktop/EpidermisMovie/StartBlock.png", fmt = "png", top = TRUE )
rgl.clear()

placeCells(theData=out[[228]], myAlpha=out[[228]]$alpha)
rgl.snapshot( filename = "~/Desktop/EpidermisMovie/FinalFrame.png", fmt = "png", top = TRUE )
rgl.clear()
rgl.close()

#---Surface Plot Movie Picture Time---#
PlotSurface <- function(theData){
    myData <- acast(theData, x~z, value.var="EGF")
    y <- exp(myData)*4
    ylim <- range(y)
    ylen <- ylim[2] - ylim[1] + 1

    colorlut <- heat.colors(ylen) # height color lookup table

    col <- colorlut[ y - ylim[1] + 1 ] # assign colors to heights for each point

    rgl.surface(seq(1,xDim), seq(1,zDim), myData,color=col, back = "lines")
}

df <- read.csv("~/Desktop/EpidermisMovie/EpiVisEGF.parsed.txt", header=F, sep="\t")
colnames(df)<-c("x","z","EGF","time")
df$time <- round(df$time,2)
out <- split( df , f = df$time )

initializeRGL()

for(i in 1:length(out)){
    PlotSurface(out[[i]])
    #Sys.sleep(0.1)
    f=paste("~/Desktop/EpidermisMovie/imgsSurfaceWounding/3D.wounding.egf.",unique(out[[i]]$time),".png",sep="")
    rgl.snapshot( filename = f, fmt = "png", top = TRUE )
    rgl.clear()
}