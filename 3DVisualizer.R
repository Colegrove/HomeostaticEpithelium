library(rgl)

xDim = 5
yDim = 20
zDim = 5

SQtoX <- function(i){ return( floor(i/(yDim*zDim))+1 )  }
SQtoY <- function(i){ return( floor((i/zDim)%%yDim)+1 )  }
SQtoZ <- function(i){ return( (i%%zDim)+1 ) }

initialize <- function(){
  #df <- data.frame(matrix(NA, ncol = 3,nrow = 0))
  l <- vector("list", xDim+1*yDim+1*zDim+1)
  i = 1
  for(x in 1:xDim){
    for(y in 1:yDim){
      for(z in 1:zDim){
        cellObj <- rgl.spheres(x=x, y=y, z=z, r=0.5, color="white", alpha=0.1)
        l[i] <- cellObj
        i = i+1
      }
    }
  }
  return(l)
}

redraw <- function(mySpheres, i, LineData){
	if(LineData$V4 != 0.0){
		try(rgl.pop(type="shape",id=mySpheres[i]), silent = T)
		if(LineData$V1==1.0 && LineData$V2==1.0 && LineData$V3==1.0){
			mySpheres[i] <- rgl.spheres(x=SQtoX(i), y=SQtoY(i), z=SQtoZ(i), r=0.5, color=rgb(LineData$V1,LineData$V2,LineData$V3,LineData$V4),alpha=0.1)
		} else{
			mySpheres[i] <- rgl.spheres(x=SQtoX(i), y=SQtoY(i), z=SQtoZ(i), r=0.5, color=rgb(LineData$V1,LineData$V2,LineData$V3,LineData$V4))
		}
	} else{
		try(rgl.pop(type="shape",id=mySpheres[i]), silent = T)
	}
	return(mySpheres)
}

#rgl.open()  depreciated changed to open3d() 25AUG23HLC
open3d()
view3d(35,30)
#rgl.bg(color="black") depreciated changed to bg3d
bg3d(color='black')

mySpheres <- initialize()
i = 1
iteration = 0
f <- file("stdin")
open(f)

prevState <- vector("list", length(mySpheres))
while(length(line <- readLines(f,n=1)) > 0) {
	write(line, stderr())
	if(grepl(line, "Done")==F){
        LineData <- as.data.frame(read.table(text=line, sep = "\t",header = F,colClasses = "numeric"))
        if(isTRUE(all.equal(LineData,prevState[[i]]))==FALSE){
        	mySpheres <- redraw(mySpheres,i,LineData)
        }
        prevState[[i]] <- LineData
        i=i+1
	} else {
		rm(i)
		i=1
		iteration = iteration + 1
	}
}
