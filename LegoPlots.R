library(rgl)

## Draws a single "column" or "stack".
## X and Y coordinates determine the area of the column
## The Z coordinate determines the height of the column
## We include "lit=FALSE" arguments to remove the nasty shiny surfaces caused by lighting
stackplot.3d<-function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1=c(rep(0,4),rep(c(0,0,z,z),4))
  y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2=c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z,4),c(y[1],y[1],y[2],y[2]),
            col=rep(topcol,each=4),alpha=1,lit=FALSE) 
  ## This line adds black edges to the column
  rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)
}



context3d<-function(df, alpha=1,scalexy=10,scalez=1,gap=0.2){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Define dimensions of the plot 
  dimensions=c(length(colnames(df)),length(rownames(df)))
  
  ## Scale column area and the gap between columns 
  y=seq(1,dimensions[1])*scalexy
  x=seq(1,dimensions[2])*scalexy
  gap=gap*scalexy
  
  ## Scale z coordinate
  z <- list()
  c <- 1
  for(i in 1:length(y)){
    for(j in 1:length(x)){
      z[c]=df[j,i]*scalez
      c = c + 1
    }
  }
  z=unlist(z)

  ## Set up colour palette
  rbPal <- colorRampPalette(c('blue','red'))
  colors <- rbPal(100)[as.numeric(cut(z,breaks = 100))]
  
  ## Plot each of the columns
  for(i in 1:dimensions[1]){
    for(j in 1:dimensions[2]){
      it=(i-1)*dimensions[2]+j # Variable to work out which column to plot; counts from 1:96
      stackplot.3d(c(gap+x[j],x[j]+scalexy),
                   c(-gap-y[i],-y[i]-scalexy),
                   z[it],
                   alpha=alpha,
                   topcol=colors[it],
                   sidecol=colors[it])
    }
  }
  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=45,phi=35,fov=0)
  #axes3d("y-+",labels=TRUE)
  axes3d(edges="bbox")
}

div3d <- read.csv("~/Desktop/Division3D_2Yrs.txt", sep="\t", header=F)
div3d <- div3d[,1:20]
colnames(div3d) <- seq(1, 20, 1)

stepSize <- 20
to.average <- split(div3d, (seq(nrow(div3d))-1) %/% stepSize)
ave.divs <- data.frame(matrix(NA, nrow=0, ncol=20))
for(i in 1:length(to.average)){
  ave.div <- colMeans(to.average[[i]])
  ave.divs <- rbind(ave.divs,ave.div)
}

context3d(ave.divs, alpha=0.4, scalexy = 1, scalez=100)
par3d("windowRect"= c(0,0,1200,1200))
rgl.viewpoint(theta=145,phi=35,fov=0,  zoom = .9 )
rgl.snapshot(
  "/Users/schencro/Desktop/Darryl_collab/Framework/Homeostatic_Epidermis/3DDivision",
  fmt = "png",
  top=TRUE
  )


