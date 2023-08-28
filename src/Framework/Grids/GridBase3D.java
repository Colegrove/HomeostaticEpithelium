package Framework.Grids;

import Framework.Tools.Utils;

/**
 * Created by bravorr on 5/17/17.
 */
abstract class GridBase3D extends GridBase{
    public final int xDim;
    public final int yDim;
    public final int zDim;
    public final int length;
    public boolean wrapX;
    public boolean wrapY;
    public boolean wrapZ;
    GridBase3D(int x,int y,int z,boolean wrapX,boolean wrapY,boolean wrapZ){
        xDim=x;
        yDim=y;
        zDim=z;
        length=x*y*z;
    }
    public int SQwrapI(int x, int y, int z){
        //wraps Coords to proper index
        if(In(x,y,z)) { return I(x,y,z);}
        return I(Utils.ModWrap(x,xDim),Utils.ModWrap(y,yDim),Utils.ModWrap(z,zDim));
    }

    /**
     * gets the index of the voxel at the specified coordinates
     */
    public int I(int x, int y, int z){
        //gets grid index from location
        return x*yDim*zDim+y*zDim+z;
    }

    /**
     * gets the xDim component of the voxel at the specified index
     */
    public int ItoX(int i){
        return i/(yDim*zDim);
    }

    /**
     * gets the yDim component of the voxel at the specified index
     */
    public int ItoY(int i){return (i/zDim)%yDim;}

    /**
     * gets the z component of the voxel at the specified index
     */
    public int ItoZ(int i){
        return i%zDim;
    }

    /**
     * gets the index of the square that contains the specified coordinates
     */
    public int I(double x, double y, double z){
        //gets grid index from location
        return (int)Math.floor(x)*yDim*zDim+(int)Math.floor(y)*zDim+(int)Math.floor(z);
    }

    /**
     * returns whether the specified coordinates are inside the grid bounds
     */
    public boolean In(int x, int y, int z){
        if(x>=0&&x<xDim&&y>=0&&y<yDim&&z>=0&&z<zDim){
            return true;
        }
        return false;
    }

    /**
     * returns whether the specified coordinates are inside the grid bounds
     */
    public boolean In(double x, double y, double z){
        int xInt=(int)Math.floor(x);
        int yInt=(int)Math.floor(y);
        int zInt=(int)Math.floor(z);
        return In(xInt,yInt,zInt);
    }
    /**
     * writes into ret the coordinates of SQs with the centerX and centerY coordinates added
     * returns the number of squares that were written into ret. any coordinates that fall out of bounds are not written
     * @param SQs list of coordinates of the form [xDim,yDim,xDim,yDim,...]
     * @param ret list into which the displaced coordinates will be written
     * @param centerX xDim displacement of coordinates
     * @param centerY yDim displacement of coordinates
     * @param centerZ z displacement of coordinates
     * @param wrapX whether to wrap the coordinates that fall out of bounds in the X dimension
     * @param wrapY whether to wrap the coordinates that fall out of bounds in the Y dimension
     * @param wrapZ whether to wrap the coordinates that fall out of bounds in the Z dimension
     * @return the number of coordinates written into the ret array
     */
    public int SQstoLocalIs(int[] SQs, int[] ret, int centerX, int centerY, int centerZ, boolean wrapX,boolean wrapY,boolean wrapZ) {
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt = 0;
        for (int i = 0; i < SQs.length / 3; i++) {
            int x = SQs[i * 3] + centerX;
            int y = SQs[i * 3 + 1] + centerY;
            int z = SQs[i * 3 + 2] + centerZ;
            if (!Utils.InDim(xDim, x)) {
                if (wrapX) {
                    x = Utils.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Utils.InDim(yDim, y)) {
                if (wrapY) {
                    y = Utils.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            if (!Utils.InDim(zDim, z)) {
                if (wrapZ) {
                    z = Utils.ModWrap(z, zDim);
                } else {
                    continue;
                }
            }
            ret[ptCt]= I(x,y,z);
            ptCt++;
        }
        return ptCt;
    }
    public int SQstoLocalIs(int[] SQs, int[] ret, int centerX, int centerY, int centerZ) {
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt = 0;
        for (int i = 0; i < SQs.length / 3; i++) {
            int x = SQs[i * 3] + centerX;
            int y = SQs[i * 3 + 1] + centerY;
            int z = SQs[i * 3 + 2] + centerZ;
            if (!Utils.InDim(xDim, x)) {
                if (wrapX) {
                    x = Utils.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Utils.InDim(yDim, y)) {
                if (wrapY) {
                    y = Utils.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            if (!Utils.InDim(zDim, z)) {
                if (wrapZ) {
                    z = Utils.ModWrap(z, yDim);
                } else {
                    continue;
                }
            }
            ret[ptCt]= I(x,y,z);
            ptCt++;
        }
        return ptCt;
    }
    public double DistSq(double x1, double y1, double z1, double x2, double y2, double z2, boolean wrapX, boolean wrapY, boolean wrapZ){
        return Utils.DistSq3D(x1,y1,z1,x2,y2,z2,xDim,yDim,zDim,wrapX,wrapY,wrapZ);
    }
    public double DistSq(double x1, double y1, double z1, double x2, double y2, double z2){
        return Utils.DistSq3D(x1,y1,z1,x2,y2,z2,xDim,yDim,zDim,wrapX,wrapY,wrapZ);
    }
}
