package Framework.Grids;
import Framework.Tools.FileIO;
import Framework.Tools.Utils;

import java.util.Arrays;

/**
 * GridDiff2 class facilitates 2D diffusion with two arrays of doubles called fields
 * the intended usage is that during a diffusion step, the current values will be read, and the next values will be written to
 * after updates, SwapNextCurr is called to set the next field as the current field.
 */
public class GridDiff2 extends GridBase2D {
    double[] currField;
    double[] nextField;
    double[] scratch;
    public GridDiff2(int xDim, int yDim){
        super(xDim,yDim,false,false);
        currField =new double[this.xDim * this.yDim];
        nextField =new double[this.xDim * this.yDim];
    }
    public GridDiff2(int xDim, int yDim,boolean wrapX,boolean wrapY){
        super(xDim,yDim,wrapX,wrapY);
        currField =new double[this.xDim * this.yDim];
        nextField =new double[this.xDim * this.yDim];
    }

    /**
     * gets the current field value at the specified index
     */
    public double GetCurr(int i){return currField[i];}

    public void ADIDiffSwap(boolean xAxis,double diffRate){
        diffRate=diffRate/2;
        if(scratch==null){
            scratch=new double[Math.max(xDim,yDim)*2];
        }
        Utils.DiffusionADI2(xAxis, currField, nextField,scratch,xDim,yDim,diffRate);
        SwapNextCurr();
    }

    public void ADIDoubleDiffSwap(double diffRate){
        if(scratch==null){
            scratch=new double[Math.max(xDim,yDim)*2];
        }
        Utils.DiffusionADI2(true, currField, nextField,scratch,xDim,yDim,diffRate);
        SwapNextCurr();
        Utils.DiffusionADI2(false, currField, nextField,scratch,xDim,yDim,diffRate);
        SwapNextCurr();
    }
    public double[] GetCurrField(){
        return this.currField;
    }
    public double[] GetNextField(){
        return this.nextField;
    }

    /**
     * gets the current field value at the specified coordinates
     */
    public double GetCurr(int x, int y) { return currField[x*yDim+y]; }

    public double GradientX(int centerX,int centerY){
        double x1=centerX+1<xDim?GetCurr(centerX+1,centerY):GetCurr(centerX,centerY);
        double x2=centerX-1>=0?GetCurr(centerX-1,centerY):GetCurr(centerX,centerY);
        return x1-x2;
    }
    public double GradientY(int centerX,int centerY){
        double y1=centerY+1<yDim?GetCurr(centerX,centerY+1):GetCurr(centerX,centerY);
        double y2=centerY-1>=0?GetCurr(centerX,centerY-1):GetCurr(centerX,centerY);
        return y1-y2;
    }

    /**
     * sets the current field value at the specified index
     */
    public void SetCurr(int i,double val){
        currField[i]=val;}

    /**
     * sets the current field value at the specified coordinates
     */
    public void SetCurr(int x, int y, double val){ currField[x*yDim+y]=val; }

    /**
     * adds to the current field value at the specified coordinates
     */
    public void AddCurr(int x, int y, double val){ currField[x*yDim+y]+=val; }

    /**
     * adds to the current field value at the specified index
     */
    public void AddCurr(int i,double val){
        currField[i]+=val;}

    /**
     * gets the next field value at the specified coordinates
     */
    public double GetNext(int x,int y){ return nextField[x*yDim+y]; }

    /**
     * gets the next field value at the specified index
     */
    public double GetNext(int i){return nextField[i];}

    /**
     * sets the next field value at the specified coordinates
     */
    public void SetNext(int x,int y,double val){ nextField[x*yDim+y]=val; }

    /**
     * sets the next field value at the specified index
     */
    public void SetNext(int i,double val){
        nextField[i]=val;}

    /**
     * sets the next field value at the specified coordinates
     */
    public void AddNext(int x,int y,double val){ nextField[x*yDim+y]+=val; }

    /**
     * adds to the next field value at the specified index
     */
    public void AddNext(int i,double val){
        nextField[i]+=val;}

    /**
     * copies the current field into the next field
     */
    public void CurrIntoNext(){ System.arraycopy(currField, 0, nextField, 0, currField.length); }

    public void CurrAddToNext(){
        for (int i = 0; i < length; i++) {
            currField[i]+=nextField[i];
        }
    }

    /**
     * Bounds all values in the current field between min and max
     */
    public void BoundAllCurr(double min,double max){
        for(int i=0;i<length;i++){
            currField[i]= Utils.BoundVal(currField[i],min,max);
        }
    }
    /**
     * Bounds all values in the next field between min and max
     */
    public void BoundAllNext(double min,double max){
        for(int i=0;i<length;i++){
            nextField[i]= Utils.BoundVal(nextField[i],min,max);
        }
    }
    /**
     * Swaps the next and current field
     */
    public void SwapNextCurr(){
        double[]temp= currField;
        currField = nextField;
        nextField =temp;

    }

    /**
     * Swaps the next and current field, and increments the tick
     */
    public void SwapInc(){
        SwapNextCurr();
        IncTick();
    }
    /**
     * Runs diffusion on the current field, putting the result into the next field
     * @param diffRate rate of diffusion
     */
    public void Diffuse(double diffRate){
        //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,false,0.0,wrapX,wrapY);
    }
    public void Diffuse(double diffRate,boolean wrapX,boolean wrapY){
        //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,false,0.0,wrapX,wrapY);
    }
    public void Diffuse(double diffRate,double boundaryValue){
        //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,true,boundaryValue,wrapX,wrapY);
    }
    public void Diffuse(double diffRate,double boundaryValue,boolean wrapX,boolean wrapY){
        //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,true,boundaryValue,wrapX,wrapY);
    }

    public void Convection(double xVel,double yVel){
        Utils.Advection2(currField,nextField,xDim,yDim,xVel,yVel,false,0.0);
    }
    public void Convection(double xVel,double yVel,double boundaryValue){
        Utils.Advection2(currField,nextField,xDim,yDim,xVel,yVel,true,boundaryValue);
    }
    public void ConvSwap(double xVel,double yVel){
        Utils.Advection2(currField,nextField,xDim,yDim,xVel,yVel,false,0.0);
        SwapNextCurr();
    }
    public void ConvSwap(double xVel,double yVel,double boundaryValue){
        Utils.Advection2(currField,nextField,xDim,yDim,xVel,yVel,true,boundaryValue);
        SwapNextCurr();
    }
    public void ConvectionInhomogeneous(double[] xVels,double[] yVels,double boundaryValue){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,true,boundaryValue);
    }
    public void ConvectionInhomogeneous(double[] xVels,double[] yVels){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,false,0.0);
    }
    public void ConvInhomogeneousSwap(double[] xVels,double[] yVels,double boundaryValue){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,true,boundaryValue);
        SwapNextCurr();
    }
    public void ConvInhomogeneousSwap(double[] xVels,double[] yVels){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,false,0.0);
        SwapNextCurr();
    }
    public void ConvVolumeConservingSwap(double[] xVels,double[]yVels){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,false,0.0);
        Utils.ConservativeTransportStep2(currField,nextField,nextField,xDim,yDim,xVels,yVels,false,0.0);
        SwapNextCurr();
    }
    public void ConvVolumeConservingSwap(double[] xVels,double[]yVels,double boundaryValue){
        Utils.Advection2(currField,nextField,xDim,yDim,xVels,yVels,true,boundaryValue);
        Utils.ConservativeTransportStep2(currField,nextField,nextField,xDim,yDim,xVels,yVels,true,boundaryValue);
        SwapNextCurr();
    }

    /**
     * Runs diffusion on the current field, putting the result into the next field, then swaps current and next
     * @param diffRate rate of diffusion
     */
    public void DiffSwap(double diffRate){
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,false,0.0,wrapX,wrapY);
        SwapNextCurr();
    }
    public void DiffSwap(double diffRate,boolean wrapX,boolean wrapY){
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,false,0.0,wrapX,wrapY);
        SwapNextCurr();
    }
    /**
     * Runs diffusion on the current field, putting the result into the next field, then swaps current and next
     * @param diffRate rate of diffusion
     * @param boundaryValue value that diffuses in from the boundary
     */
    public void DiffSwap(double diffRate,double boundaryValue){
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,true,boundaryValue,wrapX,wrapY);
        SwapNextCurr();
    }
    public void DiffSwap(double diffRate,double boundaryValue,boolean wrapX,boolean wrapY){
        Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,true,boundaryValue,wrapX,wrapY);
        SwapNextCurr();
    }


    public void DiffuseInhomogeneous(double[] diffRates,boolean wrapX,boolean wrapY){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,false,0.0,wrapX,wrapY);
    }
    public void DiffuseInhomogeneous(double[] diffRates){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,false,0.0,wrapX,wrapY);
    }
    public void DiffuseInhomogeneous(double[] diffRates,double boundaryValue){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,true,boundaryValue,wrapX,wrapY);
    }
    public void DiffInhomogeneousSwap(double[] diffRates,boolean wrapX,boolean wrapY){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,false,0.0,wrapX,wrapY);
        SwapNextCurr();
    }
    public void DiffInhomogeneousSwap(double[] diffRates){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,false,0.0,wrapX,wrapY);
        SwapNextCurr();
    }
    public void DiffInhomogeneousSwap(double[] diffRates,double boundaryValue){
        Utils.Diffusion2inhomogeneous(currField, nextField,diffRates,xDim,yDim,true,boundaryValue,wrapX,wrapY);
        SwapNextCurr();
    }

    /**
     * Runs diffusion on the current field, putting the result into the next field, then swaps current and next, and increments the tick
     * @param diffRate rate of diffusion
     * @param boundaryCond whether a boundary condition value will diffuse in from the field boundaries
     * @param boundaryValue only applies when boundaryCond is true, the boundary condition value
     * @param wrapX whether to wrap the field over the left and right boundaries
     */
  //  public void DiffSwapInc(double diffRate,boolean boundaryCond,double boundaryValue,boolean wrapX,boolean wrapY){
  //      //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
  //      Utils.Diffusion2(currField, nextField,xDim,yDim,diffRate,boundaryCond,boundaryValue,wrapX,wrapY);
  //      SwapNextCurr();
  //      IncTick();
  //  }
  //  public void DiffSwapInc(double diffRate,boolean boundaryCond,double boundaryValue){
  //      //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
  //      Utils.Diffusion2(currField, nextField,xDim,yDim,diffRate,boundaryCond,boundaryValue,wrapX,wrapY);
  //      SwapNextCurr();
  //      IncTick();
  //  }
  //  public void DiffSwapInc1(double diffRate,boolean boundaryCond,double boundaryValue,boolean wrapX,boolean wrapY){
  //      //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
  //      Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,boundaryCond,boundaryValue,wrapX,wrapY);
  //      SwapNextCurr();
  //      IncTick();
  //  }
  //  public void DiffSwapInc1(double diffRate,boolean boundaryCond,double boundaryValue){
  //      //NOTE: EXPLICIT DIFFUSION WILL ONLY BE STABLE IF diffRate <= 1/4
  //      Utils.Diffusion(currField, nextField,xDim,yDim,diffRate,boundaryCond,boundaryValue,wrapX,wrapY);
  //      SwapNextCurr();
  //      IncTick();
  //  }

    /**
     * returns the maximum difference between the current field and the next field
     * @param scaled divides the differences by the current field value (unstable at low concentrations)
     */
    public double MaxDiff(boolean scaled){
        double maxDiff=0;
        if(!scaled){
            for(int i = 0; i< currField.length; i++){
                maxDiff=Math.max(maxDiff,Math.abs(currField[i]- nextField[i]));
            }
        }else{
            for(int i = 0; i< currField.length; i++){
                maxDiff=Math.max(maxDiff,Math.abs((currField[i]- nextField[i])/ currField[i]));
            }
        }
        return maxDiff;
    }

    /**
     * sets all squares in current the field to the specified value
     */
    public void SetAllCurr(double val){
        Arrays.fill(currField,val);
    }

    /**
     * sets all squares in the next field to the specified value
     */
    public void SetAllNext(double val){
        Arrays.fill(nextField,val);
    }

    /**
     * adds specified value to all entries of the curr field
     */
    public void AddAllCurr(double val){
        for (int i = 0; i < length; i++) {
            currField[i]+=val;
        }
    }

    /**
     * adds specified value to all entries of the next field
     */
    public void AddAllNext(double val){
        for (int i = 0; i < length; i++) {
            nextField[i]+=val;
        }
    }

    /**
     * gets the average value of all squares in the current field
     */
    public double AvgCurr(){
        double tot=0;
        for(int i=0;i<length;i++){
            tot+= currField[i];
        }
        return tot/length;
    }
    public double MaxCurr(){
        double max=Double.MIN_VALUE;
        for (int i = 0; i < length; i++) {
            max=Math.max(GetCurr(i),max);
        }
        return max;
    }
    public double MinCurr(){
        double min=Double.MAX_VALUE;
        for (int i = 0; i < length; i++) {
            min=Math.max(GetCurr(i),min);
        }
        return min;
    }

    /**
     * gets the average value of all squares in the next field
     */
    public double AvgNext(){
        double tot=0;
        for(int i=0;i<length;i++){
            tot+= nextField[i];
        }
        return tot/length;
    }

    public void SetOuterLayerCurr(double val){
        for (int x = 0; x < xDim; x++) {
            SetCurr(x,0,val);
            SetCurr(x,yDim-1,val);
        }
        for (int y = 1; y < yDim; y++) {
            SetCurr(0,y,val);
            SetCurr(xDim-1,y,val);
        }
    }
    public void SetOuterLayerNext(double val){
        for (int x = 0; x < xDim; x++) {
            SetNext(x,0,val);
            SetNext(x,yDim-1,val);
        }
        for (int y = 1; y < yDim; y++) {
            SetNext(0,y,val);
            SetNext(xDim-1,y,val);
        }
    }

    /**
     * writes into ret the coordinates of SQs with the centerX and centerY coordinates added
     * returns the number of squares that were written into ret. any coordinates that fall out of bounds are not written
     * @param SQs list of coordinates of the form [xDim,yDim,xDim,yDim,...]
     * @param ret list into which the displaced coordinates will be written
     * @param centerX xDim displacement of coordinates
     * @param centerY yDim displacement of coordinates
     * @param wrapX whether to wrap the coordinates that fall out of bounds in the X direction
     * @param wrapY whether to wrap the coordinates that fall out of bounds in the Y direction
     * @return the number of coordinates written into the ret array
     */
    public int SQsToLocalIs(int[] SQs, int[] ret, int centerX, int centerY, boolean wrapX, boolean wrapY){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<SQs.length/2;i++) {
            int x = SQs[i * 2] + centerX;
            int y = SQs[i * 2 + 1] + centerY;
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
            ret[ptCt]= I(x,y);
            ptCt++;
        }
        return ptCt;
    }
    public int SQsToLocalIs(int[] SQs, int[] ret, int centerX, int centerY){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<SQs.length/2;i++) {
            int x = SQs[i * 2] + centerX;
            int y = SQs[i * 2 + 1] + centerY;
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
            ret[ptCt]= I(x,y);
            ptCt++;
        }
        return ptCt;
    }
    void CurrToMatrixCSV(FileIO out){
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                if(y==yDim-1){
                    out.Write(GetCurr(x,y)+"");
                }
                else{
                    out.Write(GetCurr(x,y)+",");
                }
            }
            out.Write("\n");
        }
    }
}
