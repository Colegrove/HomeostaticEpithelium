package Framework.Grids;
import Framework.Tools.Utils;
//import AgentFramework.Utils;


import java.util.Arrays;

/**
 * GridDiff3 class facilitates 2D diffusion with two arrays of doubles called fields
 * the intended usage is that during a diffusion step, the current values will be read, and the next values will be written to
 * after updates, SwapNextCurr is called to set the next field as the current field.
 */
public class GridDiff3 extends GridBase3D {
    public double[] currField;
    public double[] nextField;
    public double[] scratch;
    //public double[] middleField;

    public GridDiff3(int xDim, int yDim, int zDim,boolean wrapX,boolean wrapY,boolean wrapZ){
        super(xDim,yDim,zDim,wrapX,wrapY,wrapZ);
        int numElements = this.xDim * this.yDim * this.zDim;

        currField = new double[numElements];
        nextField = new double[numElements];
        scratch=null;
        //middleField = new double[numElements];
    }

    public GridDiff3(int xDim, int yDim, int zDim){
        super(xDim,yDim,zDim,false,false,false);

        int numElements = this.xDim * this.yDim * this.zDim;
        currField = new double[numElements];
        nextField = new double[numElements];
        scratch=null;
        //middleField = new double[numElements];
    }

    /**
     * gets the current field value at the specified index
     */
    public double GetCurr(int i){return currField[i];}

    /**
     * gets the current field value at the specified coordinates
     */
    public double GetCurr(int x,int y,int z) {
        return currField[x*yDim*zDim+y*zDim+z];
    }

    /**
     * sets the current field value at the specified index
     */
    public void SetCurr(int i,double val){
        currField[i]=val;}
    /**
     * sets the current field value at the specified coordinates
     */
    public void SetCurr(int x,int y,int z,double val){
        currField[x*yDim*zDim+y*zDim+z]=val;
    }

    //FIXME- copied and pasted by Jake Hogan7/13 from GridDirr2.java
    public double[] GetCurrField(){ return this.currField; }
    public double[] GetNextField(){ return this.nextField; }

    /**
     * adds to the current field value at the specified coordinates
     */

    public void AddCurr(int x,int y,int z,double val){
        currField[x*yDim*zDim+y*zDim+z]+=val;
    }

    /**
     * adds to the current field value at the specified index
     */
    public void AddCurr(int i,double val){
        currField[i]+=val;
    }

    /**
     * sets the current field value at the specified index
     */
    //public void SetMiddle(int i,double val){
    //    middleField[i]=val;}

    /**
     * gets the next field value at the specified index
     */
    //public double GetMiddle(int i){return middleField[i];}

    /**
     * gets to the current field value at the specified coordinates
     */
    public double GetNext(int x,int y,int z){
        return nextField[x*yDim*zDim+y*zDim+z];
    }

    /**
     * gets the next field value at the specified index
     */
    public double GetNext(int i){return nextField[i];}

    /**
     * gets to the current field value at the specified coordinates
     */
    public void SetNext(int x,int y,int z,double val){
        nextField[x*yDim*zDim+y*zDim+z]=val;
    }

    /**
     * sets the current field value at the specified index
     */
    public void SetNext(int i,double val){
        nextField[i]=val;}



    /**
     * adds to the next field value at the specified index
     */
    public void AddNext(int x,int y,int z,double val){
        nextField[x*yDim*zDim+y*zDim+z]+=val;
    }
    public void AddNext(int i, double val) {
        nextField[i] += val;
    }

    /**
     * copies the current field into the next field
     */
    public void NextCopyCurr(){
        System.arraycopy(currField,0, nextField,0, currField.length);
    }

    /**
     * swaps the next and current field
     */
    public void SwapNextCurr(){
        double[] temp= currField;
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
     * Runs diffusion on the current field, putting the results into the next field
     */
    public void Diffuse(double diffRate){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,false,0.0,wrapX,wrapY,wrapZ);
    }
    public void Diffuse(double diffRate,boolean wrapX,boolean wrapY,boolean wrapZ){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,false,0.0,wrapX,wrapY,wrapZ);
    }
    public void Diffuse(double diffRate,double boundaryValue){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,true,boundaryValue,wrapX,wrapY,wrapZ);
    }
    public void Diffuse(double diffRate,double boundaryValue,boolean wrapX,boolean wrapY,boolean wrapZ){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,true,boundaryValue,wrapX,wrapY,wrapZ);
    }
    /**
     * Runs diffusion on the current field, putting the results into the next field, then swaps them
     * @param diffRate rate of diffusion
     */
    public void DiffSwap(double diffRate){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,false,0.0,wrapX,wrapY,wrapZ);
        SwapNextCurr();
    }
    public void DiffSwap(double diffRate,boolean wrapX,boolean wrapY,boolean wrapZ){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,false,0.0,wrapX,wrapY,wrapZ);
        SwapNextCurr();
    }
    public void DiffSwap(double diffRate,double boundaryValue){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,true,boundaryValue,wrapX,wrapY,wrapZ);
        SwapNextCurr();
    }
    public void DiffSwap(double diffRate,double boundaryValue,boolean wrapX,boolean wrapY,boolean wrapZ){
        Utils.Diffusion3(currField, nextField,xDim,yDim,zDim,diffRate,true,boundaryValue,wrapX,wrapY,wrapZ);
        SwapNextCurr();
    }
    public void ADITripleDiffSwap(final double diffRate){
        if(scratch==null){
            scratch=new double[Math.max(Math.max(xDim,yDim),zDim)*2];
        }
        Utils.DiffusionADI3(0,currField,nextField,scratch,xDim,yDim,zDim,diffRate/3);
        SwapNextCurr();
        Utils.DiffusionADI3(1,currField,nextField,scratch,xDim,yDim,zDim,diffRate/3);
        SwapNextCurr();
        Utils.DiffusionADI3(2,currField,nextField,scratch,xDim,yDim,zDim,diffRate/3);
        SwapNextCurr();
    }

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
    public void AddAllCurr(double val){
        for (int i = 0; i < length; i++) {
            currField[i]+=val;
        }
    }

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

    /**
     * Copies the values currently contained in currField into nextField
     */
    public void CurrIntoNext() {
        System.arraycopy(currField,0,nextField,0,length);
    }
}
