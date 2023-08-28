package Epidermis_Model;

/**
 * Created by schencro on 8/14/17.
 */

public class LossReplace {
    /*
    The loss replacement rate of the progenitor pool is defined as 0.51 per week within a basal cell density of 10,000 in a mm^2 area.
    In this mm^2 basal layer 1/3 (33.3%) are defined as progenitor cells.
    0.51pw * .333 = .16983pw (adjusted for a system where the basal layer is 100% possible progenitors).
    .16983pw / 7d = 0.02426143pd. This is my target average.
     */
    private int DivisionBasal = 0;
    private int DivisionTissue = 0;
    private int DeathBasal = 0;
    private int DeathTissue = 0;
    private double[] LossRateBasal;
    private double[] LossRateTissue;
    private double[] BirthRateBasal;
    private double[] BirthRateTissue;
    private double[] OutputRateArray;
    private double[] tmpArray;
    EpidermisGrid theGrid;

    public LossReplace(EpidermisGrid GridObj, int ModelTime, int Step){ // Initialize the LossReplacement Fetcher
        this.theGrid = GridObj;
        this.LossRateBasal = new double[ModelTime+1];
        this.LossRateTissue = new double[ModelTime+1];
        this.BirthRateBasal = new double[ModelTime+1];
        this.BirthRateTissue = new double[ModelTime+1];
        this.OutputRateArray = new double[ModelTime/Step + 1];
    }

    public double RecordBasalRate(String Option){
        double Rate=0.0;
        if(Option=="Birth"){
            Rate=(DivisionBasal*1.0)/ (theGrid.xDim*theGrid.zDim);
            BirthRateBasal[theGrid.GetTick()] = Rate; // Gives a per day rate
        } else if(Option=="Death"){
            Rate=(DeathBasal*1.0)/ (theGrid.xDim*theGrid.zDim);
            LossRateBasal[theGrid.GetTick()] = Rate; // Gives a per day rate
        }
        ResetBasalCounters();
        return Rate;
    }

    public double RecordTissueRate(String Option){
        double Rate=0.0;
        if(Option=="Birth"){
            Rate=(DivisionTissue*1.0)/ (theGrid.xDim*theGrid.zDim); // Gives a per day rate
            BirthRateTissue[theGrid.GetTick()] = Rate;
        } else if(Option=="Death"){
            Rate=(DeathTissue*1.0)/ (theGrid.xDim*theGrid.zDim); // Gives a per day rate
            LossRateTissue[theGrid.GetTick()] = Rate;
        }
        ResetTissueCounters();
        return Rate;
    }

    public double GetBasalRate(String Option, int Step){
        double Rate=0.0;
        for (int i = theGrid.GetTick()-Step; i < theGrid.GetTick(); i++) {
            if(Option=="Birth"){
                Rate += BirthRateBasal[i];
            } else if(Option=="Death") {
                Rate += LossRateBasal[i];
            }
        }
        OutputRateArray[theGrid.GetTick()/7] = Rate;
        return Rate/Step; // Returns average daily rate over a week
    }

    public double GetTissueRate(String Option){
        return 0.0;
    }

    public void ResetBasalCounters(){
        DivisionBasal=0;
        DeathBasal=0;
    }

    public void ResetTissueCounters(){
        DivisionTissue=0;
        DeathTissue=0;
    }

    public void RecordDivideBasal(){
        DivisionTissue++;
    }

    public void RecordDivideTissue(){
        DivisionTissue++;
    }

    public void RecordLossBasal() {
        DeathBasal++;
    }

    public void RecordLossTissue() {
        DeathTissue++;
    }

    public int GetBasalDivisions(){ return DivisionBasal; }

    public int GetTissueDivisions(){ return DivisionTissue; }

    public int GetBasalDeath(){ return DeathBasal; }

    public int GetTissueDeath(){ return DeathTissue; }

    public double[] GetBirthRateBasal() { return BirthRateBasal; }

    public double[] GetBirthRateTissue() { return BirthRateTissue; }

    public double[] GetDeathRateBasal() { return LossRateBasal; }

    public double[] GetDeathRateTissue() { return LossRateTissue; }

    public double[] GetOutputArray() { return OutputRateArray; }

}
