package Epidermis_Model;


import Epidermis_Model.Genome.GenomeTracker;
import Framework.Grids.Grid3;
import Framework.Grids.GridDiff3;
import Framework.Gui.GuiGridVis;
import Framework.Tools.FileIO;
import Framework.Tools.Utils;

//import com.sun.javafx.util.Utils; commented out 25AUG23HLC unclear when this is used, may need to get javafx and point

import static Epidermis_Model.EpidermisConst.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.List;

/**
 * Created by schencro on 3/31/17.
 */


// Grid specific parameters
class EpidermisGrid extends Grid3<EpidermisCell> {
    final Random RN=new Random(EpidermisConst.Replicate*EpidermisConst.RecordTime);
    static final int[] moveHood={1,0,0, -1,0,0, 0,0,1, 0,0,-1, 0,-1,0};
    //static final int[] moveHood={1,0,0, -1,0,0, 0,0,1, 0,0,-1,0,1,0};
    static final int[] divHood={1,0,0, -1,0,0, 0,0,1, 0,0,-1, 0,1,0};
    static final int[] inBounds= new int[5];
    static final int MOVE=1;
    static final int DIV=2;
    static final double EGF_DIFFUSION_RATE=0.02739064; //keratinocyte growth factor
    static final double DECAY_RATE=0.0007718750; //chemical decay rate of growth factors
    static final double SOURCE_EGF=1; //constant level at basement
    static final int AIR_HEIGHT=15; //air, keratinocyte death! (threshold level for placement of keratinocytes essentially)
    static final int CHEMICAL_STEPS=100; // number of times diffusion is looped every tick
    public int[] divisions = new int[ModelTime*ySize];
    public int divs = 0;
    static double[][][][] ImageArray = new double[EpidermisConst.ySize][EpidermisConst.xSize][EpidermisConst.zSize][7];
    //static int[][][] divArray = new int[EpidermisConst.xSize][EpidermisConst.zSize][ModelTime];

    boolean running;
    int xDim;
    int yDim;
    int zDim;
    long popSum=0;
    int[] MeanProlif = new int[EpidermisConst.xSize * EpidermisConst.ySize * EpidermisConst.zSize];
    int[] MeanDeath = new int[EpidermisConst.xSize * EpidermisConst.ySize * EpidermisConst.zSize];
    GenomeTracker<EpidermisCellGenome> GenomeStore;
    LossReplace Turnover;
    GridDiff3 EGF;

    static boolean corrected = true; // changes to false during correction timestep and changes back to true afterword
                                    // so no further correction occurs 29Aug23HLC
    static List<int[]> InjectionSites = new ArrayList<>(); // create a list of coordinates for center of injection site 29Aug23HLC
    static List<double[]> CorrectionProbs = new ArrayList<>(); // create a list of coordinates containing correction probabilities 14Nov23HLC
    static List<int[]> CorrectionPoints = new ArrayList<>(); // create a list of coordinates for cells to correct 14Nov23HLC

    public EpidermisGrid(int x, int y, int z) {
        super(x,y,z,EpidermisCell.class);
        running = false;
        xDim = x;
        yDim = y;
        zDim = z;
        EGF = new GridDiff3(x, y, z);
        GenomeStore = new GenomeTracker<>(new EpidermisCellGenome(0f,0f,1f,"", this, -1), true, true);
        Turnover = new LossReplace(this, ModelTime, 7);
        PlaceCells();
    }

    public void PlaceCells() {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < AIR_HEIGHT; y++) {
                for (int z = 0; z < xDim; z++) {
                    if (GetAgent(x, y, z) == null) {
                        EpidermisCell c = NewAgentSQ(x, y, z);
                        c.init(KERATINOCYTE, GenomeStore.NewProgenitor(), tp53CloneTracker); // Initializes cell types; Uniform Start
                    }
                }
            }
        }
    }

    public void confluenceCheck(){
        int basalCorrCount = 0;
        for (int x = 0; x < EpidermisConst.xSize; x++) {
            for (int z = 0; z < EpidermisConst.zSize; z++) {
                EpidermisCell c = GetAgent(x, 0, z);
                if(c !=null ){
                    String genomeString = c.myGenome.FullLineageInfoStr("");
                    if (genomeString.contains("FAcorrection")) {
                        basalCorrCount++;
                    }
                }
            }
        }
        if(stopAtConfluence) {
            if (basalCorrCount >= xSize * xSize * 0.8) {
                confluence = true;
            }
        }
        if(stopAtLoss) {
            if(stopLossTime == 0) {
                if (basalCorrCount == 0 && GetTick() > correctionTime) {
                    loss = true;
                }
            } else if (stopLossTime <= GetTick()/(364*4.5)){
                if (basalCorrCount == 0 && GetTick() > correctionTime) {
                    loss = true;
                }
            }
        }
    }

    // For each microneedle injection site, determine which surrounding cells are corrected based on bivariate normal distribution
    public void diffusionCorrectionPoints(){
        // open correction probabilities file and add coordinates and probabilites to a list
        try{
            File myFile = new File(corrProbFile);
            Scanner myReader = new Scanner(myFile);
            while (myReader.hasNextLine()){
                String data = myReader.nextLine();
                String[] values = data.split(" ");
                int pointx = Integer.parseInt(values[0]);
                int pointz = Integer.parseInt(values[1]);
                double corrProb = Double.parseDouble(values[2]);
                CorrectionProbs.add(new double[]{pointx,0,pointz,corrProb}); // adjusted by 49
            }
        } catch (FileNotFoundException e){
            System.out.print("Error: Correction probabilities file not found.");
            e.printStackTrace();
        }

        // Draw an n number of cells with replacement to be corrected based on their probability of corrections.

        // sort the correction Probabilites into ascending order
        CorrectionProbs.sort(Comparator.comparingDouble(arr -> arr[3]));
        // compute the cumulative probabilities
        for (int i = 1; i < CorrectionProbs.size(); i++) {
            double[] prob = CorrectionProbs.get(i);
            double[] prob_minus = CorrectionProbs.get(i-1);
            prob[3] += prob_minus[3];
        }
        // How many cells are to be corrected at each microneedle
        int n = EpidermisConst.dose;

        // Adjust generic coordinates based on each injection site
        for (int[] site:InjectionSites) {

            // Find the samples
            int pointsCounter = 0;
            while(pointsCounter < n){
                double rUniform = RN.nextDouble(); // random uniform number 0-1
                for(double[] prob: CorrectionProbs){
                    if(rUniform <= prob[3]){
                        int pointx = site[0] + ((int) prob[0]);
                        int pointz = site[2] + ((int) prob[2]);

                        // Ensure equal number of samples at each injection needle
                        // check if sample falls out of range of grid, if so then sample another
                        if(pointx > xSize || pointx < 0 || pointz > xSize || pointz < 0){
                            break;
                        }

                        int[] pointArray = {pointx,0,pointz,site[3]};
                        // Since we are sampling with replacement, we don't want to sample the same cell twice
                        // check if this one already pulled, if so then sample another

                        boolean alreadyExists = false;
                        for(int[] existingPoint : CorrectionPoints){
                            if(Arrays.equals(existingPoint, pointArray)){
                                alreadyExists = true;
                                break;
                            }
                        }
                        if(alreadyExists){
                            break;
                        }

                        CorrectionPoints.add(pointArray);
                        pointsCounter++;
                        break;
                    }
                }
            }
        }

        // Iterate through all cells and see if they are corrected or not based on their correction probability.
        // Adjust generic coordinates based on each injection site
//        for (int[] site:InjectionSites) {
//            for(double[] point:CorrectionProbs){
//                if(RN.nextDouble() < point[3]){
//                    int pointx = site[0] + ((int) point[0]);
//                    int pointz = site[2] + ((int) point[2]);
//                    CorrectionPoints.add(new int[]{pointx, 0, pointz,site[3]});
//                }
//
//            }
//        }
    }

    public void GenerateCorrectionPoints(){ // Adds as an int array coordinates of the desired correction locations. 30Aug23HLC
        int planeSize = xSize;
        int maxCoordinate = microNeedles-1;
        int margin = (int) (planeSize*margin_ratio);
        int distanceX, distanceZ;
        if(microNeedles == 1){
            int centerX = (planeSize-1)/2;
            int centerZ = (planeSize-1)/2;
            InjectionSites.add(new int[]{centerX,0,centerZ,0});
            return;
        }
        if(needleSpacing != null){ // set spacing if predefined.
            distanceX = needleSpacing;
            distanceZ = needleSpacing;
        }
        else{ // Calculate default spacing based on available space and microneedle count.
            distanceX = ((planeSize-1) - 2 * margin) / maxCoordinate;
            distanceZ = ((planeSize-1) - 2 * margin) / maxCoordinate;
        }
        // Calculate initial offset to center the points
        int offsetX = (planeSize - (microNeedles - 1) * distanceX) / 2;
        int offsetZ = (planeSize - (microNeedles - 1) * distanceZ) / 2;

        // Check if any points exceed the specified margins
        boolean xExceedsMargin = false;
        boolean zExceedsMargin = false;
        for (int i=0; i<microNeedles; i++){
            if (offsetX + i * distanceX < margin || offsetX + i * distanceX > planeSize - margin){
                xExceedsMargin = true;
                break;
            }
        }
        for (int j=0; j<microNeedles; j++){
            if (offsetZ + j * distanceZ < margin || offsetZ + j * distanceZ > planeSize - margin){
                zExceedsMargin = true;
            }
        }

        if (xExceedsMargin || zExceedsMargin){
            throw new IllegalArgumentException("Microneedle spacing is not possible due to insufficient space within the tissue section.\n");
        }
        for (int i=0; i < microNeedles; i++){
            for (int j=0; j < microNeedles; j++){
                int pointx = offsetX + i * distanceX;
                int pointz = offsetZ + j * distanceZ;
                InjectionSites.add(new int[]{pointx,0,pointz,i*microNeedles+j});
            }
        }

        // print coordinates of correction sites
//        for(int i=0; i< InjectionSites.size(); i++){
//            System.out.print("Point i= " + i + " x = " + InjectionSites.get(i)[0] + " z= " + InjectionSites.get(i)[2] + "\n");
//        }
    }

    public void RunStep() {
        for (int i = 0; i < CHEMICAL_STEPS; i++) {
            ChemicalLoop();
        }
        for (EpidermisCell c: this) {
            c.CellStep();
            MeanProlif(c);
//            MeanDeath();
        }
        corrected = true; // after first step completed - corrected -> true allows for no further correction 29aug23HLC
        popSum+=GetPop();
        CleanShuffInc(RN); // Special Sauce

        Turnover.RecordBasalRate("Death");
        Turnover.RecordBasalRate("Birth");
        Turnover.RecordTissueRate("Birth");
        Turnover.RecordTissueRate("Death");
    }

    public void DrawChemicals(GuiGridVis chemVis, boolean egf, boolean bfgf) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                    if (egf) {
                        chemVis.SetColorHeat(x, y, EGF.GetCurr(x, y, zDim/2) / SOURCE_EGF, "rgb");
                }
            }
        }
    }

    public void ActivityHeatMap(GuiGridVis heatVis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw, int[] MeanLife, String heatColor) {
        for(int i=0; i<MeanLife.length; i++){
            if(MeanLife[i]!=0) {
                heatVis.SetColorHeat(ItoX(i), ItoY(i), MeanLife[i] / (float)EpidermisConst.VisUpdate, heatColor);
            } else {
                heatVis.SetColor(ItoX(i),ItoY(i), 0.0f, 0.0f, 0.0f);
            }
        }
    }

    public void LayerVis(GuiGridVis heatVis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw, int[] MeanLife, String heatColor) {
        int[] MeanLayer = new int[EpidermisConst.ySize];
        for(int i=0; i<MeanLife.length; i++){
            int y = ItoY(i);
            MeanLayer[y] += MeanLife[i];
            MeanLife[i] = 0;
        }
        for(int y = 0; y<MeanLayer.length; y++) {
            float LayerAvg = MeanLayer[y] / (EpidermisConst.xSize * (float) EpidermisConst.VisUpdate);
                if(LayerAvg!=0) {
                    for (int x = 0; x < EpidermisConst.xSize; x++) {
                        heatVis.SetColorHeat(x, y, LayerAvg, heatColor);
                    }
                } else {
                    for (int x = 0; x < EpidermisConst.xSize; x++) {
                        heatVis.SetColor(x, y, 0.0f, 0.0f, 0.0f);
                    }
                }
        }
    }

    public void MeanProlif(EpidermisCell c){
            if(c.Action == DIVIDE){
                MeanProlif[c.Isq()] += 1;
            }

    }

    public float GetMeanAge(EpidermisGrid Epidermis){
        float Age = 0;
        int aliveCells = 0;
        for (EpidermisCell c: this) {
            if(c!=null){
                Age += c.Age();
                aliveCells += 1;
            }
        }
        return Age/aliveCells;
    }

    public void DrawCellActivity(GuiGridVis vis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw) {
        long time = System.currentTimeMillis();
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                    EpidermisCell c = Epidermis.GetAgent(x, y, zDim/2);
                    if (c != null) {
                        CellDraw.DrawCellonGrid(vis, c);
                    } else {
                        CellDraw.DrawEmptyCell(vis, x, y);
                    }
            }
        }
    }

    public void DrawCellPops(GuiGridVis vis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw){
        long time = System.currentTimeMillis();
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                EpidermisCell c = Epidermis.GetAgent(x, y, zDim/2);
                if (c != null) {
                    CellDraw.DrawCellonGridPop(vis, c);
                } else {
                    CellDraw.DrawEmptyCell(vis, x, y);
                }
            }
        }
    }

    public void DrawCellPopsBottom(GuiGridVis vis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw){
        for(int x=0; x < xDim; x++) {
            for(int z=0; z<zDim; z++) {
                EpidermisCell c = Epidermis.GetAgent(x, 0, z);
                if (c != null){
                    CellDraw.DrawCellonGridPopZ(vis, c);
                } else {
                    CellDraw.DrawEmptyCell(vis, x, z);
                }
            }
        }
    }

    public void DrawCellPopsBottomActivity(GuiGridVis vis, EpidermisGrid Epidermis, EpidermisCellVis CellDraw) {
        long time = System.currentTimeMillis();
        for (int x = 0; x < xDim; x++) {
            for (int z = 0; z < zDim; z++) {
                EpidermisCell c = Epidermis.GetAgent(x, 0, z);
                if (c != null) {
                    CellDraw.DrawCellonGrid3D(vis, c);
                } else {
                    CellDraw.DrawEmptyCell(vis, x, z);
                }
            }
        }
    }

    public void rglVisualization(){
        for(int i=0; i < (EpidermisConst.ySize*EpidermisConst.xSize*EpidermisConst.zSize);i++) {
            EpidermisCell c = GetAgent(i);
            if (c != null) {
                String outLine = i + "\t" + c.myGenome.h + "\t" + c.myGenome.s + "\t" + c.myGenome.v + "\t" + 0.8;
                System.out.println(outLine);
            }
        }
    }

    public void EGFrglVisualization(){
        double egfCol = 0;
        for (int x = 0; x < xDim; x++) {
            for (int z = 0; z < zDim; z++) {
                egfCol = 0;
                for (int y = 0; y < yDim; y++) {
                    egfCol += EGF.GetCurr(x,y,z);
                }
                String outLine = x + "\t" + z + "\t" + egfCol;
                System.out.println(outLine);
            }
        }
    }


    // commented out appears to not be used - 25AUG23 HLC
    // Utils.HSBtoRGB error wants 4 doubles, but input is 3 floats
    // Reworked code to allow proper input to HSBtoRGB function. 25AUG23HLC

    public void BuildMathematicaArray(){
        for(int i=0; i < (EpidermisConst.ySize*EpidermisConst.xSize*EpidermisConst.zSize);i++){
            EpidermisCell c = GetAgent(i);
            if (c != null){
                if(c.myGenome.h ==1.0){ //?
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][0] = 1.0;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][1] = 1.0;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][2] = 1.0;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][3] = 0.1;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = 0;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][5] = -1;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][6] = -1;
                } else {
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][0] = Utils.GetHSBtoRGB(c.myGenome.h,c.myGenome.s,c.myGenome.v)[0];
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][1] = Utils.GetHSBtoRGB(c.myGenome.h,c.myGenome.s,c.myGenome.v)[1];
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][2] = Utils.GetHSBtoRGB(c.myGenome.h,c.myGenome.s,c.myGenome.v)[2];
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][3] = 0.80;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][5] = c.myGenome.injSite;
                    ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][6] = c.tp53Clone;
                    //String genomeString = c.myGenome.GenomeInfoStr();
                    String genomeString = c.myGenome.FullLineageInfoStr("");
//                    if(genomeString.length() >0){
//                        System.out.print(genomeString + " : ");
//                        System.out.print(lineageString + "\n");
//                    }

                    if(genomeString.contains("FAcorrection")){
                        //System.out.print(mutString);
                        //System.out.println(genomeString);
                        ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = 666; // arbitrary number (>71) to mark corrected clones 31OCT23HLC
                        if (genomeString.contains(".68.")){
                            ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = 777; // arbitrary number (>71) to mark corrected + tp53 clones 08FEB24HLC
                        }

                    }
                    else {
                        String[] parts = c.myGenome.GenomeInfoStr().split("\\.");
                        if (parts.length >= 3) {
                            // Extract the numbers between the first two "."
                            Double result = Double.parseDouble(parts[1]);
                            ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = result;
                            if(genomeString.contains(".68.")){
                                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = 68; // keep track of tp53 mutants 08Feb24HLC
                            }
                        }
                        else{
                            ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = 0; // reference cells
                        }
                    }
//                    if(ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] == 68){
//                        if(ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][0] == 1.0 && ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][1] ==1.0){
//                            if(c.Ysq() == 0) {
//                                System.out.print("\n" + ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][0] + "\n");
//                                System.out.print("\nWeird case - Timestep: " + GetTick() + " X:" + c.Xsq() + " Y: " + c.Ysq() + " Z: " + c.Zsq() + "\n");
//                                System.out.print("GenomeInfo String: " + c.myGenome.GenomeInfoStr() + " Full Lineage String: " + genomeString + "\n");
//                            }
//                        }
//                    }
                }
            } else { // note: null (dead) cells are returned as 0,0,0,0,-1 - 29AUG23HLC
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][0] = 0.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][1] = 0.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][2] = 0.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][3] = 0.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] = -1.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][5] = -1.0;
                ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][6] = -1.0;
            }
            if(ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][4] == 68 && c.Ysq() == 0){
                if(ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][6] == -1) {
                    System.out.print(">>>>>> tp53 mutant clone: " + ImageArray[ItoY(i)][ItoX(i)][ItoZ(i)][6] + ":\n");
                    System.out.print(c.myGenome.FullLineageInfoStr("") + "\n");
                }
                }
        }
    }

    // Inflicting a wound to simulate wound repair...
    public void inflict_wound(){
        for (int x = xDim/5; x < (xDim/5)*4; x++) {
            for (int z = zDim/5; z < (zDim/5)*4; z++) {
                for (int y = 0; y < yDim; y++) {
                    EpidermisCell c = GetAgent(x,y,z);
                    if(c!=null){
                        c.itDead();
                    }
                }
            }
        }
//        for(int i=0; i < (EpidermisConst.ySize*EpidermisConst.xSize*EpidermisConst.zSize);i++) {
//            EpidermisCell c = GetAgent(i);
//            if(c!=null){
//                c.itDead();
//            }
//        }
    }

    public void DamageTissueWithUV(double FractionOfDeadCells){
        for (int i = 0; i < (xDim*zDim*yDim); i++) {
            EpidermisCell c = GetAgent(i);
            if(c!=null) {
                String thisGenome = c.myGenome.GenomeInfoStr();
                // If P53 Mutation present standard death function
                if(!thisGenome.contains(".68.") && RN.nextDouble()<FractionOfDeadCells) {
                    c.itDead();
                }
            }
        }
    }

    public boolean checkWoundHeal(int AvgHeight){
        if(GetAgent(xDim/2, 0, zDim/2)!=null){
            return true;
        } else {
            return false;
        }
    }

    public double GetMeanCellHeight(){
        int allColumns = 0;
        for (int x = 0; x < xDim; x++) {
            for (int z = 0; z < zDim; z++) {
                int column = 0;
                for (int y = 0; y < yDim; y++) {
                    EpidermisCell c = GetAgent(x, y, z);
                    if (c != null) {
                        column++;
                    }
                }
                allColumns += column;
            }
        }
        return (allColumns*1.0)/(xDim*zDim);
    }

    public void GetCellPositions(FileIO PositionOut){
        for(int i=0; i < (EpidermisConst.ySize*EpidermisConst.xSize*EpidermisConst.zSize);i++) {
            EpidermisCell c = GetAgent(i);
            if(c!=null){
                String OutString = ItoX(i) + "," + ItoY(i) + "," + ItoZ(i) + "," + c.myGenome.IDGetter() + "\n";
                PositionOut.Write(OutString);
            }
        }
    }

    public String GetDivisionProportion(){
        double[] OutProportions = new double[yDim];
        for (int i = 0; i < (ModelTime-1)*yDim; i+=yDim) {
            for (int y = 0; y < yDim; y++) {
                OutProportions[y]+=divisions[i+y];
            }
        }
        StringBuilder OutNums = new StringBuilder();
        for (int y = 0; y < yDim; y++) {
            String OutNess=OutProportions[y]/divs + "\t";
            OutNums.append(OutNess);
        }
        return OutNums.toString();
    }

    public void ChemicalLoop(){
        //DIFFUSION
        EGF.Diffuse(EGF_DIFFUSION_RATE);
        //CELL CONSUMPTION
        for (EpidermisCell c: this) {
                EGF.AddNext(c.Xsq(),c.Ysq(),c.Zsq(), c.KERATINO_EGF_CONSPUMPTION*EGF.GetCurr(c.Xsq(), c.Ysq(), c.Zsq()));
        }

        //DECAY RATE
        for(int i=0;i<EGF.length;i++){
            EGF.SetNext(ItoX(i),ItoY(i),ItoZ(i), EGF.GetNext(ItoX(i), ItoY(i), ItoZ(i))*(1.0-DECAY_RATE));
        }

        //SOURCE ADDITION
        for(int x=0;x<xDim;x++) {
            for(int z=0;z<zDim;z++) {
                EGF.SetNext(x, 0, z, SOURCE_EGF);
            }
        }

        //SWAP CURRENT FOR NEXT
        EGF.SwapNextCurr();
    }

    public void GetEGFVal(){
        StringBuilder EGFCons = new StringBuilder();
        for (int y=0; y < yDim; y++) {
            String out = String.valueOf(EGF.GetCurr(xDim/2, y, zDim/2)) + "\t";
            EGFCons.append(out);
        }
        System.out.println(EGFCons.toString());
    }

}
