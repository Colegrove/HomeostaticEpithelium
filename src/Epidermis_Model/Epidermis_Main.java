package Epidermis_Model;

import Framework.Gui.GuiGridVis;
import Framework.Gui.GuiLabel;
import Framework.Gui.GuiWindow;
import Framework.Tools.FileIO;
import Framework.Tools.*;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

import static Epidermis_Model.EpidermisConst.*;


/**
 * Created by schencro on 3/24/17.
 */

//Holds Constants for rest of model
class EpidermisConst {
    static int xSize = 20; // CHANGE

    //static int xSize = 100; // keratinocyte modal cell size = 15µm (Proc. Natl. Acad. Sci. USA Vol.82,pp.5390-5394,August1985; YANN BARRANDON and HOWARD GREEN) == volume == 1766.25µm^3

    // (Sampled area = 1mm-2mm^2); Sampled volume = 4.4*10^8µm^3; Total cells needed for 2mm^2 area with depth of 140µm= 249115cells (xSize = 12456, ySize = 20);
    // For 1mm^2 area with depth of 140µm = 62279cells (xSize = 3114, ySize = 20);
    // Takes forever to reach even a year. Cutting the smallest biopsy into a quarter (1/4) = 15570cells (xSize = 1038, ySize = 20)
    // Above numbers are for 2D, for 3D the xSize = 100
    static final int ySize = 20;
    static int zSize = xSize;

    static final int KERATINOCYTE = 0; //setting types into a binary 0 or 1
    static final int DIVIDE = 2; // Attribute if cell is dividing
    static final int STATIONARY = 3; // Attribute if cell is stationary
    static final int MOVING = 4; //Attribute if cell is moving

    static int years = 50; // time in years.
    //static int RecordTime = years * 365;
    static int RecordTime = 365;

    //static int[] RecordTimeArray = new int[120];
    static int[] RecordTimeArray = new int[1001];

//    static int[] RecordTimeArray = {
//            2, 365, 730, 1095, 1460, 1825, 2190, 2555, 2920, 3285, 3650, 4015, 4380, 4745,
//            5110, 5475, 5840, 6205, 6570, 6935, 7300
//    };

    // record time scaled by 7 timesteps per day and by 4.5 timesteps per day
    // record 2 days, 3 weeks, 6 weeks, and 12 weeks
//    static int[] RecordTimeArray = {
//            2, 95, 189, 378, 147, 294, 588
//
//    };
    static int ModelTime = years * 365 + 10; // Time in days + 10 days after time for recording! e.v. 65 years = 23725
    //static int ModelTime = 50; // Time in days + 10 days after time for recording! e.v. 65 years = 23725


    static final int VisUpdate = 7; // Timestep interval to update Division and Death, etc.
    static int Replicate = 1; // Replicate number to be multiplied by the RecordTime to set the seed

    /**
     * Correction parameters 30Aug23HLC
     */

    static boolean stopAtConfluence = false; // determines if the model should stop running if confluence is reached.
    static boolean stopAtLoss = false; // determines if the model should stop running if corrected cells are lost.
    static int stopLossTime = 0; // 0 will stop at any time of loss, 1 will only stop simulations lost within 1 year.
    static boolean confluence = false; // if model reaches confluence, boolean changes true. 14Nov23HLC
    static boolean loss = false; // if model loses all corrected cells, boolean changes true. 14Nov23HLC
    static int correctionTime = 10; // Timestep of correction (days)
    static int microNeedles = 1; // Microneedle injection sites: square value for total number of needles
    static double margin_ratio = 0.15; // Margins on edge of tissue where microneedles are not injected.
    static Integer needleSpacing = 50; // A defined distance between adjacent microneedles in cellular distance (null for default spacing)
    static boolean CorrectedBlockChanges = true; // Whether to run corrected cell's selective advantage
    static boolean allowVerticalDivision = true; // If blocking probability prevents division, then divide up.
    static double CorrectedBlockProbability = 1.0; // Corrected cell's blocking probability
    static boolean CorrectedGrowthChanges = false; // Whether to run corrected cell's selective advantage
    static double CorrectedGrowthIncrease = 1.0; // Acts on proliferation scale factor. e.g. value of 1.1 provides a 10% increase to proliferation scale factor
    //static double tp53BlockProb = 0.41; // Gives tp53 a blocking probability advantage compared to background FA
    static double tp53BlockProb = 0.01; // Gives tp53 a blocking probability advantage when it's on a corrected mutant background
    static double FAtp53BlockProb = 0.01; // Gives tp53 a blocking probabiliy advantage when it's on a FA mutant background
    static int transGeneDiffusion = 2; // Diagonals of bivariate normal distribution covariance matrix acting on the diffusion of transgene.
    static int dose = 30; // How many cells to be corrected at each microneedle
    //static String corrProbFile = "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/corrProbs/corrProbs_sigma_" + transGeneDiffusion + ".txt";
    static String corrProbFile = "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/corrProbs/corrProbs_sigma_" + transGeneDiffusion + ".txt";
    static boolean divRate = false; // if true output all division information into a file.
    static int timepoint_collection = 0; // 0 for every six months, 1 for every year, 2 for every timestep (1000), 3 for martincorena timepoints, 4 cell division

    /**
     * Mutation Rate Set Options
     */
    static int MutRateSet = 5; // Select which mutation rate is required.
    static double mutRateAdjustment = 1; // reduces the mutation rate by a factor of
    static int corrMutRate = 2; // Corrected cells reduce mutation rate by factor of

    /**
     * TP53 Fitness Change Options
     */
    static final int SunDays = 7; // Number of days with high sun exposure in a year.
    static int SunDaysFreqency = 12; // Number of days between sun exposures within a year.
    static double SunDaysDeathProb = 0.1; // Fraction of cells that die.
    static final boolean PrintPopsForODE = false;
    static final boolean PrintSunDays = false;

    static double tp53GrowthIncreaseFA = 1.1; // Acts on proliferation scale factor on FA background. e.g. value of 1.1 provides a 10% increase to proliferation scale factor
    static double tp53GrowthIncrease = 1.01; // Acts on proliferation scale factor on normal (FA corrected) background.
    static int tp53CloneTracker = -1; // counts tp53 clones


    /**
     * NOTCH Fitness Change Options
     */
    static double NOTCHBlockProbability=0.0;

    /**
     * Booleans for Run Options
     */
    static final boolean GuiOn = false; // use for visualization, set to false for jar file / multiple runs
    static final boolean JarFile = true; // Set to true if running from command line as jar file!!!!!!!!
    static final boolean RecordParents = true; // use when you want parents information
    static final boolean RecordLineages = true; // use when you want
    static final boolean RecordPopSizes = true; // Use to record clone population sizes
    static final boolean RecordAllPopSizes = true; // Use this to record all population sizes for each time step.
    static final boolean get_r_lambda = true; // use when you want the r_lambda value for the visualization
    static final boolean writeValues = true; // use this when you want the data to be saved!
    static final boolean sliceOnly = false; // use this when you want slice of the 3D model data to be output!!!!!!!!!!!!!!
    static final boolean GetImageData = true; // Use for 3D data for visualization
    static final boolean GetEGFSum = false; // Use for 3D data for visualization of EGF concentrations
    static final boolean Wounding = false; // Use to do wounding
    static final boolean PFiftyThree = false; // Whether to perform P53 Fitness testing through turnind Random Death Prob off.
    static final boolean PFiftyThreeSunDays = false; // Whether to include a sun days UV damage rate.
    static final boolean p53Growth = false; // 27Oct23 HLC gives p53 mutations a growth rate advantage.
    static final boolean p53Blocking = true;
    static final boolean NOTCH1FitnessChanges = false; // Whether to run NOTCH1 Fitness Changes.

}

public class Epidermis_Main {

    static GuiLabel LabelGuiSet(String text, int compX, int compY) {
        GuiLabel ret= new GuiLabel(text, compX, compY);
        ret.SetColor(Color.white,Color.black);
        return ret;
    }

    public static void main (String[] args){
        /*
        Initialization
         */
        GuiWindow MainGUI=null;
        GuiGridVis ActivityVis = null;
        GuiGridVis EGFVis = null;
        GuiGridVis DivVis = null;
        GuiGridVis DivLayerVis = null;
        GuiGridVis DeathVis = null;
        GuiGridVis DeathLayerVis = null;
        GuiGridVis ClonalVis = null;
        GuiGridVis BottomVis = null;
        GuiGridVis BottomVisMove = null;
        GuiLabel YearLab = null;
        GuiLabel rLambda_Label = null;
        GuiLabel OldestCell = null;
        GuiLabel HealLab = null;
        GuiLabel HeightLab = null;
        GuiLabel NullLabel = null;
        EpidermisCellVis CellDraw = null;
        ArrayList<Float> r_lambda_WriteValue = new ArrayList();
        int r_lambda_index = 0;
        ArrayList<Float> meanCellAge = new ArrayList();
        int meanCellAgeIndex = 0;

        String ParentFile = System.getProperty("user.dir") + "/TestOutput/ParentFile.csv";
        String PopSizes = System.getProperty("user.dir") + "/TestOutput/PopSizes.csv";
        String MutationFile = System.getProperty("user.dir") + "/TestOutput/MutationFile.csv";
        String r_lambda_file = System.getProperty("user.dir") + "/TestOutput/R_Lambda_Values.csv";
        String PositionFile = System.getProperty("user.dir") + "/TestOutput/PositionList.csv";
        String Image_file = System.getProperty("user.dir") + "/TestOutput/VisFiles/VisFile.txt";
        String divisionsFile = System.getProperty("user.dir") + "/TestOutput/divFile.txt";


        /*
        Sets up Data Files if on cluster or if ran locally
         */
        if(EpidermisConst.JarFile){
            ParentFile = args[0];
            PopSizes = args[1];
            MutationFile = args[2];
            r_lambda_file = args[3];
            Image_file = args[4]; // 12SEP23HLC
            EpidermisConst.xSize = Integer.parseInt(args[5]);
            EpidermisConst.zSize = Integer.parseInt(args[5]);
            int Time = Integer.parseInt(args[6]);
            EpidermisConst.years = Time;
            EpidermisConst.ModelTime = (int) (Time * 364 * 4.5);
            //EpidermisConst.ModelTime = 590;
            EpidermisConst.RecordTime = Time * 364;
            EpidermisConst.MutRateSet = Integer.parseInt(args[7]);
            //EpidermisCellGenome.MutRateSet = EpidermisConst.MutRateSet;
            //EpidermisConst.SunDaysFreqency = Integer.parseInt(args[7]);
            //EpidermisConst.SunDaysDeathProb = Double.parseDouble(args[8]);
            //EpidermisConst.Replicate = Integer.parseInt(args[9]);
            //EpidermisConst.NOTCHBlockProbability = Double.parseDouble(args[10]);
            EpidermisConst.microNeedles = Integer.parseInt(args[8]); // 12SEP23HLC
            EpidermisConst.needleSpacing = Integer.parseInt(args[9]); // 12SEP23HLC
            EpidermisConst.CorrectedBlockProbability = Double.parseDouble(args[10]); // 12SEP23HLC
            EpidermisConst.CorrectedGrowthIncrease = Double.parseDouble(args[11]); // 03OCT23HLC
            EpidermisConst.CorrectedBlockChanges = Boolean.parseBoolean(args[12]); //16OCT23HLC
            EpidermisConst.allowVerticalDivision = Boolean.parseBoolean(args[13]); // 16OCT23HLC
            EpidermisConst.CorrectedGrowthChanges = Boolean.parseBoolean(args[14]); //16OCT23HLC
            EpidermisConst.correctionTime = Integer.parseInt(args[15]); // 27OCT23 HLC
            EpidermisConst.tp53GrowthIncreaseFA = Double.parseDouble(args[16]); // 30Oct23HLC
            EpidermisConst.tp53GrowthIncrease = Double.parseDouble(args[17]); // 06NOV23HLC
            EpidermisConst.tp53BlockProb = Double.parseDouble(args[18]);
            EpidermisConst.transGeneDiffusion = Integer.parseInt(args[19]); // 13NOV23HLC
            corrProbFile = "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/corrProbs/corrProbs_sigma_" + transGeneDiffusion + ".txt";
            EpidermisConst.dose = Integer.parseInt(args[20]); // 13NOV23HLC
            stopAtConfluence = Boolean.parseBoolean(args[21]); // 10JAN24HLC
            stopAtLoss = Boolean.parseBoolean(args[22]); //10JAN24HLC
            EpidermisConst.mutRateAdjustment = Double.parseDouble(args[23]);
            EpidermisCellGenome.mutRateAdjustment = EpidermisConst.mutRateAdjustment;
            timepoint_collection = Integer.parseInt(args[24]);
            stopLossTime = Integer.parseInt(args[25]);
            EpidermisConst.corrMutRate = Integer.parseInt(args[26]);
            EpidermisCellGenome.corrMutRate = EpidermisConst.corrMutRate;
            EpidermisConst.FAtp53BlockProb = Double.parseDouble(args[27]);

//            PositionFile = args[7];
        } else {
            EpidermisCellGenome.MutRateSet = EpidermisConst.MutRateSet;
        }
        // update mutation rates if needed
        EpidermisCellGenome.PoissonDists_corr = EpidermisCellGenome.BuildPoissons_corr(corrMutRate);
        EpidermisCellGenome.PoissonDists = EpidermisCellGenome.BuildPoissons(mutRateAdjustment);

        if(EpidermisConst.GuiOn == false && EpidermisConst.GetImageData == false){
            System.out.println("xSize and zSize: " + EpidermisConst.xSize);
            System.out.println("Years: " + EpidermisConst.years);
        }

        // Record Time array to record every 6 months with 4.5 timestep scaling.
        if(timepoint_collection == 0){
            double rawValue = 182*4.5;
            // every 6 months
            EpidermisConst.RecordTimeArray[0] = 2;
            for(int i=1; i<=EpidermisConst.years*2; i++){
                EpidermisConst.RecordTimeArray[i] = (int) rawValue*i;
            }
        } else if (timepoint_collection == 1) {
             //every year
            double rawValue = 364*4.5;
            EpidermisConst.RecordTimeArray[0] = 2;
            for(int i=1; i<=EpidermisConst.years; i++){
                EpidermisConst.RecordTimeArray[i] = (int) rawValue*i;
            }
        } else if (timepoint_collection == 2){
            // every timestep 1000 frames
            for(int i=1; i<=1000; i++){
                EpidermisConst.RecordTimeArray[i] = i;
            }
        } else if (timepoint_collection == 3){
            // martincorena timepoints only
            EpidermisConst.RecordTimeArray[0] = 35217;
            EpidermisConst.RecordTimeArray[1] = 41769;
            EpidermisConst.RecordTimeArray[2] = 61425;
            EpidermisConst.RecordTimeArray[3] = 74529;
            EpidermisConst.RecordTimeArray[4] = 81081;
            EpidermisConst.RecordTimeArray[5] = 87633;
        } else if (timepoint_collection == 4){
            // every timepoint
            ModelTime = 100;
        }


        for(int i=0; i< RecordTimeArray.length; i++){
            System.out.print("Time Array: " + RecordTimeArray[i] + " Final run time: " + ModelTime + "\n");
        }
        System.out.print("MutrateAdjustment: " + mutRateAdjustment + "\n");
        System.out.print("PossionDists[3] " + EpidermisCellGenome.PoissonDists[3] + "\n");
        System.out.print("Corr MutrateAdjustment: " + corrMutRate + "\n");
        System.out.print("PossionDists[3] " + EpidermisCellGenome.PoissonDists_corr[3] + "\n");

        final EpidermisGrid Epidermis = new EpidermisGrid(EpidermisConst.xSize, EpidermisConst.ySize, EpidermisConst.zSize); // Initializes and sets up the program for running
        Runtime rt = Runtime.getRuntime();

        // Sets up GUI
        if(EpidermisConst.GuiOn) {
            CellDraw = new EpidermisCellVis();
            MainGUI = new GuiWindow("Homeostatic Epidermis Model", true);
            MainGUI.panel.setOpaque(true);
            MainGUI.panel.setBackground(Color.black);
            ClonalVis = new GuiGridVis(EpidermisConst.xSize*5, EpidermisConst.ySize*5, 1, 2, 1);
            DivVis = new GuiGridVis(EpidermisConst.xSize, EpidermisConst.ySize, 3, 1, 1);
            DivLayerVis = new GuiGridVis(EpidermisConst.xSize, EpidermisConst.ySize, 3, 1, 1);
            DeathVis = new GuiGridVis(EpidermisConst.xSize, EpidermisConst.ySize, 3, 1, 1);
            DeathLayerVis = new GuiGridVis(EpidermisConst.xSize, EpidermisConst.ySize, 3, 1, 1);
            ActivityVis = new GuiGridVis(EpidermisConst.xSize * 5, EpidermisConst.ySize * 5, 1, 2, 1); // Main Epidermis visualization window
            BottomVis = new GuiGridVis(EpidermisConst.xSize*6, EpidermisConst.zSize*6, 1,1,1);
            BottomVisMove = new GuiGridVis(EpidermisConst.xSize*6, EpidermisConst.zSize*6, 1,1,1);
            EGFVis = new GuiGridVis(EpidermisConst.xSize, EpidermisConst.ySize, 5, 2, 1);
            YearLab = LabelGuiSet("Age (Yrs.): ", 1, 1);
            MainGUI.AddCol(YearLab, 0);
            HealLab = LabelGuiSet("Heal Time (Days): ", 1, 1);
            MainGUI.AddCol(HealLab, 0);
            HeightLab = LabelGuiSet("Height: ", 1, 1);
            MainGUI.AddCol(HeightLab, 0);
            OldestCell = LabelGuiSet("Oldest Cell: ", 1, 1);
            rLambda_Label = LabelGuiSet("rLambda: ", 1, 1);
            MainGUI.AddCol(rLambda_Label, 1);
            MainGUI.AddCol(OldestCell, 1);
            NullLabel = LabelGuiSet(" ", 1, 1);
            MainGUI.AddCol(NullLabel,1);
            MainGUI.AddCol(BottomVis, 0);
            MainGUI.AddCol(BottomVisMove, 1);
            MainGUI.AddCol(LabelGuiSet("Population", 2, 1), 0);
            MainGUI.AddCol(ClonalVis, 0);
            MainGUI.AddCol(LabelGuiSet("Division (per week)", 1, 1), 0);
            MainGUI.AddCol(DivVis, 0);
            MainGUI.AddCol(LabelGuiSet("Division Layers (per week)", 1, 1), 1);
            MainGUI.AddCol(LabelGuiSet("Death (per week)", 1, 1), 0);
            MainGUI.AddCol(DeathVis,0);
            MainGUI.AddCol(DivLayerVis, 1);
            MainGUI.AddCol(LabelGuiSet("Death Layer (per week)", 1, 1), 1);
            MainGUI.AddCol(DeathLayerVis, 1);
            MainGUI.AddCol(LabelGuiSet("Epidermis", 2, 1), 0);
            MainGUI.AddCol(ActivityVis, 0); // Main Epidermis visualization window
            MainGUI.AddCol(LabelGuiSet("EGF", 2, 1), 0);
            MainGUI.AddCol(EGFVis, 0);

            MainGUI.RunGui();
        }
        int woundTick = 0;
        boolean Healed = true;
        double avgHeight=0;
        int tickSum=0;
        int wounded=0;
        int SunDayCounter=0;
        int[] SunTimes=new int[EpidermisConst.SunDays];

        // method to generate a list of coordinates for correction based on microneedles and spacing 29aug23HLC
        Epidermis.GenerateCorrectionPoints();
        // method to generate a list of cells to be corrected after transgene diffusion from above needle injection sites 14Nov23HLC
        Epidermis.diffusionCorrectionPoints();

        TickRateTimer tickIt = new TickRateTimer();
        while(Epidermis.GetTick() < EpidermisConst.ModelTime){
            //System.out.print("Model is starting timestep: " + Epidermis.GetTick() + "\n");

//            tickIt.TickPause(60); // Adjusting a frame rate

            if(Epidermis.GetTick() == EpidermisConst.correctionTime){ // correction timestep - triggers boolean for correction 29aug23HLC
                EpidermisGrid.corrected = false;
            }

            // Main Running of the steps within the model
            Epidermis.RunStep();

            /*
            All Injuries Occuring Here!
             */
            if(EpidermisConst.Wounding) {
                int healTick = 0;
                if (Healed && Epidermis.GetTick() % 100 == 0 && wounded < 1) {
                    Epidermis.inflict_wound();
                    woundTick = Epidermis.GetTick();
                    Healed = false;
                    wounded++;
                }
            }

//            if(!Healed && Epidermis.GetTick()%50!=0) {
//                Healed = Epidermis.checkWoundHeal((int) avgHeight);
//                healTick = Epidermis.GetTick();
////                if (Healed && HealLab != null) {
////                    if (HealLab != null) {
////                        HealLab.setText("Heal Time (Days): " + new DecimalFormat("#.0").format((healTick - woundTick)));
////                    }
////                }
//            }
            /*
            Get the Diffusion Values for examining 2D versus 3D differences
             */
//            if(Epidermis.GetTick()<=365){Epidermis.GetEGFVal();}

            /*
            Output Time Options
             */
            if(ActivityVis==null){
                if(Epidermis.GetTick()%365==0){
//                    System.out.println(new DecimalFormat("#.0").format((Epidermis.GetTick() / 365f)));
                }
            }

            if(EpidermisConst.PrintPopsForODE){
                System.out.println(Epidermis.GetPop());

            }

//            System.out.println(Epidermis.Turnover.GetBasalRate("Death",Epidermis.GetTick()));
//            if(Epidermis.GetTick()==EpidermisConst.ModelTime-1){
//            System.out.println(Epidermis.GetDivisionProportion());
//            }

            /*
            All Visualization Components are here
             */
            if(Epidermis.GetTick()%7==0){
                if(rLambda_Label!=null){rLambda_Label.SetText("µrLambda(w): " + new DecimalFormat("#.000").format( Epidermis.Turnover.GetBasalRate("Death",7) ));}
                if(HeightLab!=null){HeightLab.SetText("H: " + new DecimalFormat("#.00").format(Epidermis.GetMeanCellHeight()));}
            }
            if(ActivityVis!=null){YearLab.SetText("Age(Y): " + new DecimalFormat("#.00").format((Epidermis.GetTick() / 365f)));}
            if(DivVis!=null&Epidermis.GetTick()%EpidermisConst.VisUpdate==0){Epidermis.ActivityHeatMap(DivVis, Epidermis, CellDraw, Epidermis.MeanProlif, "gbr");}
            if(DivLayerVis!=null&Epidermis.GetTick()%EpidermisConst.VisUpdate==0){Epidermis.LayerVis(DivLayerVis, Epidermis, CellDraw, Epidermis.MeanProlif, "gbr");}
            if(DeathVis!=null&Epidermis.GetTick()%EpidermisConst.VisUpdate==0){Epidermis.ActivityHeatMap(DeathVis, Epidermis, CellDraw, Epidermis.MeanDeath, "rbg");}
            if(DeathLayerVis!=null&Epidermis.GetTick()%EpidermisConst.VisUpdate==0){Epidermis.LayerVis(DeathLayerVis, Epidermis, CellDraw, Epidermis.MeanDeath, "rbg");}
            if(ClonalVis!=null){Epidermis.DrawCellPops(ClonalVis, Epidermis, CellDraw);} // 3D Good
            if(OldestCell!=null){OldestCell.SetText("µCellAge: " + new DecimalFormat("#.00").format(Epidermis.GetMeanAge(Epidermis)));}
            if(ActivityVis!=null){Epidermis.DrawCellActivity(ActivityVis, Epidermis, CellDraw);}
            if(BottomVis!=null){Epidermis.DrawCellPopsBottom(BottomVis, Epidermis, CellDraw);}
            if(BottomVisMove!=null){Epidermis.DrawCellPopsBottomActivity(BottomVisMove, Epidermis, CellDraw);}
            if(EGFVis!=null){Epidermis.DrawChemicals(EGFVis, true, false);} // 3D Good

            // check if corrected cells have reached confluence (or loss) of the basal layer. 14Nov23HLC
            Epidermis.confluenceCheck();

            for (int i=0; i< EpidermisConst.RecordTimeArray.length; i++) {
                int Time = EpidermisConst.RecordTimeArray[i];
                // Use this to get the information for 3D visualizations for OpenGL
                if (EpidermisConst.GetImageData && (Epidermis.GetTick() == Time || EpidermisConst.confluence || EpidermisConst.loss)) {
                    Epidermis.BuildMathematicaArray();
                    FileIO VisOut = new FileIO(Image_file + "." + Epidermis.GetTick() + ".txt", "w");
                    for (int x = 0; x < EpidermisConst.xSize; x++) {
                        for (int y = 0; y < EpidermisConst.ySize; y++) {
                            for (int z = 0; z < EpidermisConst.zSize; z++) {
                                //if (Epidermis.ImageArray[y][x][z][0] != 0.0f && Epidermis.ImageArray[y][x][z][1] != 0.0f && Epidermis.ImageArray[y][x][z][2] != 0.0f && Epidermis.ImageArray[y][x][z][3] != 0.0f){
                                String outLine =
                                        x + "\t" + z + "\t" + y + "\t" +
                                                Epidermis.ImageArray[y][x][z][0] + "\t" + Epidermis.ImageArray[y][x][z][1] +
                                                "\t" + Epidermis.ImageArray[y][x][z][2] + "\t" + Epidermis.ImageArray[y][x][z][3] +
                                                "\t" + Epidermis.ImageArray[y][x][z][4] + "\t" + Epidermis.ImageArray[y][x][z][5] +
                                                "\t" + Epidermis.ImageArray[y][x][z][6] + "\n";
                                VisOut.Write(outLine);
//                                if(Epidermis.ImageArray[y][x][z][4] == 68){
//                                    System.out.print(">>>>>> tp53 mutant clone: " + Epidermis.ImageArray[y][x][z][6] + "\n");
//                                    System.out.print()
//                                }
                            }
                            //}
                        }
                    }
                    if(EpidermisConst.confluence){
                        VisOut.Write("Confluence");
                    }
                    if(EpidermisConst.loss){
                        VisOut.Write("Loss");
                    }
                    VisOut.Close();
                    System.out.println("Vis timepoint saved. Year: " + Epidermis.GetTick()/4.5/364);
                }
                if(EpidermisConst.confluence || EpidermisConst.loss){
                    break;
                }
            }


//            if(EpidermisConst.GetImageData==true && (Epidermis.GetTick() / 365f == 25 || Epidermis.GetTick() / 365f == 50 || Epidermis.GetTick() / 365f == 75)){
//                System.out.println(new DecimalFormat("#.0").format((Epidermis.GetTick() / 365f)));
//                Epidermis.rglVisualization();
//            }

            if(EpidermisConst.Wounding) {
                if (EpidermisConst.GetImageData == true && (Epidermis.GetTick() % 25f == 0)) {
                    System.out.println(new DecimalFormat("#.0").format((Epidermis.GetTick() / 365f)));
                    Epidermis.rglVisualization();
                }

                if (EpidermisConst.GetEGFSum == true && (Epidermis.GetTick() % 25f == 0)) {
                    System.out.println(new DecimalFormat("#.0").format((Epidermis.GetTick() / 365f)));
                    Epidermis.EGFrglVisualization();
                }
            }


            // Use this to get the information for 3D visualizations
//            if(EpidermisConst.GetImageData && EpidermisConst.RecordTime == Epidermis.GetTick()){
//                Epidermis.BuildMathematicaArray();
//                FileIO VisOut = new FileIO(Image_file, "w");
//                String open="{\n";
//                String closer="}\n";
//                for(int y=EpidermisConst.ySize-1; y >= 0;y--){
//                    for(int x=0; x < EpidermisConst.xSize;x++){
//                        for(int z=0; z < EpidermisConst.zSize;z++){
//                            String outLine =
//                                    Epidermis.ImageArray[y][x][z][0] + "\t" + Epidermis.ImageArray[y][x][z][1] +
//                                    "\t" + Epidermis.ImageArray[y][x][z][2] + "\t" + Epidermis.ImageArray[y][x][z][3] +
//                                    "\n";
//                            VisOut.Write(outLine);
//                        }
//                    }
//                }
//
//                VisOut.Close();
//                /* Use this code snippit to get the threeD vis on mathematica
//                file=Import["VisFile(2).txt","Data"]
//                matrix = ArrayReshape[file,{19,14,14,4}]
//                Image3D[matrix, ImageSize->Large,ColorSpace->"RGB", Axes->True,Boxed->False, Method-> {"InterpolateValues" -> False},Background->Black]
//                 */
//            }

            if(EpidermisConst.PFiftyThreeSunDays && Epidermis.GetTick()%365.0==0){
                for (int i = 0; i < EpidermisConst.SunDays; i++) {
                    SunTimes[i] = Epidermis.GetTick() + i * EpidermisConst.SunDaysFreqency;
                }
            }

            if(EpidermisConst.PFiftyThreeSunDays && Epidermis.GetTick()==SunTimes[SunDayCounter]) {
                SunDayCounter++;
                Epidermis.DamageTissueWithUV(EpidermisConst.SunDaysDeathProb);
                if(EpidermisConst.PrintSunDays) {
                    System.out.println(Epidermis.GetTick());
                }
                if (SunDayCounter > EpidermisConst.SunDays-1) {
                    SunDayCounter = 0;
                }
            }

            if(EpidermisConst.RecordAllPopSizes && EpidermisConst.RecordTime != Epidermis.GetTick() && (Epidermis.GetTick()%1.)==0){ //Updated to take every timepoint 20Sep23HLC
                Epidermis.GenomeStore.RecordClonePops();
            }

            /*
            All Model Data Recording Is Below This line
             */
            if(EpidermisConst.writeValues==true) {
                /*
                This section of the code is responsible for recording the full modeled dimensions.
                 */
                if (EpidermisConst.RecordParents == true && EpidermisConst.RecordTime == Epidermis.GetTick()) {
                    FileIO ParentOut = new FileIO(ParentFile, "w");
                    Epidermis.GenomeStore.WriteParentIDs(ParentOut, "\n");
                    ParentOut.Close();
                    System.out.println("Parents written to file.");
                }
                if (EpidermisConst.RecordLineages == true && EpidermisConst.RecordTime == Epidermis.GetTick()) {
                    FileIO MutsOut = new FileIO(MutationFile, "w");
                    Epidermis.GenomeStore.WriteAllLineageInfoLiving(MutsOut, ",", "\n");
                    MutsOut.Close();
                    System.out.println("Lineage genomes written to file.");
                }
                if (EpidermisConst.RecordPopSizes == true && EpidermisConst.RecordTime == Epidermis.GetTick()) {
                    FileIO PopSizeOut = new FileIO(PopSizes, "w");
                    Epidermis.GenomeStore.RecordClonePops();
                    Epidermis.GenomeStore.WriteClonePops(PopSizeOut, ",", "\n");
                    PopSizeOut.Close();
                    System.out.println("Population sizes written to file.");
                }
                if (EpidermisConst.get_r_lambda == true && EpidermisConst.RecordTime == Epidermis.GetTick()) {
                    FileIO RLambdaWriter = new FileIO(r_lambda_file, "w");
                    float r_lamb_print = 0;
                    for (int i = 0; i < Epidermis.Turnover.GetDeathRateBasal().length ; i++) {
                        RLambdaWriter.Write(Epidermis.Turnover.GetDeathRateBasal()[i] + "\n");
                        r_lamb_print+=Epidermis.Turnover.GetDeathRateBasal()[i];
                    }
                    RLambdaWriter.Close();
                    System.out.println("Mean weekly rLambda: " + new DecimalFormat("#.000").format(r_lamb_print / Epidermis.Turnover.GetDeathRateBasal().length) + "\n");
                }
                if (EpidermisConst.get_r_lambda == true && EpidermisConst.RecordTime == Epidermis.GetTick()) {
                    float MeanWeekPrint = 0;
                    for (int i = 0; i < meanCellAge.size(); i++) {
                        MeanWeekPrint += meanCellAge.get(i);
                    }
                }
                if (EpidermisConst.sliceOnly==true && EpidermisConst.RecordTime == Epidermis.GetTick()){
                    FileIO PositionOut = new FileIO(PositionFile, "w");
                    Epidermis.GetCellPositions(PositionOut);
                    PositionOut.Close();
                    System.out.println("Position Information Saved to File");
                }
            }
            // end simulation if corrected cells have reached confluence or loss 14NOV23HLC
            if(EpidermisConst.confluence){
                System.out.print("Corrected cells reached basal layer confluence at day: " + Epidermis.GetTick());
                break;
            }
            if(EpidermisConst.loss){
                System.out.print("Corrected cells lost from basal layer at day: " + Epidermis.GetTick());
                break;
            }
        }
//        if(EpidermisConst.divRate){
//            FileIO divOut = new FileIO(divisionsFile, "w");
//            for (int x = 0; x < EpidermisConst.xSize; x++) {
//                for (int z = 0; z < EpidermisConst.zSize; z++) {
//                    for(int t=0; t < EpidermisConst.ModelTime; t++){
//                        String outLine = x + "\t" + z + "\t" + t + "\t" + EpidermisGrid.divArray[x][z][t] + "\n";
//                        divOut.Write(outLine);
//                    }
//                }
//            }
//            divOut.Close();
//        }




//        System.out.println(java.util.Arrays.toString(EpidermisCell.dipshit));
//        System.out.println(java.util.Arrays.toString(EpidermisCell.dipshitDiv));
//
//        Utils.PrintMemoryUsage();
        if(timepoint_collection == 4) {
            FileIO divTrackerFile = new FileIO(divisionsFile, "w");
            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < xSize; j++) {
                    for (int k = 0; k < ModelTime; k++) {
                        String outLine =
                                i + "\t" + j + "\t" + k + "\t" + Epidermis.divTrack[i][j][k] + "\n";
                        divTrackerFile.Write(outLine);
                        //divTrackerFile.Write(String.valueOf(Epidermis.total_divisions));
                    }
                }
            }

            divTrackerFile.Close();
            System.out.println("divisions written to file.");
        }
    }
}
