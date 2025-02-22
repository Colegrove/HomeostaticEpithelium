package Epidermis_Model;

import Epidermis_Model.Genome.GenomeInfo;
import Framework.Tools.FileIO;
import Framework.Tools.Utils;
import cern.jet.random.Poisson;

import java.util.ArrayList;
import java.util.Random;

import static Epidermis_Model.EpidermisCell.RNEngine;


/**
 * Created by schencro on 3/25/17.
 */


public class EpidermisCellGenome extends GenomeInfo<EpidermisCellGenome> {
    /*
    New Information To Keep Inside the Model!!!!! Official Information
     */

    //private static final String BaseIndexFile= System.getProperty("user.dir") + "/src/Epidermis_Model/Global_Info/BaseIndexes.csv";
    private static final String BaseIndexFile= "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/src/Epidermis_Model/Global_Info/BaseIndexes.csv";
    static final int GenomeComponents = 71;
    static final double HumanGenome = 3200000000.0;
    static final String[] GeneNames = new String[]{"Genome","ADAM29","ADAMTS18","AJUBA","AKT1","AKT2","APOB","ARID1A","ARID2","AURKA","BAI3","BRAF","CASP8","CCND1","CDH1","CDKN2A","CR2","CREBBP","CUL3","DICER1","EGFR","EPHA2","ERBB2","ERBB3","ERBB4","EZH2","FAT1","FAT4","FBXW7","FGFR1","FGFR2","FGFR3","FLG2","GRIN2A","GRM3","HRAS","IRF6","KCNH5","KEAP1","KRAS","MET","MUC17","NF1","NFE2L2","NOTCH1","NOTCH2","NOTCH3","NRAS","NSD1","PCED1B","PIK3CA","PLCB1","PPP1R3A","PREX2","PTCH1","PTEN","PTPRT","RB1","RBM10","SALL1","SCN11A","SCN1A","SETD2","SMAD4","SMO","SOX2","SPHKAP","SUFU","TP53","TP63","TRIOBP"};
    static final int dumb = GeneNames.length;

    static final long[] GeneLengths = new long[]{3191618082l,2463l,3666l,1617l,1443l,1444l,13692l,6205l,5506l,1212l,4569l,2301l,1440l,888l,2649l,514l,3279l,7329l,2307l,5767l,3633l,2931l,3695l,4023l,3927l,2239l,13767l,14944l,2124l,2616l,2614l,2443l,7176l,4391l,2781l,633l,1404l,1976l,1875l,687l,4185l,13482l,8520l,1771l,7668l,7416l,6966l,570l,7351l,1299l,3207l,3750l,3369l,4821l,4337l,1212l,4383l,2787l,2951l,3985l,5376l,6026l,7695l,1659l,2364l,954l,5103l,1455l,1289l,2041l,7098};
    //static final double[] ExpectedMuts = new double[]{1.02E+01,1.2298960944e-05,1.96790865336e-05,5.17789144098534e-06,4.0831751376e-06,3.4550775719999997e-06,5.31298384596e-05,1.09091587995e-05,1.4262587681399999e-05,5.799936312e-06,3.26067900354e-05,7.5642604596e-06,3.1009560480000003e-06,2.632816548e-06,4.9328483301e-06,1.4758253784000002e-06,1.01311941753e-05,1.31902673427e-05,7.5599847855e-06,1.7638435243800002e-05,1.5155579382299998e-05,8.232701246999999e-06,8.4182405355e-06,1.0923286209299999e-05,2.01405160341e-05,7.383735465300001e-06,5.52827183922e-05,7.44672207456e-05,4.5076216176e-06,6.3409807368e-06,6.0913895615999996e-06,6.023301272100001e-06,2.7098910352799996e-05,1.43892038115e-05,1.96441633269e-05,1.9557932031000003e-06,4.7323285776e-06,8.1609406632e-06,4.252267125e-06,2.5443069732000003e-06,1.1842056792e-05,4.2777651231e-05,1.9588205052e-05,4.5322769646000005e-06,2.17774750284e-05,2.8159527204e-05,2.39869062126e-05,1.577318022e-06,2.11065999156e-05,4.15960481251698e-06,1.3912263288900002e-05,1.5799708125e-05,1.98371928474e-05,3.12956873424e-05,9.8357773446e-06,2.3451916392e-06,2.41207004817e-05,7.1323861662e-06,7.7055865731e-06,1.6665454107e-05,1.33106416128e-05,2.2693447177199998e-05,1.28883485745e-05,5.3329515560999995e-06,5.4104442479999995e-06,4.0392948618e-06,2.79768591711e-05,2.5580702745e-06,2.9903398857e-06,6.2211890403e-06,1.60235412246e-05};
    //static int MutRateSet = 6; // Sets which genome to select.
    static int MutRateSet = EpidermisConst.MutRateSet; // Take from grid
    static double mutRateAdjustment = EpidermisConst.mutRateAdjustment;
    static int corrMutRate = EpidermisConst.corrMutRate;
    static final double[] ExpectedMuts = SetAdjustedMuts(MutRateSet);

    /*
    Use these different Expected Mutation arrays for different mutation tests.
     */
    // Set 1: Orginal. Mean Mutation Rate = 3.2e10-9
//    static final double[] ExpectedMuts = new double[]{3.5,1.2298960944e-05,1.96790865336e-05,5.17789144098534e-06,4.0831751376e-06,3.4550775719999997e-06,5.31298384596e-05,1.09091587995e-05,1.4262587681399999e-05,5.799936312e-06,3.26067900354e-05,7.5642604596e-06,3.1009560480000003e-06,2.632816548e-06,4.9328483301e-06,1.4758253784000002e-06,1.01311941753e-05,1.31902673427e-05,7.5599847855e-06,1.7638435243800002e-05,1.5155579382299998e-05,8.232701246999999e-06,8.4182405355e-06,1.0923286209299999e-05,2.01405160341e-05,7.383735465300001e-06,5.52827183922e-05,7.44672207456e-05,4.5076216176e-06,6.3409807368e-06,6.0913895615999996e-06,6.023301272100001e-06,2.7098910352799996e-05,1.43892038115e-05,1.96441633269e-05,1.9557932031000003e-06,4.7323285776e-06,8.1609406632e-06,4.252267125e-06,2.5443069732000003e-06,1.1842056792e-05,4.2777651231e-05,1.9588205052e-05,4.5322769646000005e-06,2.17774750284e-05,2.8159527204e-05,2.39869062126e-05,1.577318022e-06,2.11065999156e-05,4.15960481251698e-06,1.3912263288900002e-05,1.5799708125e-05,1.98371928474e-05,3.12956873424e-05,9.8357773446e-06,2.3451916392e-06,2.41207004817e-05,7.1323861662e-06,7.7055865731e-06,1.6665454107e-05,1.33106416128e-05,2.2693447177199998e-05,1.28883485745e-05,5.3329515560999995e-06,5.4104442479999995e-06,4.0392948618e-06,2.79768591711e-05,2.5580702745e-06,2.9903398857e-06,6.2211890403e-06,1.60235412246e-05};
    // Set 2: Cancer. Mean Mutation Rate = 3.2e10-6
//    static final double[] ExpectedMuts = new double[]{3.5,0.0118224,0.01891656,0.004977258,0.00392496,0.0033258,0.05107116,0.01048983,0.01371492,0.0055752,0.0104615,0.00727116,0.0029808,0.0025308,0.00474171,0.00110124,0.00973863,0.01267917,0.00726705,0.01696086,0.01456833,0.0079137,0.00805482,0.01051569,0.00435319,0.00710397,0.05314062,0.07159134,0.00433296,0.00568287,0.00553056,0.00573777,0.02604888,0.01384425,0.0179256,0.00152361,0.00454896,0.00744375,0.0040875,0.0020292,0.01135056,0.0411201,0.0188292,0.0043542,0.02093364,0.0270684,0.02305746,0.01478214,0.0015162,0.02010384,0.003998428,0.01337319,0.0142641,0.01906854,0.03008304,0.00946338,0.00225432,0.02318607,0.00685602,0.00701043,0.01449612,0.01279488,0.0218286,0.01238895,0.00512631,0.0052008,0.00388278,0.02689281,0.00245895,0.00228798,0.00598599,0.01540266};
    // Set 3: Mean Mutation Rate = 3.2e10-10
//    static final double[] ExpectedMuts = new double[]{3.5,1.18E-06,1.89E-06,4.98E-07,3.92E-07,3.33E-07,5.11E-06,1.05E-06,1.37E-06,5.58E-07,1.05E-06,7.27E-07,2.98E-07,2.53E-07,4.74E-07,1.10E-07,9.74E-07,1.27E-06,7.27E-07,1.70E-06,1.46E-06,7.91E-07,8.05E-07,1.05E-06,4.35E-07,7.10E-07,5.31E-06,7.16E-06,4.33E-07,5.68E-07,5.53E-07,5.74E-07,2.60E-06,1.38E-06,1.79E-06,1.52E-07,4.55E-07,7.44E-07,4.09E-07,2.03E-07,1.14E-06,4.11E-06,1.88E-06,4.35E-07,2.09E-06,2.71E-06,2.31E-06,1.48E-06,1.52E-07,2.01E-06,4.00E-07,1.34E-06,1.43E-06,1.91E-06,3.01E-06,9.46E-07,2.25E-07,2.32E-06,6.86E-07,7.01E-07,1.45E-06,1.28E-06,2.18E-06,1.24E-06,5.13E-07,5.20E-07,3.88E-07,2.69E-06,2.46E-07,2.29E-07,5.99E-07,1.54E-06};
    // Set 4: Mean Mutation Rate = 3.2e10-7
//    static final double[] ExpectedMuts = new double[]{3.5,0.00118224,0.001891656,0.000497726,0.000392496,0.00033258,0.005107116,0.001048983,0.001371492,0.00055752,0.00104615,0.000727116,0.00029808,0.00025308,0.000474171,0.000110124,0.000973863,0.001267917,0.000726705,0.001696086,0.001456833,0.00079137,0.000805482,0.001051569,0.000435319,0.000710397,0.005314062,0.007159134,0.000433296,0.000568287,0.000553056,0.000573777,0.002604888,0.001384425,0.00179256,0.000152361,0.000454896,0.000744375,0.00040875,0.00020292,0.001135056,0.00411201,0.00188292,0.00043542,0.002093364,0.00270684,0.002305746,0.001478214,0.00015162,0.002010384,0.000399843,0.001337319,0.00142641,0.001906854,0.003008304,0.000946338,0.000225432,0.002318607,0.000685602,0.000701043,0.001449612,0.001279488,0.00218286,0.001238895,0.000512631,0.00052008,0.000388278,0.002689281,0.000245895,0.000228798,0.000598599,0.001540266};
    // Set 5: Mean Mutation Rate = 3.2e10-5
//    static final double[] ExpectedMuts = new double[]{3.5,0.118224,0.1891656,0.049772582,0.0392496,0.033258,0.5107116,0.1048983,0.1371492,0.055752,0.104615,0.0727116,0.029808,0.025308,0.0474171,0.0110124,0.0973863,0.1267917,0.0726705,0.1696086,0.1456833,0.079137,0.0805482,0.1051569,0.0435319,0.0710397,0.5314062,0.7159134,0.0433296,0.0568287,0.0553056,0.0573777,0.2604888,0.1384425,0.179256,0.0152361,0.0454896,0.0744375,0.040875,0.020292,0.1135056,0.411201,0.188292,0.043542,0.2093364,0.270684,0.2305746,0.1478214,0.015162,0.2010384,0.039984282,0.1337319,0.142641,0.1906854,0.3008304,0.0946338,0.0225432,0.2318607,0.0685602,0.0701043,0.1449612,0.1279488,0.218286,0.1238895,0.0512631,0.052008,0.0388278,0.2689281,0.0245895,0.0228798,0.0598599,0.1540266};
    // Set 6: Mean Mutation Rate = 3.2e10-8
//    static final double[] ExpectedMuts = new double[]{3.5,0.000118224,0.000189166,4.98E-05,3.92E-05,3.33E-05,0.000510712,0.000104898,0.000137149,5.58E-05,0.000104615,7.27E-05,2.98E-05,2.53E-05,4.74E-05,1.10E-05,9.74E-05,0.000126792,7.27E-05,0.000169609,0.000145683,7.91E-05,8.05E-05,0.000105157,4.35E-05,7.10E-05,0.000531406,0.000715913,4.33E-05,5.68E-05,5.53E-05,5.74E-05,0.000260489,0.000138443,0.000179256,1.52E-05,4.55E-05,7.44E-05,4.09E-05,2.03E-05,0.000113506,0.000411201,0.000188292,4.35E-05,0.000209336,0.000270684,0.000230575,0.000147821,1.52E-05,0.000201038,4.00E-05,0.000133732,0.000142641,0.000190685,0.00030083,9.46E-05,2.25E-05,0.000231861,6.86E-05,7.01E-05,0.000144961,0.000127949,0.000218286,0.00012389,5.13E-05,5.20E-05,3.88E-05,0.000268928,2.46E-05,2.29E-05,5.99E-05,0.000154027};
    /*
    End Different Expected Mutation arrays
     */

    private static final double[] BaseMutProb = new double[]{1.0/6,3/6.,1.0/6,1.0/6}; // A>,C>,G>,T>
    private static final long[][][] BaseIndex = ParseBaseIndexes();
    static Poisson[] PoissonDists = BuildPoissons(mutRateAdjustment);
    static Poisson[] PoissonDists_corr = BuildPoissons_corr(corrMutRate);
    private static final String[] Base = new String[]{"A","C","G","T"};
    private static final boolean QuickMut = false;
    private static final double QuickMutRate = 0.1;
    static final Random RN=new Random();
    EpidermisGrid theGrid;
    String PrivateGenome;
    float h;
    float s;
    float v;
    int injSite;
    int tp53Clone;
    /*
    End New Information To Keep Inside the Model!!!!!
     */

    EpidermisCellGenome(float h, float s, float v, String PrivateGenome, EpidermisGrid theGrid, int injSite) {
        this.h = h;
        this.s = s;
        this.v = v;
        this.injSite = injSite;
        this.PrivateGenome = PrivateGenome;
        this.theGrid = theGrid;
        this.tp53Clone = tp53Clone;
    }

    //
    @Override
    public EpidermisCellGenome _RunPossibleMutation() {
        StringBuilder MutsObtained = new StringBuilder();

        if (QuickMut == false) {
            for (int j = 0; j < ExpectedMuts.length; j++) {
                if (j != 0) {
//                    Poisson poisson_dist = new Poisson(ExpectedMuts[j], RNEngine); // Setup the Poisson distributions for each gene.
//                    int mutations = poisson_dist.nextInt(); // Gets how many mutations will occur for each gene
                    int mutations = PoissonDists[j].nextInt();
                    for (int hits = 0; hits < mutations; hits++) {
                        int MutatedBaseKind = Utils.RandomVariable(BaseMutProb, RN);
                        long mutIndex = BaseIndex[j - 1][MutatedBaseKind][RN.nextInt(BaseIndex[j - 1][MutatedBaseKind].length)];
                        String MutOut = "";
                        if (j == ExpectedMuts.length - 1) {
                            MutOut = theGrid.GetTick() + "." + j + "." + Base[MutatedBaseKind] + "." + mutIndex;
                        } else {
                            MutOut = theGrid.GetTick() + "." + j + "." + Base[MutatedBaseKind] + "." + mutIndex + ",";
                        }
                        MutsObtained.append(MutOut);
                    }
                }
                //            else {
                //                if(EpidermisConst.GuiOn == true) {
                //                    Poisson poisson_dist = new Poisson(ExpectedMuts[j], RNEngine); // Setup the Poisson distributions for each gene.
                //                    int mutations = poisson_dist.nextInt(); // Gets how many mutations will occur for the Genome
                //                    for (int hits = 0; hits < mutations; hits++) {
                //                        long mutIndex = RN.nextLong();
                //                        String MutOut = "";
                //                        if(j==ExpectedMuts.length-1){
                //                            MutOut = j + "." + ".N." + "." + mutIndex;
                //                        } else {
                //                            MutOut = j + "." + ".N." + "." + mutIndex + ",";
                //                        }
                //                        MutsObtained.append(MutOut);
                //                    }
                //                }
                //            }}


            }
            String PrivGenome = MutsObtained.toString();

            if (PrivGenome.length() > 0) {
                if(PrivGenome.contains(".68.")){
                    EpidermisConst.tp53CloneTracker +=1; // if tp53 mutation, keep track of clone number. will be subtracted if not in basal layer
//
//                    System.out.print("Tp53 clone: " + EpidermisConst.tp53Clone + " ");
                }

                if (h == 0f && s == 0f && v == 1f) {
                    return new EpidermisCellGenome(RN.nextFloat(), 1f , 0.75f, PrivGenome, theGrid,-1);
                } else {
                    return new EpidermisCellGenome(h, RN.nextFloat()*0.3f+0.6f, RN.nextFloat()*0.55f+0.3f, PrivGenome, theGrid,-1);
                }
            } else {
                return null; // If No Mutation Occurs
            }
        }else {
            if(RN.nextDouble()<QuickMutRate) {
                String EmptyGenome = "";
                if (h == 0f && s == 0f && v == 1f) {
                    return new EpidermisCellGenome(RN.nextFloat(), 1f, 0.75f, EmptyGenome, theGrid,-1);
                } else {
                    return new EpidermisCellGenome(h, RN.nextFloat() * 0.3f + 0.6f, RN.nextFloat() * 0.55f + 0.3f, EmptyGenome, theGrid,-1);
                }
            } else {
                return null;
            }
        }
    }

    @Override
    public EpidermisCellGenome _RunPossibleMutation_corrected() {
        StringBuilder MutsObtained = new StringBuilder();
        if (QuickMut == false) {
            for (int j = 0; j < ExpectedMuts.length; j++) {
                if (j != 0) {
//                    Poisson poisson_dist = new Poisson(ExpectedMuts[j], RNEngine); // Setup the Poisson distributions for each gene.
//                    int mutations = poisson_dist.nextInt(); // Gets how many mutations will occur for each gene
                    int mutations = PoissonDists_corr[j].nextInt();
                    for (int hits = 0; hits < mutations; hits++) {
                        int MutatedBaseKind = Utils.RandomVariable(BaseMutProb, RN);
                        long mutIndex = BaseIndex[j - 1][MutatedBaseKind][RN.nextInt(BaseIndex[j - 1][MutatedBaseKind].length)];
                        String MutOut = "";
                        if (j == ExpectedMuts.length - 1) {
                            MutOut = theGrid.GetTick() + "." + j + "." + Base[MutatedBaseKind] + "." + mutIndex;
                        } else {
                            MutOut = theGrid.GetTick() + "." + j + "." + Base[MutatedBaseKind] + "." + mutIndex + ",";
                        }
                        MutsObtained.append(MutOut);
                    }
                }
            }
            String PrivGenome = MutsObtained.toString();

            if (PrivGenome.length() > 0) {
                if(PrivGenome.contains(".68.")){
                    EpidermisConst.tp53CloneTracker +=1; // if tp53 mutation, keep track of clone number. will be subtracted if not in basal layer
//
//                    System.out.print("Tp53 clone: " + EpidermisConst.tp53Clone + " ");
                }

                if (h == 0f && s == 0f && v == 1f) {
                    return new EpidermisCellGenome(RN.nextFloat(), 1f , 0.75f, PrivGenome, theGrid,-1);
                } else {
                    return new EpidermisCellGenome(h, RN.nextFloat()*0.3f+0.6f, RN.nextFloat()*0.55f+0.3f, PrivGenome, theGrid,-1);
                }
            } else {
                return null; // If No Mutation Occurs
            }
        }else {
            if(RN.nextDouble()<QuickMutRate) {
                String EmptyGenome = "";
                if (h == 0f && s == 0f && v == 1f) {
                    return new EpidermisCellGenome(RN.nextFloat(), 1f, 0.75f, EmptyGenome, theGrid,-1);
                } else {
                    return new EpidermisCellGenome(h, RN.nextFloat() * 0.3f + 0.6f, RN.nextFloat() * 0.55f + 0.3f, EmptyGenome, theGrid,-1);
                }
            } else {
                return null;
            }
        }
    }
    // Function to create a new genome with a specific color profile for corrected cells in the grid 28Aug23HLC
    @Override
    public EpidermisCellGenome _RunPossibleCorrection(int siteID){
        String CorrectedGenome = "FAcorrection,";

        if (h == 0f && s == 0f && v == 1f) {
            return new EpidermisCellGenome(RN.nextFloat(), 1f, 0.75f, CorrectedGenome, theGrid, siteID);
        } else {
            return new EpidermisCellGenome(h, RN.nextFloat() * 0.3f + 0.6f, RN.nextFloat() * 0.55f + 0.3f, CorrectedGenome, theGrid, siteID);
        }
        //return new EpidermisCellGenome(.4722f, 0.6f,.67f, EmptyGenome,theGrid); //green corrected color
        //return new EpidermisCellGenome(.861f, 0.6f,.67f, EmptyGenome,theGrid);
    }

    @Override
    public String GenomeInfoStr() {
        return PrivateGenome;
    }

    // Parses Base Mutation Information
    private static long[][][] ParseBaseIndexes(){
        FileIO reader = new FileIO(BaseIndexFile, "r");
        ArrayList<long[]> data = new ArrayList<> (reader.ReadLongDelimit(","));
        long[][][] BaseIndexes = new long[data.size()/4][4][];
        for (int i = 0; i < data.size(); i+=4) {
            BaseIndexes[i/4][0] = data.get(i);
            BaseIndexes[i/4][1] = data.get(i+1);
            BaseIndexes[i/4][2] = data.get(i+2);
            BaseIndexes[i/4][3] = data.get(i+3);
        }
        return BaseIndexes;
    }

    static Poisson[] BuildPoissons(double mutRateModify){
        Poisson[] PoissonDists = new Poisson[ExpectedMuts.length];
        for (int i = 0; i < ExpectedMuts.length; i++) {
            Poisson poisson_dist = new Poisson(ExpectedMuts[i]/mutRateAdjustment, RNEngine);
            PoissonDists[i] = poisson_dist;
        }
        return PoissonDists;
    }
    static Poisson[] BuildPoissons_corr(int mutRateModify){
        Poisson[] PoissonDists = new Poisson[ExpectedMuts.length];
        for (int i = 0; i < ExpectedMuts.length; i++) {
            Poisson poisson_dist = new Poisson(ExpectedMuts[i]/corrMutRate, RNEngine);
            PoissonDists[i] = poisson_dist;
        }
        return PoissonDists;
    }

    private static double[] SetAdjustedMuts(int Selected){
        System.out.println("Mutation Rate Selected: "+ Selected);
        double[] ExpectedMuts;
        if (Selected==0){
            // Set 1: Original. Mean Mutation Rate = 3.2e10-9
            ExpectedMuts = new double[]{3.5,1.2298960944e-05,1.96790865336e-05,5.17789144098534e-06,4.0831751376e-06,3.4550775719999997e-06,5.31298384596e-05,1.09091587995e-05,1.4262587681399999e-05,5.799936312e-06,3.26067900354e-05,7.5642604596e-06,3.1009560480000003e-06,2.632816548e-06,4.9328483301e-06,1.4758253784000002e-06,1.01311941753e-05,1.31902673427e-05,7.5599847855e-06,1.7638435243800002e-05,1.5155579382299998e-05,8.232701246999999e-06,8.4182405355e-06,1.0923286209299999e-05,2.01405160341e-05,7.383735465300001e-06,5.52827183922e-05,7.44672207456e-05,4.5076216176e-06,6.3409807368e-06,6.0913895615999996e-06,6.023301272100001e-06,2.7098910352799996e-05,1.43892038115e-05,1.96441633269e-05,1.9557932031000003e-06,4.7323285776e-06,8.1609406632e-06,4.252267125e-06,2.5443069732000003e-06,1.1842056792e-05,4.2777651231e-05,1.9588205052e-05,4.5322769646000005e-06,2.17774750284e-05,2.8159527204e-05,2.39869062126e-05,1.577318022e-06,2.11065999156e-05,4.15960481251698e-06,1.3912263288900002e-05,1.5799708125e-05,1.98371928474e-05,3.12956873424e-05,9.8357773446e-06,2.3451916392e-06,2.41207004817e-05,7.1323861662e-06,7.7055865731e-06,1.6665454107e-05,1.33106416128e-05,2.2693447177199998e-05,1.28883485745e-05,5.3329515560999995e-06,5.4104442479999995e-06,4.0392948618e-06,2.79768591711e-05,2.5580702745e-06,2.9903398857e-06,6.2211890403e-06,1.60235412246e-05};
        } else if (Selected==1){
            // Set 2: Cancer. Mean Mutation Rate = 3.2e10-6
            ExpectedMuts = new double[]{3.5,0.0118224,0.01891656,0.004977258,0.00392496,0.0033258,0.05107116,0.01048983,0.01371492,0.0055752,0.0104615,0.00727116,0.0029808,0.0025308,0.00474171,0.00110124,0.00973863,0.01267917,0.00726705,0.01696086,0.01456833,0.0079137,0.00805482,0.01051569,0.00435319,0.00710397,0.05314062,0.07159134,0.00433296,0.00568287,0.00553056,0.00573777,0.02604888,0.01384425,0.0179256,0.00152361,0.00454896,0.00744375,0.0040875,0.0020292,0.01135056,0.0411201,0.0188292,0.0043542,0.02093364,0.0270684,0.02305746,0.01478214,0.0015162,0.02010384,0.003998428,0.01337319,0.0142641,0.01906854,0.03008304,0.00946338,0.00225432,0.02318607,0.00685602,0.00701043,0.01449612,0.01279488,0.0218286,0.01238895,0.00512631,0.0052008,0.00388278,0.02689281,0.00245895,0.00228798,0.00598599,0.01540266};
        } else if (Selected==2){
            // Set 3: Mean Mutation Rate = 3.2e10-10
            ExpectedMuts = new double[]{3.5,1.18E-06,1.89E-06,4.98E-07,3.92E-07,3.33E-07,5.11E-06,1.05E-06,1.37E-06,5.58E-07,1.05E-06,7.27E-07,2.98E-07,2.53E-07,4.74E-07,1.10E-07,9.74E-07,1.27E-06,7.27E-07,1.70E-06,1.46E-06,7.91E-07,8.05E-07,1.05E-06,4.35E-07,7.10E-07,5.31E-06,7.16E-06,4.33E-07,5.68E-07,5.53E-07,5.74E-07,2.60E-06,1.38E-06,1.79E-06,1.52E-07,4.55E-07,7.44E-07,4.09E-07,2.03E-07,1.14E-06,4.11E-06,1.88E-06,4.35E-07,2.09E-06,2.71E-06,2.31E-06,1.52E-07,2.01E-06,4.00E-07,1.34E-06,1.43E-06,1.91E-06,3.01E-06,9.46E-07,2.25E-07,2.32E-06,6.86E-07,7.01E-07,1.45E-06,1.28E-06,2.18E-06,1.24E-06,5.13E-07,5.20E-07,3.88E-07,2.69E-06,2.46E-07,2.29E-07,5.99E-07,1.54E-06};
        } else if (Selected==3){
            // Set 4: Mean Mutation Rate = 3.2e10-7
            ExpectedMuts = new double[]{3.5,0.00118224,0.001891656,0.000497726,0.000392496,0.00033258,0.005107116,0.001048983,0.001371492,0.00055752,0.00104615,0.000727116,0.00029808,0.00025308,0.000474171,0.000110124,0.000973863,0.001267917,0.000726705,0.001696086,0.001456833,0.00079137,0.000805482,0.001051569,0.000435319,0.000710397,0.005314062,0.007159134,0.000433296,0.000568287,0.000553056,0.000573777,0.002604888,0.001384425,0.00179256,0.000152361,0.000454896,0.000744375,0.00040875,0.00020292,0.001135056,0.00411201,0.00188292,0.00043542,0.002093364,0.00270684,0.002305746,0.00015162,0.002010384,0.000399843,0.001337319,0.00142641,0.001906854,0.003008304,0.000946338,0.000225432,0.002318607,0.000685602,0.000701043,0.001449612,0.001279488,0.00218286,0.001238895,0.000512631,0.00052008,0.000388278,0.002689281,0.000245895,0.000228798,0.000598599,0.001540266};
        } else if (Selected==4){
            // Set 5: Mean Mutation Rate = 3.2e10-5
            ExpectedMuts = new double[]{3.5,0.118224,0.1891656,0.049772582,0.0392496,0.033258,0.5107116,0.1048983,0.1371492,0.055752,0.104615,0.0727116,0.029808,0.025308,0.0474171,0.0110124,0.0973863,0.1267917,0.0726705,0.1696086,0.1456833,0.079137,0.0805482,0.1051569,0.0435319,0.0710397,0.5314062,0.7159134,0.0433296,0.0568287,0.0553056,0.0573777,0.2604888,0.1384425,0.179256,0.0152361,0.0454896,0.0744375,0.040875,0.020292,0.1135056,0.411201,0.188292,0.043542,0.2093364,0.270684,0.2305746,0.015162,0.2010384,0.039984282,0.1337319,0.142641,0.1906854,0.3008304,0.0946338,0.0225432,0.2318607,0.0685602,0.0701043,0.1449612,0.1279488,0.218286,0.1238895,0.0512631,0.052008,0.0388278,0.2689281,0.0245895,0.0228798,0.0598599,0.1540266};
        } else if (Selected==5){
            // Set 6: Mean Mutation Rate = 3.2e10-8
            ExpectedMuts = new double[]{3.5,0.000118224,0.000189166,4.98E-05,3.92E-05,3.33E-05,0.000510712,0.000104898,0.000137149,5.58E-05,0.000104615,7.27E-05,2.98E-05,2.53E-05,4.74E-05,1.10E-05,9.74E-05,0.000126792,7.27E-05,0.000169609,0.000145683,7.91E-05,8.05E-05,0.000105157,4.35E-05,7.10E-05,0.000531406,0.000715913,4.33E-05,5.68E-05,5.53E-05,5.74E-05,0.000260489,0.000138443,0.000179256,1.52E-05,4.55E-05,7.44E-05,4.09E-05,2.03E-05,0.000113506,0.000411201,0.000188292,4.35E-05,0.000209336,0.000270684,0.000230575,1.52E-05,0.000201038,4.00E-05,0.000133732,0.000142641,0.000190685,0.00030083,9.46E-05,2.25E-05,0.000231861,6.86E-05,7.01E-05,0.000144961,0.000127949,0.000218286,0.00012389,5.13E-05,5.20E-05,3.88E-05,0.000268928,2.46E-05,2.29E-05,5.99E-05,0.000154027};
        }
        // Added 25Aug23HLC no mutations occurring
        else if(Selected==6){
            // Set 7: Mean Mutation Rate = 0
            ExpectedMuts = new double[]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        }
        else if(Selected==7) {
            // Set 7: Mean mutation rate = 3.2e10-9, with updated tp53 based on hg38 gene length
            ExpectedMuts = new double[]{3.5,1.2298960944e-05,1.96790865336e-05,5.17789144098534e-06,4.0831751376e-06,3.4550775719999997e-06,5.31298384596e-05,1.09091587995e-05,1.4262587681399999e-05,5.799936312e-06,3.26067900354e-05,7.5642604596e-06,3.1009560480000003e-06,2.632816548e-06,4.9328483301e-06,1.4758253784000002e-06,1.01311941753e-05,1.31902673427e-05,7.5599847855e-06,1.7638435243800002e-05,1.5155579382299998e-05,8.232701246999999e-06,8.4182405355e-06,1.0923286209299999e-05,2.01405160341e-05,7.383735465300001e-06,5.52827183922e-05,7.44672207456e-05,4.5076216176e-06,6.3409807368e-06,6.0913895615999996e-06,6.023301272100001e-06,2.7098910352799996e-05,1.43892038115e-05,1.96441633269e-05,1.9557932031000003e-06,4.7323285776e-06,8.1609406632e-06,4.252267125e-06,2.5443069732000003e-06,1.1842056792e-05,4.2777651231e-05,1.9588205052e-05,4.5322769646000005e-06,2.17774750284e-05,2.8159527204e-05,2.39869062126e-05,1.577318022e-06,2.11065999156e-05,4.15960481251698e-06,1.3912263288900002e-05,1.5799708125e-05,1.98371928474e-05,3.12956873424e-05,9.8357773446e-06,2.3451916392e-06,2.41207004817e-05,7.1323861662e-06,7.7055865731e-06,1.6665454107e-05,1.33106416128e-05,2.2693447177199998e-05,1.28883485745e-05,5.3329515560999995e-06,5.4104442479999995e-06,4.0392948618e-06,2.79768591711e-05,2.5580702745e-06,3.7824e-06,6.2211890403e-06,1.60235412246e-05};
        }
        else {
            ExpectedMuts = new double[]{0.0};
            System.out.println("You didn't select a mutation rate. This will throw an error.");
        }
        return ExpectedMuts;
    }

    // Parses Base Mutation Function Information
//    public String long[][] ParseMutationInfo(){
//        FileIO reader = new FileIO(BaseIndexFile, "h");
//        ArrayList<long[]> data = new ArrayList<>(reader.ReadBinString(",")));
//    }
}
