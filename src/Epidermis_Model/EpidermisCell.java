package Epidermis_Model;

import Framework.Grids.AgentSQ3unstackable;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;
import com.sun.org.apache.bcel.internal.generic.NEW;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import static Epidermis_Model.EpidermisCellGenome.ExpectedMuts;
import static Epidermis_Model.EpidermisCellGenome.GeneLengths;
import static Epidermis_Model.EpidermisCellGenome.RN;
import static Epidermis_Model.EpidermisConst.*;

/**
 * Created by schencro on 3/31/17.
 */


class EpidermisCell extends AgentSQ3unstackable<EpidermisGrid> {
    /**
     * parameters that may be changed for cell behavior
     **/
    double prolif_scale_factor = 0.02870462; //Correction for appropriate proliferation rate (Default = 0.15-0.2 with KERATINO_APOPTOSIS_EGF=0.01)
    double KERATINO_EGF_CONSPUMPTION = -0.006954863; //consumption rate by keratinocytes
    double KERATINO_APOPTOSIS_EGF = 0.005939094; //level at which apoptosis occurs by chance (above this and no apoptosis)
    double DEATH_PROB = 0.0038163034; //Overall Death Probability
    double MOVEPROBABILITY = 0.81996739; //RN float has to be greater than this to move...
    double DIVISIONLOCPROB = 0.2518617; // Probability of dividing up vs side to side
    double MUTATION_SCALE_FACTOR = 1; // Adjusts the mutation rates while keeping relative mutation rates the same.
    static int[] dipshit = new int[5];
    static int[] dipshitDiv = new int[5];
    int myType; //cell type
    int Action; //cells action
    static public RandomEngine RNEngine = new DRand();
    int OldLocation;
    int NewLocation;

    /**
     * Parameters for cell specific tracking and genome information
     **/
    // Clonal dynamic tracking
    EpidermisCellGenome myGenome; // Creating genome class within each cell

    public void init(int cellType, EpidermisCellGenome myGenome) { //This initilizes an agent with whatever is inside of this function...
        this.myType = cellType;
        this.Action = STATIONARY;
        // Storing Genome Reference to Parent and Itself if mutation happened
        this.myGenome = myGenome;
    }

    // Gets where a cell is dividing if it's a basal cell and is proliferating
    public int ProlifLoc(){
        double divideWhere = G().RN.nextDouble();
        double OtherOptionProb=(1-DIVISIONLOCPROB)/4.0;
        if(divideWhere<=DIVISIONLOCPROB){
            dipshitDiv[4] ++;
            return 4; // Dividing up
        } else if(divideWhere <= (DIVISIONLOCPROB+OtherOptionProb)){
            dipshitDiv[0] ++;
            return 0; // Dividing right
        } else if (divideWhere <= (DIVISIONLOCPROB+OtherOptionProb*2)){
            dipshitDiv[1] ++;
            return 1; // Dividing Left
        } else if (divideWhere <= (DIVISIONLOCPROB+OtherOptionProb*3)) {
            dipshitDiv[3] ++;
            return 3; // Dividing front
        } else {
            dipshitDiv[2] ++;
            return 2; // Dividing back
        }
    }


    //Checking if a cell is going to proliferate...
    public boolean CheckProliferate() {
        int x = Xsq();
        int y = Ysq();
        int z = Zsq();
        int i = Isq();
        int iDivLoc;

        // TODO Implement change to acquire inBounds for the cell.

        // If EGF is low then next double is likely to be higher...Results in no proliferation
        if (myType == KERATINOCYTE && G().RN.nextDouble() > G().EGF.GetCurr(x, y, z) * prolif_scale_factor) {
            return false;
        }

        GetCoords(G().divHood, G().DIV);
        iDivLoc = ProlifLoc(); // Where the new cell is going to be (which index) if basal cell

        boolean Pushed = CellPush(iDivLoc);

        if(Pushed!=false && y==0 && (iDivLoc==0 || iDivLoc==1 || iDivLoc==3 || iDivLoc==2)){
            G().Turnover.RecordLossBasal(); // Record Cell Loss from Pushing
        }
        if(Pushed==false){
            return false; // Only false if NOTCH mutation
        }

        EpidermisCell newCell = G().NewAgentI(G().inBounds[iDivLoc]);

        newCell.init(myType, myGenome.NewChild().PossiblyMutate()); // initializes a new skin cell, pass the cellID for a new value each time.
        myGenome = myGenome.PossiblyMutate(); // Check if this duaghter cell, i.e. the progenitor gets mutations during this proliferation step.

        if(newCell.Ysq()==0){
            G().Turnover.RecordDivideBasal();
            G().Turnover.RecordDivideTissue();
        } else {
            G().Turnover.RecordDivideTissue();
        }

        G().divisions[G().GetTick()*ySize+Ysq()]++;
        G().divs++;

        return true;
    }

    public boolean CellPush(int iDivLoc){
        int i = G().inBounds[iDivLoc];
        EpidermisCell c=G().GetAgent(i);
        if(c!=null){
            if(EpidermisConst.NOTCH1FitnessChanges) {
                String cellToMoveGenome = c.myGenome.GenomeInfoStr();
                if (cellToMoveGenome.contains(".44.")) {
                    if (EpidermisConst.NOTCHBlockProbability > RN.nextDouble()) {
                        return false;
                    }
                }
            }
            int x = G().ItoX(i);
            int y = G().ItoY(i);
            int z = G().ItoZ(i);
            //look up for empty square
            int colTop=y;
//            EpidermisCell c=G().ItoAgent(i);
            while(c!=null){
                colTop++;
                c=G().GetAgent(x,colTop,z);
            }
            //move column of cells up
            for(;colTop>y;colTop--){
                c=(G().GetAgent(x,colTop-1,z));
                c.MoveSQ(x,colTop,z);
            }
            if(c.Ysq()>= G().yDim-2){c.itDead();}
            return true;

        } else{
            return false;
        }
    }

    // Sets the coordinates for a cell that is moving.
    public int GetCoords(int[] hood, int ACTION) {
        int iMoveCoord=-1;  //when it's time to move, it is the index of coordinate that is picked from Coords array above. -1 == Not Moving
        int finalCount=0;
        int inBoundsCount = G().SQstoLocalIs(hood, G().inBounds, Xsq(),Ysq(), Zsq(), true, false, true); // Gets all inbound indices
        if(ACTION==G().MOVE) {
            for (int i = 0; i < inBoundsCount; i++) {
                if (G().GetAgent(G().inBounds[i]) == null) {
                    G().inBounds[finalCount] = G().inBounds[i];
                    finalCount++;
                }
            }
            if(finalCount>0&&myType==KERATINOCYTE) {
                iMoveCoord=G().RN.nextInt(finalCount);
            }
            // Return statement is only used for movement. Not division.
            return iMoveCoord;
        } else if(ACTION==G().DIV){
            return -1; // Not used for division but setting G().inBounds.
        } else{
            throw new RuntimeException("EpidermisCell.GetCoords() is not working.");
        }
    }

    public void itDead(){
        myGenome.DisposeClone(); // Decrements Population
        Dispose();

        G().MeanDeath[Isq()] += 1;

        if(Ysq()==0){
            G().Turnover.RecordLossBasal();
        }
        G().Turnover.RecordLossTissue();
    }

    public void CellStep(){

        if(!EpidermisGrid.corrected){
            Correct(); // Add step to correct specific cells 28Aug23HLC
        }

        int x=Xsq();int y=Ysq();int z=Zsq(); // Get discrete x and y coordinates
        Action = STATIONARY;
        if (y>=G().AIR_HEIGHT){
            itDead();
            return;
        }
        if (G().EGF.GetCurr(x, y, z) < KERATINO_APOPTOSIS_EGF && G().RN.nextDouble() < (Math.pow(1.0 - G().EGF.GetCurr(x, y, z) / KERATINO_APOPTOSIS_EGF, 5))) {
            //DEATH FROM LACK OF NUTRIENTS KERATINOCYTE
            itDead();
            return;
        }

        // If doing P53 Fitness
        if(EpidermisConst.PFiftyThree) {
            String thisGenome = myGenome.GenomeInfoStr();
            // If P53 Mutation present standard death function
            if(thisGenome.contains(".68.")) { //Implement Check for mutation in cell
                return;
            // If P53 Mutation NOT present
            } else {
                if(RN.nextDouble() < (DEATH_PROB)) {
                    //Random Fucked
                    itDead();
                    return;
                }
            }
        } else {
            if(RN.nextDouble() < DEATH_PROB){
                //Random Fucked
                itDead();
                return;
            }
        }

        if (G().RN.nextFloat() >= MOVEPROBABILITY) {
            int iMoveCoord = GetCoords(G().moveHood, G().MOVE); // -1 if not moving
            if (iMoveCoord != -1) {
                dipshit[ DirectionTracker(G().inBounds[iMoveCoord]) ] ++;
                MoveI(G().inBounds[iMoveCoord]); // We are moving
                Action = MOVING;
                if (Ysq() != 0 && y == 0) {
                    if (Ysq() > y) {
                        throw new RuntimeException("Cell is Moving Up.");
                    }
                }
            }
        }

        boolean divided = CheckProliferate();
        if(divided){
            Action = DIVIDE;
        }

    }


    public int DirectionTracker(int NextMoveIndex){
        int x=G().ItoX(NextMoveIndex);
        int z=G().ItoZ(NextMoveIndex);
        int dx = Xsq()-x;
        int dz = Zsq()-z;
        if (dx == 1) { return 0; }
        if (dx == -1) { return 1; }
        if (dz == 1) { return 2; }
        if (dz == -1) { return 3; }
        return 4;
    }


    public void Correct() { //method to use the list of coordinate points generated in EpidermisGrid and perform the correction 29aug23HLC

        for(int i=0; i<EpidermisGrid.CorrectionPoints.size(); i++){
            int[] coords = EpidermisGrid.CorrectionPoints.get(i);

            // find the cells associated with the coordinates in the list
            if(coords[0] != this.Xsq() || coords[1] != this.Ysq() || coords[2] != this.Zsq()){
                continue;
            }
            else {
                itDead(); // kill and remove the current cell
                EpidermisCell newCell = G().NewAgentPT(coords[0],coords[1],coords[2]); // adds a new cell
                newCell.init(myType, myGenome.NewChild().PerformCorrection()); // initializes a new skin cell
                //myGenome = myGenome.PerformCorrection(); // Appears to be redunant from line above: cells still become corrected
                                                        // keeping line in case this may be for genome tracking purposes
            }
        }
    }

//    // Builds my genome information for data analysis
//    String ToString(){
//        String cellInfo="{["+createStrID()+"];[";
//        for(int iGene=0;iGene<myGenome.genomelength;iGene++){
//            cellInfo+="[";
//            for(int iMut=0;iMut<myGenome.mut_pos[iGene].size();iMut++){
//                cellInfo+=myGenome.mut_pos_getter(iGene,iMut)+",";
//            }
//            cellInfo+="],";
//        }
//        cellInfo+="]}";
//        return cellInfo;
//    }

}
