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

    int divisionDirection; // added to this level for corrected cells to access if needs to change during division. 16OCT23HLC
    int tp53Clone; // keeps track of order that tp53 mutations arise. This allows tp53 clones to gain new mutations

    /**
     * Parameters for cell specific tracking and genome information
     **/
    // Clonal dynamic tracking
    EpidermisCellGenome myGenome; // Creating genome class within each cell

    public void init(int cellType, EpidermisCellGenome myGenome, int tp53Clone) { //This initilizes an agent with whatever is inside of this function...
        this.myType = cellType;
        this.Action = STATIONARY;
        // Storing Genome Reference to Parent and Itself if mutation happened
        this.myGenome = myGenome;
        this.tp53Clone = -1;
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
        //int iDivLoc;

        // TODO Implement change to acquire inBounds for the cell.


        EpidermisCell c = G().GetAgent(i); //03OCT23 HLC
        String cellToGrowGenome = c.myGenome.FullLineageInfoStr(""); // 03OCT23HLC get full genome info
        // If EGF is low then next double is likely to be higher...Results in no proliferation
        // 03OCT23HLC if a corrected cell and using corrected growth advantage, modify proliferation scaling factor


        if (cellToGrowGenome.contains("FAcorrection") && CorrectedGrowthChanges){
            if(cellToGrowGenome.contains(".68") && p53Growth){
                if (myType == KERATINOCYTE && G().RN.nextDouble() > G().EGF.GetCurr(x, y, z) * (prolif_scale_factor * CorrectedGrowthIncrease * tp53GrowthIncrease)) {
                    //System.out.println("FA + tp53");
                    //System.out.print("dividing corrected cell\n");
                    return false;
                }
            }
            else if (myType == KERATINOCYTE && G().RN.nextDouble() > G().EGF.GetCurr(x, y, z) * (prolif_scale_factor * CorrectedGrowthIncrease)) {
                //System.out.print("dividing corrected cell\n");
                return false;
            }
        }
        else if (cellToGrowGenome.contains(".68.") && p53Growth){
            if (myType == KERATINOCYTE && G().RN.nextDouble() > G().EGF.GetCurr(x, y, z) * (prolif_scale_factor * tp53GrowthIncreaseFA)) {
                //System.out.print("dividing corrected cell\n");
                return false;
            }
        }
        else if (myType == KERATINOCYTE && G().RN.nextDouble() > G().EGF.GetCurr(x, y, z) * prolif_scale_factor) {
            //System.out.print("dividing non-corrected cell\n");
            return false;
        }

        GetCoords(G().divHood, G().DIV);
        divisionDirection = ProlifLoc(); // Where the new cell is going to be (which index) if basal cell

        boolean Pushed = CellPush();

        if(Pushed!=false && y==0 && (divisionDirection==0 || divisionDirection==1 || divisionDirection==3 || divisionDirection==2)){
            G().Turnover.RecordLossBasal(); // Record Cell Loss from Pushing
        }
        if(Pushed==false){
            return false; // Only false if NOTCH mutation
        }

        //cell division
        //System.out.print("dividing cell - x position = " + x + " y position = " + y + " time = " + G().GetTick() + "");

        EpidermisCell newCell = G().NewAgentI(G().inBounds[divisionDirection]);

        int currentTp53Clone = tp53Clone;
        int currentInjSite = myGenome.injSite;


        // do daughter cells acquire a mutation
        // based on if cell is corrected or not, perform potential to mutate based on specific mutation rates
        if(newCell.myGenome.injSite == 0){
            newCell.init(myType, myGenome.NewChild().PossiblyMutate_corrected(), tp53Clone); // initializes a new skin cell, pass the cellID for a new value each time.
        }else{
            newCell.init(myType, myGenome.NewChild().PossiblyMutate(), tp53Clone); // initializes a new skin cell, pass the cellID for a new value each time.
        }
        if(myGenome.injSite == 0){
            myGenome = myGenome.PossiblyMutate_corrected(); // Check if this corrected daughter cell, i.e. the progenitor gets mutations during this proliferation step.
        }else{
            myGenome = myGenome.PossiblyMutate(); // Check if this daughter cell, i.e. the progenitor gets mutations during this proliferation step.
        }


        // keep prior injection site info if cells gain new mutation
        newCell.myGenome.injSite = currentInjSite;
        myGenome.injSite = currentInjSite;


        if(newCell.myGenome.FullLineageInfoStr("").contains(".68.")){
                // check to see what time step the mutation occurs at
                // if the mutation occurs at current timestep, then subtract clone count
                boolean currentMutation = false;
                //String[] elements = newCell.myGenome.PrivateGenome.split(",");
                String[] elements = newCell.myGenome.FullLineageInfoStr("").split(",");
                for (String element : elements) {
                    if (element.contains(".68.")) {
                        String[] parts = element.split("\\.");
                        if (parts.length >= 2) {
                            String startingNumber = parts[0];
                            int startingNumberInt = Integer.parseInt(startingNumber);
                            if(startingNumberInt == G().GetTick()){
                                currentMutation = true;
                            }}
                    }
                }

                if(newCell.Ysq() > 0 && currentMutation){
                    tp53CloneTracker -= 1;
                } else if(newCell.Ysq() == 0 && currentMutation) { // basal layer cell
                    newCell.tp53Clone = tp53CloneTracker;
                } else if(newCell.Ysq() == 0){ // new cell but no new mutations
                    newCell.tp53Clone = this.tp53Clone;
                }
        }

        if(myGenome.FullLineageInfoStr("").contains(".68.")){
            boolean currentMutation = false;
            boolean currentNontp53Mutation = false;
            // If the tp53 mutation is not in the basal layer
            // check to see what time step the mutation occurs at
            // if the mutation occurs at current timestep, then subtract clone count
            String[] elements = myGenome.FullLineageInfoStr("").split(",");
            for (String element : elements) {
                if (element.contains(".68.")) {
                    String[] parts = element.split("\\.");
                    if (parts.length >= 2) {
                        String startingNumber = parts[0];
                        int startingNumberInt = Integer.parseInt(startingNumber);
                        if (startingNumberInt == G().GetTick()) {
                            currentMutation = true;
                        }
                    }
                } else {
                    String[] parts = element.split("\\.");
                    if (parts.length >= 2) {
                        String startingNumber = parts[0];
                        int startingNumberInt = Integer.parseInt(startingNumber);
                        if (startingNumberInt == G().GetTick()) {
                            currentNontp53Mutation = true;
                        }
                    }}}
            if(Ysq() > 0 && currentMutation){
                tp53CloneTracker -= 1;
            } else if(Ysq() == 0 && currentMutation){
                tp53Clone = tp53CloneTracker;
            } else if(Ysq() == 0 && currentNontp53Mutation){
                this.tp53Clone = currentTp53Clone;
            }
        }

        if(newCell.Ysq()==0){
            G().Turnover.RecordDivideBasal();
            G().Turnover.RecordDivideTissue();
        } else {
            G().Turnover.RecordDivideTissue();
        }

        G().divisions[G().GetTick()*ySize+Ysq()]++;
        G().divs++;

        // if tracking divisions, add division event to an array
//        if(divRate){
//            EpidermisGrid.divArray[c.Xsq()][c.Zsq()][G().GetTick()] = 1;
//        }
        return true;
    }

    public boolean CellPush(){
        int i = G().inBounds[divisionDirection];
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

            if(CorrectedBlockChanges || p53Blocking) { // if cell is corrected apply the blocking probability 31Aug23HLC
                //String cellToMoveGenome = c.myGenome.GenomeInfoStr();
                String cellToMoveGenome = c.myGenome.FullLineageInfoStr("");
                int blockSwitch = 0;
                if(cellToMoveGenome.contains("FAcorrection") && cellToMoveGenome.contains(".68.")){
                    blockSwitch = 1; // correction + tp53
                }
                else if(cellToMoveGenome.contains(".68.")){
                    blockSwitch = 2; // FA background + tp53
                }
                else if(cellToMoveGenome.contains("FAcorrection")){
                    blockSwitch = 3; // correction only
                }
                switch(blockSwitch){
                    case 1: // correction + tp53
                        if (CorrectedBlockProbability + tp53BlockProb > RN.nextDouble()) {
                            // if cell is blocked via blocking probability, change direction to divide up.
                            if (allowVerticalDivision) {
                                divisionDirection = 4;
                                i = G().inBounds[divisionDirection];
                                c = G().GetAgent(i);
                                if (c == null) {
                                    return false;
                                }
                            } else {
                                return false;
                            }
                        }
                        break;
                    case 2: // FA background + tp53
                        if (FAtp53BlockProb > RN.nextDouble()) {
                            // if cell is blocked via blocking probability, change direction to divide up.
                            if(allowVerticalDivision){
                                divisionDirection = 4;
                                i = G().inBounds[divisionDirection];
                                c=G().GetAgent(i);
                                if(c==null){
                                    return false;
                                }
                            }else{
                                return false;
                            }

                        }
                        break;
                    case 3: // corrected only
                        if (CorrectedBlockProbability > RN.nextDouble()) {
                            // if cell is blocked via blocking probability, change direction to divide up.
                            if(allowVerticalDivision){
                                divisionDirection = 4;
                                i = G().inBounds[divisionDirection];
                                c=G().GetAgent(i);
                                if(c==null){
                                    return false;
                                }
                            }else{
                                return false;
                            }
                        }
                        break;
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
        //tp53Clone = -1;
        G().MeanDeath[Isq()] += 1;

        if(Ysq()==0){
            //System.out.print("\n>>> >>>> >>> >> Basal cell death");
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

            if (y != 0) {
                itDead();
            }
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
                    if(y !=0) {
                        itDead();
                    }
                    return;
                }
            }
        } else {
            if(RN.nextDouble() < DEATH_PROB){
                //Random Apoptosis

                if(y != 0) {
                    itDead();
                }
                return;
            }
        }

        if (G().RN.nextFloat() >= MOVEPROBABILITY) {
            int iMoveCoord = GetCoords(G().moveHood, G().MOVE); // -1 if not moving
            if (iMoveCoord != -1) {
                dipshit[ DirectionTracker(G().inBounds[iMoveCoord]) ] ++;
                MoveI(G().inBounds[iMoveCoord]); // We are moving
                Action = MOVING;
                if(Ysq() == 0 && y !=0){
                    System.out.print("CELL MOVED DOWN FROM POSITION - Y: " + y + "to Y: " + Ysq() + "\n");
                }
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
                newCell.init(myType, myGenome.NewChild().PerformCorrection(coords[3]), tp53Clone); // initializes a new skin cell
                //myGenome = myGenome.PerformCorrection(); // Appears to be redunant from line above: cells still become corrected

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
