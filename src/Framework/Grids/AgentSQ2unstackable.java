package Framework.Grids;

import java.util.ArrayList;

/**
 * extend the AgentSQ2unstackable class if you want agents that exist on a 2D discrete lattice
 * without the possibility of stacking multiple agents on the same grid square
 * @param <T> the extended Grid2unstackable class that the agents will live in
 * Created by rafael on 11/18/16.
 */

public class AgentSQ2unstackable<T extends Grid2> extends AgentBaseSpatial <T>{
    int xSq;
    int ySq;

    public void SwapPosition(AgentBaseSpatial other){
        if(!alive||!other.alive){
            throw new RuntimeException("attempting to move dead agent");
        }
        if(other.myGrid!=myGrid){
            throw new IllegalStateException("can't swap positions between agents on different grids!");
        }
        int iNew=other.Isq();
        int iNewOther=Isq();
        other.RemSQ();
        this.RemSQ();
        other.MoveI(iNewOther);
        this.MoveI(iNew);
    }
    void Setup(double i){
        Setup((int)i);
    }
    void Setup(double xSq,double ySq){
        Setup((int)xSq,(int)ySq);
    }
    void Setup(double xSq,double ySq,double zSq){
        throw new IllegalStateException("shouldn't be adding 2D agent to 3D grid");
    }

    @Override
    void Setup(int i) {
        xSq=myGrid.ItoX(i);
        ySq=myGrid.ItoY(i);
        iSq=i;
        AddSQ(i);
    }

    @Override
    void Setup(int x, int y) {
        this.xSq=x;
        this.ySq=y;
        iSq=myGrid.I(xSq,ySq);
        AddSQ(iSq);
    }

    @Override
    void Setup(int x, int y, int z) {
        throw new IllegalStateException("shouldn't be adding 2D agent to 3D grid");
    }

    /**
     * Moves the agent to the square with the specified index
     */
    public void MoveI(int i){
        //moves agent discretely
        if(!this.alive){
            throw new RuntimeException("Attempting to move dead agent!");
        }
        this.xSq=myGrid.ItoX(i);
        this.ySq=myGrid.ItoY(i);
        myGrid.grid[iSq]=null;
        iSq=i;
        AddSQ(i);
    }
    void AddSQ(int i){
        if(myGrid.grid[iSq]!=null){
            throw new RuntimeException("Adding multiple unstackable agents to the same square!");
        }
        myGrid.grid[iSq]=this;
    }
    void RemSQ(){
        myGrid.grid[iSq]=null;
    }

    /**
     * Moves the agent to the square at the specified coordinates
     */
    public void MoveSQ(int x, int y){
        if(!this.alive){
            throw new RuntimeException("Attempting to move dead agent!");
        }
        int iNewPos=myGrid.I(x,y);
        RemSQ();
        AddSQ(iNewPos);
        this.xSq=x;
        this.ySq=y;
        this.iSq=iNewPos;
    }

    /**
     * Gets the xDim coordinate of the square that the agent occupies
     */
    public int Xsq(){
        return xSq;
    }
    /**
     * Gets the yDim coordinate of the square that the agent occupies
     */
    public int Ysq(){
        return ySq;
    }
    /**
     * Gets the xDim coordinate agent
     */
    public double Xpt(){
        return xSq+0.5;
    }
    /**
     * Gets the yDim coordinate agent
     */
    public double Ypt(){
        return ySq+0.5;
    }

    /**
     * Deletes the agent
     */
    public void Dispose(){
        if(!this.alive){
            throw new RuntimeException("Attempting to dispose already dead agent!");
        }
        RemSQ();
        myGrid.agents.RemoveAgent(this);
    }
    public void GetAllOnSquare(ArrayList<AgentBaseSpatial> putHere){
        putHere.add(this);
    }

    /**
     * Gets the index of the square that the agent occupies
     */
    public int Isq(){
        return iSq;
    }
}
