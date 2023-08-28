package Framework.Grids;

abstract class GridBase {
    int tick;
    public GridBase(){
    }

    /**
     * gets the current grid tick
     */
    public int GetTick(){
        return tick;
    }

    /**
     * increments the current grid tick
     */
    public void IncTick(){
        tick+=1;
    }
}

