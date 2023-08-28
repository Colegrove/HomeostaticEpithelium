package Framework.Grids;

/**
 * extend the AgentSQ2 class if you want agents that exist on a 2D discrete lattice
 * with the possibility of stacking multiple agents on the same grid square
 * @param <T> the extended Grid2 class that the agents will live in
 * Created by rafael on 11/18/16.
 */
abstract public class Agent0<T extends Grid0> extends AgentBase<T>{


    /**
     * deletes the agent
     */
    public void Dispose(){
        //kills agent
        myGrid.RemoveAgent(this);
    }

    /**
     * Gets the index of the square that the agent occupies
     */
    //addCoords
}
