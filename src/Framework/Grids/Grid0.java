package Framework.Grids;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Extend the Grid2unstackable class if you want a 2D lattice with at most one agent per grid square
 * @param <T> the AgentSQ2unstackable extending agent class that will inhabit the grid
 */
public class Grid0<T extends Agent0> extends GridBase implements Iterable<T>{
    AgentList<T> agents;
    T[] grid;

    /**
     * @param agentClass pass T.class, used to instantiate agent instances within the grid as needed
     */
    public Grid0(Class<T> agentClass){
        agents=new AgentList<T>(agentClass,this);
    }

    /**
     * returns an uninitialized agent at the specified coordinates
     */
    public T NewAgent(){
        T newAgent=agents.GetNewAgent();
        newAgent.alive=true;
        return newAgent;
    }

    void RemoveAgent(T agent){
        //internal function, removes agent from world
        agents.RemoveAgent(agent);
    }

    /**
     * shuffles the agent list to randomize iteration
     * do not call this while in the middle of iteration
     * @param rn the Random number generator to be used
     */
    public void ShuffleAgents(Random rn){
        agents.ShuffleAgents(rn);
    }

    /**
     * cleans the list of agents, removing dead ones, may improve the efficiency of the agent iteration if many agents have died
     * do not call this while in the middle of iteration
     */
    public void CleanAgents(){
        agents.CleanAgents();
    }

    /**
     * calls CleanAgents, then SuffleAgents, then IncTick. useful to call at the end of a tick
     * do not call this while in the middle of iteration
     * @param rn the Random number generator to be used
     */
    public void CleanShuffInc(Random rn){
        agents.CleanAgents();
        agents.ShuffleAgents(rn);
        IncTick();
    }

    public void CleanInc(){
        agents.CleanAgents();
        IncTick();
    }
    /**
     * returns an umodifiable copy of the complete agentlist, including dead and just born agents
     */
    public ArrayList<T> AllAgents(){return (ArrayList<T>)this.agents.GetAllAgents();}

    /**
     * calls dispose on all agents in the grid
     */
    public void Reset(){
        List<T> AllAgents=this.agents.GetAllAgents();
        AllAgents.stream().filter(curr -> curr.alive).forEach(Agent0::Dispose);
        tick=0;
    }

    /**
     * returns the number of agents that are alive in the grid
     */
    public int GetPop(){
        //gets population
        return agents.pop;
    }

    @Override
    public Iterator<T> iterator(){
        return agents.iterator();
    }

}
