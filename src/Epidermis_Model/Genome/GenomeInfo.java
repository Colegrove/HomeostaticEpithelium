package Epidermis_Model.Genome;

import Epidermis_Model.EpidermisCellGenome;
import Epidermis_Model.Epidermis_Main;

/**
 * should be declared myType extends GenomeInfo <myType>
 */
public abstract class GenomeInfo <T extends GenomeInfo> {
    public int id;
    int popSize;
    T next;
    T prev;
    GenomeTracker myTracker;


    /**
     * gets the current number of clones that share this genome
     */
    public int GetClonePop(){
        return popSize;
    }

    /**
     * ignore
     */
    void _Init(GenomeTracker myTracker, int id, T next, T prev){
        this.myTracker=myTracker;
        this.id=id;
        this.next=next;
        this.prev=prev;
        this.popSize=0;
    }

    /**
     * removes clone from GenomeInfo population
     */
    public void DisposeClone(){
        myTracker.DisposeClone(this);
    }

    /**
     * adds new clone to GenomeInfo population
     */
    public T NewChild(){
        popSize++;
        return (T)this;
    }

    /**
     * returns a GenomeInfo instance that is either identical to the original or a mutant
     */
    public T PossiblyMutate(){
        T nextGenome= (T) _RunPossibleMutation();
        if(nextGenome==null){
            return (T)this;
        }
        DisposeClone();
        myTracker.AddMutant(this,nextGenome);
        return nextGenome;
    }
    /**
     * returns a GenomeInfo instance that is either identical to the original or a mutant
     */
    public T PossiblyMutate_corrected(){
        T nextGenome= (T) _RunPossibleMutation_corrected();
        if(nextGenome==null){
            return (T)this;
        }
        DisposeClone();
        myTracker.AddMutant(this,nextGenome);
        return nextGenome;
    }

    public T PerformCorrection(int injID){
        T nextGenome = (T) _RunPossibleCorrection(injID);
        DisposeClone();
        myTracker.AddMutant(this,nextGenome);
        return nextGenome;
    }
    public String FullLineageInfoStr(String delim){
        return myTracker.FullLineageInfoStr(id,delim);
    }

    /**
     * a potential mutation event, return null if the genome did not change, otherwise return a new genome with the change inside
     * do not change the calling genome!!!!!!!!!!!!!!
     */
    //public abstract T _RunPossibleMutation();

    //public abstract EpidermisCellGenome _RunPossibleCorrection();

    //
    public abstract EpidermisCellGenome _RunPossibleMutation();

    public abstract EpidermisCellGenome _RunPossibleMutation_corrected();

    // Function to create a new genome with a specific color profile for corrected cells in the grid 28Aug23HLC
    public abstract EpidermisCellGenome _RunPossibleCorrection(int siteID);

    /**
     * returns a string with info about the genome to be stored
     */
    public abstract String GenomeInfoStr();

    public int IDGetter(){ return id; }
}