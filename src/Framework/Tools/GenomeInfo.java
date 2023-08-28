package Framework.Tools;

/**
 * should be declared myType extends GenomeInfo <myType>
 */
//TODO: MAKE THE GENOMEINFO AND GENOMETRACKER LIKE CHANDLER'S, with a live_pops list
public class GenomeInfo <T extends GenomeInfo> {
    int id;
    long pop;
    T next;
    T prev;
    GenomeTracker<T> myTracker;
    String MutantInfo;

    /**
     * gets the current number of clones that share this genome
     */
    public long GetPop(){
        return pop;
    }

    /**
     * ignore
     */
    void _Init(GenomeTracker myTracker,int id,T next,T prev){
        this.myTracker=myTracker;
        this.id=id;
        this.next=next;
        this.prev=prev;
        this.pop =0;
    }

    /**
     * removes clone from GenomeInfo population
     */
    public void DecPop(){
        myTracker.SetClonePop((T) this,pop -1);
    }

    /**
     * adds new clone to GenomeInfo population
     */
    public void IncPop(){
        if(pop<=0&&this!=myTracker.progenitors){
            throw new IllegalStateException("reviving extinct clone");
        }
        if(pop<=0){
            myTracker._ReviveProgenitors();
        }
        pop++;
    }

    public void SetPop(long newPop){
        if(pop<=0&&this!=myTracker.progenitors){
            throw new IllegalStateException("reviving extinct clone");
        }
        if(pop<=0&&newPop>0){
            myTracker._ReviveProgenitors();
        }
        myTracker.SetClonePop((T) this,newPop);
    }

    public String GetMutantInfo(){
        if(myTracker.allGenomeInfos==null){
            throw new IllegalStateException("Mutant to String function was not defined in GenomeTracker constructor");
        }
        return myTracker.MutantToString.MutantToString(this);
    }

    public T NewMutantGenome(){
        T newGenome=myTracker.NewMutant((T)this);
        pop--;
        newGenome.pop++;
        return newGenome;
    }

    public String GetMutationInfo(){
        return MutantInfo;
    }
    public String FullLineageInfoStr(String delim){
        return myTracker.FullLineageInfoStr(id,delim);
    }

    public GenomeTracker<T> MyTracker(){
        return myTracker;
    }
}
