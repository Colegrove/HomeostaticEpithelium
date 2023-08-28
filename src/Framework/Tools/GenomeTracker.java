package Framework.Tools;

import Framework.Misc.GetMutationInfo;

import java.lang.reflect.Constructor;
import java.util.*;

/**
 * should be instantiated GenomeTracker<myGenomeInfoType>
 */
public class GenomeTracker <T extends GenomeInfo> implements Iterable<T>{
    ArrayList<Integer> parentIDs;
    ArrayList<String> allGenomeInfos;
    ArrayList<long[]> CloneCounts;
    final ArrayList<T> deadClones;
    final Constructor<?> builder;
    T livingCloneInfos;
    T progenitors;
    int nMutations=0;
    int nActiveClones=0;
    GetMutationInfo MutantToString;

    /**
     * creates a new genome tracker, with a progenitor genome that will be returned whenever GetProgenitor is called
     * @param progenitorGenome the starting genome
     * @param trackParents whether to track the lineages
     */
    public GenomeTracker(Class<T> progenitorGenome,boolean trackParents){
        this.builder=progenitorGenome.getDeclaredConstructors()[0];
        this.builder.setAccessible(true);
        deadClones=new ArrayList<>();
        CloneCounts=new ArrayList<>();
        try {
            progenitors = (T)builder.newInstance();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (Exception e) {
            System.err.println("could not instantiate");
            e.printStackTrace();
        }
        if(trackParents) {
            parentIDs = new ArrayList<>();
            parentIDs.add(-1);
        }
        progenitors._Init(this,0,null,null);
    }
    public GenomeTracker(T progenitorGenome, boolean trackParents, GetMutationInfo<T> MutantToString){
        this.builder=progenitorGenome.getClass().getDeclaredConstructors()[0];
        this.builder.setAccessible(true);
        this.MutantToString=MutantToString;
        deadClones=new ArrayList<>();
        CloneCounts=new ArrayList<>();
        progenitors = progenitorGenome;
        if(trackParents) {
            parentIDs = new ArrayList<>();
            parentIDs.add(-1);
        }
            allGenomeInfos = new ArrayList<>();
            allGenomeInfos.add(MutantToString.MutantToString(progenitors));
        progenitorGenome._Init(this,0,null,null);
    }

    /**
     * gets the parent ids list
     */
    public List<Integer> GetParentIDs(){
        return Collections.unmodifiableList(parentIDs);
    }

    /**
     * gets the info on all genomes
     */
    public List<String> GetGenomeInfos(){
        return Collections.unmodifiableList(allGenomeInfos);
    }

    /**
     * returns all clone counts that have been recorded
     */
    public List<long[]> GetCloneCounts(){
        return Collections.unmodifiableList(CloneCounts);
    }

    /**
     * returns the progenitor genome, and increases the population size
     */
    public T GetProgenitor(){
        return progenitors;
    }
    void _ReviveProgenitors(){
        nActiveClones++;
        if(livingCloneInfos!=null) { livingCloneInfos.prev = progenitors; }
        progenitors.next=livingCloneInfos;
        progenitors.prev=null;
        livingCloneInfos =progenitors;
    }

    /**
     * ignore
     */
    T NewMutant(T parent) {
        nMutations++;
        nActiveClones++;
        T child=null;
        if(deadClones.size()>0){
            child=deadClones.remove(deadClones.size()-1);
        }
        else{
            try {
                child=(T)builder.newInstance();
            } catch (Exception e) {
                System.err.println("could not instantiate");
                e.printStackTrace();
            }
        }
        return child;
    }

    void AddMutant(T parent,T child){
        if(allGenomeInfos!=null) { allGenomeInfos.add(MutantToString.MutantToString(child)); }
        if(parentIDs!=null) { parentIDs.add(parent.id); }
        child._Init(this, GetNumMutations(), livingCloneInfos,null);
        child.pop++;
        if(livingCloneInfos!=null) { livingCloneInfos.prev = child; }
        livingCloneInfos =child;
    }

    /**
     * returns the total number of unique clone populations that have ever existed
     */
    public int GetNumMutations(){
        return nMutations;
    }

    /**
     * returns the total number of unique clone populations that currently exist
     */
    public int GetNumLivingGenomes(){
        return nActiveClones;
    }
    void SetClonePop(T cloneInfo,long newPop){
        cloneInfo.pop =newPop;
        if(cloneInfo.pop <0){
            throw new IllegalStateException("a population size is less than 0");
        }
        if(cloneInfo.pop ==0){
            if(cloneInfo!=progenitors){
                deadClones.add(cloneInfo);
            }
            nActiveClones--;
            //remove livingCloneInfos with 0 population
            if(livingCloneInfos ==cloneInfo){
                livingCloneInfos =(T)cloneInfo.next;
            }
            if(cloneInfo.next!=null){
                cloneInfo.next.prev=cloneInfo.prev;
            }
            if(cloneInfo.prev!=null){
                cloneInfo.prev.next=cloneInfo.next;
            }
        }
    }

    /**
     * records the sizes of all clonal populations that are alive when the function is called
     */
    public void RecordClonePops(){
        long[]ret=new long[nActiveClones*2];
        int iClone=0;
//        if(progenitors.pop!=0){
//            ret[0]=progenitors.id;
//            ret[1]=progenitors.pop;
//            iClone++;
//        }
        T clone= livingCloneInfos;
        if(livingCloneInfos!=null && livingCloneInfos.prev!=null){
            System.out.println("here?");
        }
        while(clone !=null){
            ret[iClone*2]=clone.id;
            ret[iClone*2+1]=clone.pop;
            clone=(T)clone.next;
            iClone++;
        }
        CloneCounts.add(ret);
    }

    /**
     * writes out the parent ids to a file
     * @param parentIDsOut the file to write to
     * @param delim delimeter to separate each id entry
     */
    public void WriteParentIDs(FileIO parentIDsOut,String delim){
        for (int id : parentIDs) {
            if(parentIDsOut.binary){
                parentIDsOut.WriteBinInt(id);
            }
            else{
                parentIDsOut.Write(id+delim);
            }
        }
    }

    /**
     * writes the mutation info strings to a file
     * @param mutationInfoOut the file to write to
     * @param delim delimeter to separate each info entry
     */
    public void WriteMutationInfo(FileIO mutationInfoOut,String delim){
        for(String info: allGenomeInfos){
            if(mutationInfoOut.binary){
                mutationInfoOut.WriteBinString(info);
            }
            else{
                mutationInfoOut.Write(info+delim);
            }
        }
    }

    /**
     * writes the mutation info of all clone types that are currently alive in the form id delim info delim id delim info delim... to a file
     * @param mutationInfoOut file to write to
     * @param delim delimiter to separate each entry
     */
    public void WriteMutationInfoLiving(FileIO mutationInfoOut,String delim){
        if(allGenomeInfos==null){
            throw new IllegalStateException("Mutant to String function was not defined in GenomeTracker constructor");
        }
        T curr=livingCloneInfos;
        while(curr!=null) {
            mutationInfoOut.Write(curr.id+delim+MutantToString.MutantToString(curr)+delim);
            curr=(T)curr.next;
        }
    }

    /**
     * gets a concatenated string of the full lineage of the particular clone in the form info delim info delim...
     * @param id id of genome of interest
     * @param delim used to separate genome strings along lineage
     */
    public String FullLineageInfoStr(int id, String delim){
        StringBuilder sb=new StringBuilder();
        while(id!=-1){
            sb.append(allGenomeInfos.get(id)+delim);
            id=parentIDs.get(id);
        }
        return sb.toString();
    }

    /**
     * writes the full lineage of all clone types that are currently alive in the form id outerdelim info innerdelim info innerdelim... outerdelim id outerdelim info...
     * @param lineageInfoOut the file to write to
     * @param innerDelim used to separate info entries
     * @param outerDelim used to separate ids from info
     */
    public void WriteAllLineageInfoLiving(FileIO lineageInfoOut,String innerDelim,String outerDelim){
        T curr=livingCloneInfos;
        while(curr!=null){
            if(lineageInfoOut.binary){
                lineageInfoOut.WriteBinString(curr.id+innerDelim+ FullLineageInfoStr(curr.id,innerDelim)+outerDelim);
            }
            else {
                lineageInfoOut.Write(curr.id+innerDelim+ FullLineageInfoStr(curr.id, innerDelim) + outerDelim);
            }
            curr=(T)curr.next;
        }
    }

    /**
     * writes the full lineage of all clone types that have ever existed in the form outerdelim info innerdelim info innerdelim... outerdelim info...
     * @param lineageInfoOut the file to write to
     * @param innerDelim used to separate info entries
     * @param outerDelim used to separate lineages
     */
    public void WriteAllLineageInfo(FileIO lineageInfoOut,String innerDelim,String outerDelim) {
        for (int i = 0; i < parentIDs.size(); i++) {
            if(lineageInfoOut.binary){
                lineageInfoOut.WriteBinString(i+innerDelim+FullLineageInfoStr(i,innerDelim)+outerDelim);
            }
            else{
                lineageInfoOut.Write(i+innerDelim+FullLineageInfoStr(i,innerDelim)+outerDelim);
            }
        }
    }

    /**
     * writes all clone population counts stored to a file in the form outerdelim id innerdelim info innerdelim id innerdelim info... outerdelim...
     * @param cloneCountsOut file to write clonecounts to
     * @param innerDelim used to separate ids and infos
     * @param outerDelim used to separate clonecount events
     */
    public void WriteClonePops(FileIO cloneCountsOut,String innerDelim,String outerDelim){
        for (long[] counts : CloneCounts) {
            cloneCountsOut.WriteDelimit(counts,innerDelim);
            cloneCountsOut.Write(outerDelim);
        }
    }
    /**
     * writes all clone population counts stored to a file in the form outerdelim id innerdelim info innerdelim id innerdelim info... outerdelim...
     * @param cloneCountsOut file to write clonecounts to
     * @param innerDelim used to separate ids and infos
     * @param outerDelim used to separate clonecount events
     */
    public void WriteClonePopsLineageCumulative(FileIO cloneCountsOut,String innerDelim,String outerDelim){
        HashMap<Integer,Integer> popCounts;
        for (long[] counts : CloneCounts) {
            for(int i=0;i<counts.length;i+=2){
                long id=counts[i*2];
                long pop=counts[i*2+1];
                while(id!=-1){

                }
            }
            cloneCountsOut.WriteDelimit(counts,innerDelim);
            cloneCountsOut.Write(outerDelim);
        }
    }

    @Override
    public Iterator<T> iterator() {
        return null;
    }
    private class myIter implements Iterator<T>{
        T curr;

        @Override
        public boolean hasNext() {
            return curr.next!=null;
        }

        @Override
        public T next() {
            curr=(T)curr.next;
            return curr;
        }
    }
}
