package Epidermis_Model.ParamSweeper;

import Framework.Tools.FileIO;

import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

class SweepRun implements Runnable{
    //Runs parameter sweep, saves results to a file
    ParamSweeper mySweep;
    int iParamSet;
    SweepRun(ParamSweeper mySweep,int iParamSet){
        this.mySweep=mySweep;
        this.iParamSet=iParamSet;
    }

    @Override
    public void run() {
        String ret=mySweep.runner.Run(mySweep.params[iParamSet]);
        mySweep.outStrings[iParamSet]=ret;
    }
}

/**
 * ParamSweeper class uses multithreading and randomization to explore an input space
 */
public class ParamSweeper {
    ParamSweepRunFunction runner;
    ArrayList<ParamGeneratorFunction> paramGens;
    Random rn;
    double[][]params;
    String[] outStrings;
    FileIO out;

    /**
     * @param out the writer that will record model run results
     * @param runner the function that will be called with random inputs
     */
    public ParamSweeper(FileIO out, ParamSweepRunFunction runner){
        this.runner=runner;
        this.rn=new Random();
        this.paramGens=new ArrayList<>();
        this.out=out;
    }
    public void AddParam(ParamGeneratorFunction gen){
        paramGens.add(gen);
    }
    public void Sweep(int nRuns,int nThreads){
        SweepRun[] runs=new SweepRun[nRuns];
        params=new double[nRuns][paramGens.size()];
        outStrings=new String[nRuns];
        for(int i=0;i<runs.length;i++){
            runs[i]=new SweepRun(this,i);
            for(int j=0;j<paramGens.size();j++){
                params[i][j]=paramGens.get(j).GenParam(this.rn);
            }
        }
        ExecutorService exec= Executors.newFixedThreadPool(nThreads);
        for(SweepRun run:runs) {
            exec.execute(run);
        }
        exec.shutdown();
        while(!exec.isTerminated());
        //wait for all threads to finish
        for(String s:outStrings){
            out.Write(s);
        }
    }
}
