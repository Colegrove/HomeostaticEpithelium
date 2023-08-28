package Framework.Misc;

import java.util.ArrayList;

/**
 * Created by rafael on 4/15/17.
 */
public class SweepRun<T> implements Runnable{
    //Runs parameter sweep, saves results to a file
    final int iParamSet;
    final SweepRunFunction<T> runner;
    final ArrayList<T> runOut;
    public SweepRun(SweepRunFunction<T> runner, ArrayList<T> runOut, int iParamSet){
        this.iParamSet=iParamSet;
        this.runner=runner;
        this.runOut=runOut;
    }

    @Override
    public void run() {
        runOut.set(iParamSet,runner.Run(iParamSet));
    }
}
