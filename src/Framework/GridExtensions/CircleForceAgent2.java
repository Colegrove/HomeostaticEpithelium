package Framework.GridExtensions;

import Framework.Grids.AgentPT2;
import Framework.Grids.Grid2;

import java.util.ArrayList;
import java.util.Random;

import static Framework.Tools.Utils.*;

/**
 * Created by bravorr on 6/26/17.
 */
public class CircleForceAgent2<T extends Grid2> extends AgentPT2<T> {
    public double radius;
    public double xVel;
    public double yVel;
    public <Q extends CircleForceAgent2> double SumForces(double interactionRad, ArrayList<Q> scratchAgentList, OverlapForceResponse OverlapFun, boolean wrapX, boolean wrapY){
        scratchAgentList.clear();
        double sum=0;
        G().AgentsInRad(scratchAgentList,Xpt(),Ypt(),interactionRad,wrapX,wrapY);
        for (Q a : scratchAgentList) {
            if(a!=this){
                double xComp=Xdisp(a,wrapX);
                double yComp=Ydisp(a,wrapY);
                if(xComp==0&&yComp==0){
                    xComp=Math.random()-0.5;
                    yComp=Math.random()-0.5;
                }
                double dist=Norm(xComp,yComp);
                if(dist<interactionRad) {
                    double touchDist = (radius + a.radius) - dist;
                    double force=OverlapFun.CalcForce(touchDist);
                    xVel+=(xComp/dist)*force;
                    yVel+=(yComp/dist)*force;
                    if(force>0) {
                        sum += Math.abs(force);
                    }
                }
            }
        }
        return sum;
    }
    public <Q extends CircleForceAgent2> double SumForces(double interactionRad, ArrayList<Q> scratchAgentList, OverlapForceResponse OverlapFun){
        return SumForces(interactionRad, scratchAgentList, OverlapFun, G().wrapX, G().wrapY);
    }

    public void ForceMove(double friction,boolean wrapX,boolean wrapY){
        xVel*=friction;
        yVel*=friction;
        MoveSafePT(Xpt()+xVel,Ypt()+yVel,wrapX,wrapY);
    }
    public void ForceMove(double friction){
        xVel*=friction;
        yVel*=friction;
        MoveSafePT(Xpt()+xVel,Ypt()+yVel,G().wrapX,G().wrapY);
    }
    public <Q extends CircleForceAgent2> Q Divide(double divRadius, double[] scratchCoordArr, Random rn, boolean wrapX, boolean wrapY){
        if(rn!=null){
            RandomPointOnCircleEdge(divRadius,rn,scratchCoordArr);
        }
        Q child=(Q)(G().NewAgentPTSafe(Xpt()+scratchCoordArr[0],Ypt()+scratchCoordArr[1],Xpt(),Ypt(),wrapX,wrapY));
        MoveSafePT(Xpt()-scratchCoordArr[0],Ypt()-scratchCoordArr[1],wrapX,wrapY);
        return child;
    }
    public <Q extends CircleForceAgent2> Q Divide(double divRadius, double[] scratchCoordArr, Random rn){
        if(rn!=null){
            RandomPointOnCircleEdge(divRadius,rn,scratchCoordArr);
        }
        Q child=(Q)(G().NewAgentPTSafe(Xpt()+scratchCoordArr[0],Ypt()+scratchCoordArr[1],Xpt(),Ypt(),G().wrapX,G().wrapY));
        MoveSafePT(Xpt()-scratchCoordArr[0],Ypt()-scratchCoordArr[1],G().wrapX,G().wrapY);
        return child;
    }
}
