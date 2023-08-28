package Framework.GridExtensions;

import Framework.Misc.GuiComp;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created by Rafael on 8/16/2017.
 */
public class GuiColorBar implements GuiComp{
    //creates a vertical colorbar
    GuiColorBar(int pixX,int pixY,int[] labels,int compX,int compY,boolean active){
        if(active){
            if(compX<2){
                throw new IllegalArgumentException();
            }
        }
    }
    @Override
    public int compX() {
        return 0;
    }

    @Override
    public int compY() {
        return 0;
    }

    @Override
    public boolean IsActive() {
        return false;
    }

    @Override
    public void SetActive(boolean isActive) {

    }

    @Override
    public void GetComps(ArrayList<Component> putHere, ArrayList<Integer> coordsHere, ArrayList<Integer> compSizesHere) {

    }
}
