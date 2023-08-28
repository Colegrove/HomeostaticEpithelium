package Epidermis_Model;


import Framework.Gui.GuiGridVis;

/**
 * Created by schencro on 3/27/17.
 */

public class EpidermisCellVis {
    public int[] division_vis;
    public int[] division_vis_alt;
    public int[] stationary_vis;
    public int[] movement_vis;

    public void DrawCellonGridPop(GuiGridVis vis, EpidermisCell theCell){
        for(int y=0;y<5;y++){
            for(int x=0;x<5;x++){
                int cVal=stationary_vis[x+y*5];
                vis.SetColorHSV(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.myGenome.h,cVal*theCell.myGenome.s,cVal*theCell.myGenome.v);
            }
        }
    }

    public void DrawCellonGridPopZ(GuiGridVis vis, EpidermisCell theCell){
        for(int z=0;z<5;z++){
            for(int x=0;x<5;x++){
                int cVal=stationary_vis[x+z*5];
                vis.SetColorHSV(theCell.Xsq()*5+x, theCell.Zsq()*5+z,cVal*theCell.myGenome.h,cVal*theCell.myGenome.s,cVal*theCell.myGenome.v);
            }
        }
    }

    public void DrawCellonGrid(GuiGridVis vis, EpidermisCell theCell){
        switch (theCell.Action) {
            case EpidermisConst.DIVIDE:
                for(int y=0;y<5;y++){
                    for(int x=0;x<5;x++){
                        int cVal=division_vis[x+y*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*1f,cVal*0f,cVal* 0f);
                    }
                }
            break;
            case EpidermisConst.MOVING:
                for(int y=0;y<5;y++){
                    for(int x=0;x<5;x++){
                        int cVal=movement_vis[x+y*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*.5f,cVal*1f,cVal* 0f);
                    }
                }
            break;
            case EpidermisConst.STATIONARY:
                for(int y=0;y<5;y++){
                    for(int x=0;x<5;x++){
                        int cVal=stationary_vis[x+y*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*0f,cVal*0f,cVal* 1f);
                    }
                }
            break;
            default:
                System.out.println("ERROR: Unable to draw cell. Cell Action not found.");
        }
    }

    public void DrawCellonGrid3D(GuiGridVis vis, EpidermisCell theCell){
        switch (theCell.Action) {
            case EpidermisConst.DIVIDE:
                for(int z=0;z<5;z++){
                    for(int x=0;x<5;x++){
                        int cVal=division_vis[x+z*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Zsq()*5+z,cVal*1f,cVal*0f,cVal* 0f);
                    }
                }
                break;
            case EpidermisConst.MOVING:
                for(int z=0;z<5;z++){
                    for(int x=0;x<5;x++){
                        int cVal=movement_vis[x+z*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Zsq()*5+z,cVal*.5f,cVal*1f,cVal* 0f);
                    }
                }
                break;
            case EpidermisConst.STATIONARY:
                for(int z=0;z<5;z++){
                    for(int x=0;x<5;x++){
                        int cVal=stationary_vis[x+z*5];
                        //vis.SetColor(theCell.Xsq()*5+x, theCell.Ysq()*5+y,cVal*theCell.h,cVal*theCell.v,cVal*theCell.s);
                        vis.SetColor(theCell.Xsq()*5+x, theCell.Zsq()*5+z,cVal*0f,cVal*0f,cVal* 1f);
                    }
                }
                break;
            default:
                System.out.println("ERROR: Unable to draw cell. Cell Action not found.");
        }
    }

    public void DrawEmptyCell(GuiGridVis vis, int x, int y){
        for(int i=0; i<5; i++){
            for(int j=0; j<5; j++){
                vis.SetColor(i+x*5, j+y*5, 0f, 0f, 0f);
            }
        }
    }

    /**
     * Used to visualize cell status for those that divided in that time step
     **/
    public EpidermisCellVis() {
        division_vis = new int[]{
                0,1,1,1,0,
                1,0,0,0,1,
                1,0,0,0,1,
                1,0,0,0,1,
                0,1,1,1,0
        };


        /** Used to visualize cell status for those that are alive, but are quiescent no activity in that time step **/
        stationary_vis = new int[]{
                1,1,1,1,1,
                1,1,1,1,1,
                1,1,1,1,1,
                1,1,1,1,1,
                1,1,1,1,1}; // *****

        /** Used to visualize cell status for those that moved, but are not dividing **/
        movement_vis = new int[]{
                1,1,1,1,1, // *****
                1,0,1,0,1, // *-*-*
                1,1,1,1,1, // *****
                1,0,1,0,1,// *-*-*
                1,1,1,1,1}; // *****
}
}
