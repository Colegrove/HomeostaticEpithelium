package Framework.Gui;

import Framework.Grids.*;
import Framework.Misc.AgentToColorInt;
import Framework.Misc.DoubleToIntFunction;
import Framework.Tools.Utils;
import Framework.Misc.GuiComp;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.File; import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

class PaintPanel extends JPanel{
    final BufferedImage buff;
    final int xDim;
    final int yDim;
    public final int scale;
    Graphics2D g;
    PaintPanel(BufferedImage buff,int xDim,int yDim,int scale){
        this.buff = buff;
        this.xDim = xDim;
        this.yDim = yDim;
        this.scale = scale;
    }
    @Override
    public void paint(Graphics g){
        ((Graphics2D)g).drawImage(buff.getScaledInstance(scale*xDim,-scale*yDim,Image.SCALE_FAST),null,null);
        repaint();
    }

}

/**
 * a gui item that is used to efficiently visualize in 2 dimensions
 * uses an array of pixels whose color values are individually set
 */
public class GuiGridVis implements GuiComp {
    boolean active;
    final PaintPanel panel;
    final public int xDim;
    final public int yDim;
    public final int scale;
    public final int compX;
    public final int compY;
    final BufferedImage buff;
    final int[] data;

    /**
     * @param gridW width of the GuiGridVis in pixels
     * @param gridH height of the GuiGridVis in pixels
     * @param scaleFactor the width and height in screen pixels of each GuiGridVis pixel
     * @param compX width on the gui GridBagLayout
     * @param compY height on the gui GridBagLayout
     * @param active
     */
    public GuiGridVis(int gridW, int gridH, int scaleFactor, int compX, int compY, boolean active){
        this.active= active;
        xDim=gridW;
        yDim=gridH;
        this.compX=compX;
        this.compY=compY;
        scale=scaleFactor;
        if(this.active) {
            buff=new BufferedImage(xDim,yDim,BufferedImage.TYPE_INT_RGB);
            panel=new PaintPanel(buff,xDim,yDim,scale);
            panel.setVisible(true);
            panel.setPreferredSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMaximumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMinimumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            data=((DataBufferInt)buff.getRaster().getDataBuffer()).getData();
        }
        else{
           panel=null;
           buff=null;
           data=null;
        }
    }
    public GuiGridVis(int gridW, int gridH, int scaleFactor, int compX, int compY){
        this.active= true;
        xDim=gridW;
        yDim=gridH;
        this.compX=compX;
        this.compY=compY;
        scale=scaleFactor;
        if(this.active) {
            buff=new BufferedImage(xDim,yDim,BufferedImage.TYPE_INT_RGB);
            panel=new PaintPanel(buff,xDim,yDim,scale);
            panel.setVisible(true);
            panel.setPreferredSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMaximumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMinimumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            data=((DataBufferInt)buff.getRaster().getDataBuffer()).getData();
        }
        else{
            panel=null;
            buff=null;
            data=null;
        }
    }

    /**
     * @param gridW width of the GuiGridVis in pixels
     * @param gridH height of the GuiGridVis in pixels
     * @param scaleFactor the width and height in screen pixels of each GuiGridVis pixel
     */
    public GuiGridVis(int gridW, int gridH, int scaleFactor,boolean active){
        this.active=active;
        xDim=gridW;
        yDim=gridH;
        this.compX=1;
        this.compY=1;
        scale=scaleFactor;
        if(active) {
            buff=new BufferedImage(xDim,yDim,BufferedImage.TYPE_INT_RGB);
            panel=new PaintPanel(buff,xDim,yDim,scale);
            panel.setVisible(true);
            panel.setPreferredSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMaximumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMinimumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            data=((DataBufferInt)buff.getRaster().getDataBuffer()).getData();
        }
        else{
            panel=null;
            buff=null;
            data=null;
        }
    }
    public GuiGridVis(int gridW, int gridH, int scaleFactor){
        this.active=true;
        xDim=gridW;
        yDim=gridH;
        this.compX=1;
        this.compY=1;
        scale=scaleFactor;
        if(active) {
            buff=new BufferedImage(xDim,yDim,BufferedImage.TYPE_INT_RGB);
            panel=new PaintPanel(buff,xDim,yDim,scale);
            panel.setVisible(true);
            panel.setPreferredSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMaximumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            panel.setMinimumSize(new Dimension((int)Math.ceil(xDim*scaleFactor), (int)Math.ceil(yDim*scaleFactor)));
            data=((DataBufferInt)buff.getRaster().getDataBuffer()).getData();
        }
        else{
            panel=null;
            buff=null;
            data=null;
        }
    }
    /**
     * sets pixel at the specified xDim,yDim position to the rgb color value specified, bounding components between 0 and 1
     */
    public void SetColorBound(int x, int y, float r, float g, float b){
        r=Utils.BoundValF(r,0,1);
        g=Utils.BoundValF(g,0,1);
        b=Utils.BoundValF(b,0,1);
        SetColor(x,y,r,g,b);
    }
//    public void DrawStamp(float[] stampColors,double xMid,double yMid,int xSize,int ySize){
//        int xStart=
//    }

    /**
     * sets the pixel at the specified xDim,yDim position using the heat colormap, which goes from black to red to yellow to white
     */
    public void SetColorHeat(int x, int y, double val) {
        if(active) {
            float r = (float) Math.min(1, val * 4);
            float g = 0;
            float b = 0;
            if (val > 0.25) {
                g = (float) Math.min(1, (val - 0.25) * 2);
            }
            if (val > 0.75) {
                b = (float) Math.min(1, (val - 0.75) * 4);
            }
            SetColor(x, y, r, g, b);
        }
    }

    /**
     * sets the pixel at the specified xDim,yDim position using the heat colormap, which goes from black to red to yellow to white
     * @param colorOrder the order in which to fill the colors, default is rgb, but any permutation of these three characters is valid
     */
    public void SetColorHeat(int x, int y, double val,String colorOrder) {
        if(active) {
            float c1 = (float) Math.min(1, val * 4);
            float c2 = 0;
            float c3 = 0;
            if (val > 0.25) {
                c2 = (float) Math.min(1, (val - 0.25) * 2);
            }
            if (val > 0.75) {
                c3 = (float) Math.min(1, (val - 0.75) * 4);
            }
            switch (colorOrder) {
                case "rgb":
                    SetColor(x, y, c1, c2, c3);
                    break;
                case "rbg":
                    SetColor(x, y, c1, c3, c2);
                    break;
                case "grb":
                    SetColor(x, y, c2, c1, c3);
                    break;
                case "gbr":
                    SetColor(x, y, c3, c1, c2);
                    break;
                case "brg":
                    SetColor(x, y, c2, c3, c1);
                    break;
                case "bgr":
                    SetColor(x, y, c3, c2, c1);
                    break;
                default:
                    SetColor(x, y, c1, c2, c3);
                    System.out.println("Invalid colorOrder string passed to SetColorHeat:" + colorOrder + "\ncolorOrder String must be some permutation of the characters 'r','g','b'");
                    break;
            }
        }
    }
    public void SetColorHeatBound(int x,int y,double val){
        if(active) {
            val = Utils.BoundVal(val, 0, 1);
            SetColorHeat(x, y, val);
        }
    }
    public void SetColorHeatBound(int x,int y,double val,String colorOrder){
        if(active) {
            val = Utils.BoundVal(val, 0, 1);
            SetColorHeat(x, y, val, colorOrder);
        }
    }
    /** * Converts HSV values to RGB values
     * hue: 0->1 (Picks Color)
     * saturation: 0->1
     * brightness: 0->1
     **/
    public void SetColorHSV(int x,int y,float hue,float saturation,float brightness){
        if(active) {
            float scale = 255f;
            int rgb = Color.HSBtoRGB(hue, saturation, brightness);
            data[(yDim - y - 1) * xDim + x] = rgb;
        }
    }

    /**
     * gets the xDim component of the vis window
     */
    @Override
    public int compX(){return compX;}
    /**
     * gets the yDim component of the vis window
     */
    @Override
    public int compY(){return compY;}

    @Override
    public boolean IsActive() {
        return active;
    }

    @Override
    public void SetActive(boolean isActive) {
        this.active=isActive;
    }

    /**
     * called by the GuiWindow class to place the vis window
     */
    @Override
    public void GetComps(ArrayList<Component> putHere, ArrayList<Integer> coordsHere, ArrayList<Integer> compSizesHere) {
        putHere.add(panel);
        coordsHere.add(0);
        coordsHere.add(0);
        compSizesHere.add(compX);
        compSizesHere.add(compY);
    }
    public void SetColor(int i, double r, double g, double b){
        if(active) {
            SetColor(i / yDim, i % yDim, r, g, b);
        }
    }
    public void SetColor(int x, int y, double r, double g, double b){
        if(active) {
            int ri = (int) (r * 255);
            int gi = (int) (g * 255);
            int bi = (int) (b * 255);
            data[(yDim - y - 1) * xDim + x] = (((ri & 0x0ff) << 16) | ((gi & 0x0ff) << 8) | (bi & 0x0ff)) - 16777216;
        }
    }
    public void SetColor(int x,int y,int colorInt){
        if(active) {
            data[(yDim - y - 1) * xDim + x] = colorInt;
        }
    }
    public void SetColor(int i,int colorInt){
        if(active) {
            SetColor(i / yDim, i % yDim, colorInt);
        }
    }
    public void SetColor(int x, int y, Color color) {
        if(active) {
            data[(yDim - y - 1) * xDim + x] = color.getRGB();
        }
    }
    public void SetColor(int i, Color color) {
        if(active) {
            SetColor(i / yDim, i % yDim, color);
        }
    }

    public double ClickXpt(MouseEvent e){
        return e.getX()*1.0/scale;
    };
    public int ClickXsq(MouseEvent e){
        return e.getX()/scale;
    };
    public double ClickYpt(MouseEvent e){
        return (yDim-1)-e.getY()*1.0/scale;
    };
    public int ClickYsq(MouseEvent e){
        return (yDim-1)-e.getY()/scale;
    };

    /**
     * sets all pixels to the rgb color specified, bounding components between 0 and 1
     */
    public void ClearColor(double r, double g, double b){
        if(active) {
            Arrays.fill(data, Utils.ColorInt(r,g,b));
        }
    }

    public void DrawGridDiff(GridDiff2 drawMe,double min,double max,String colorOrder){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    SetColorHeat(x, y, (drawMe.GetCurr(x, y) - min) / range, colorOrder);
                }
            }
        }
    }
    public void DrawGridDiff(GridDiff2 drawMe,double min,double max){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    SetColorHeat(x, y, (drawMe.GetCurr(x, y) - min) / range);
                }
            }
        }
    }
    public void DrawGridDiff(GridDiff2 drawMe,DoubleToIntFunction ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    SetColor(x, y, ColorFn.DoubleToInt(drawMe.GetCurr(x,y)));
                }
            }
        }
    }
    public void DrawGridDiffBound(GridDiff2 drawMe,double min,double max){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    SetColorHeatBound(x, y, (drawMe.GetCurr(x, y) - min) / range);
                }
            }
        }
    }
    public void DrawGridDiffBound(GridDiff2 drawMe,double min,double max,String colorOrder){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    SetColorHeatBound(x, y, (drawMe.GetCurr(x, y) - min) / range, colorOrder);
                }
            }
        }
    }
    public void DrawGridDiffXY(GridDiff3 drawMe, double min, double max, String colorOrder){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    double sum = 0;
                    for (int z = 0; z < drawMe.zDim; z++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(x, y, sum / drawMe.zDim, colorOrder);
                }
            }
        }
    }
    public void DrawGridDiffYZ(GridDiff3 drawMe, double min, double max, String colorOrder){
        if(active) {
            double range = max - min;
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int x = 0; x < drawMe.xDim; x++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(y, z, sum / drawMe.xDim, colorOrder);
                }
            }
        }
    }
    public void DrawGridDiffXZ(GridDiff3 drawMe, double min, double max, String colorOrder){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int y = 0; y < drawMe.yDim; y++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(x, z, sum / drawMe.xDim, colorOrder);
                }
            }
        }
    }
    public void DrawGridDiffXY(GridDiff3 drawMe, DoubleToIntFunction ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    double sum = 0;
                    for (int z = 0; z < drawMe.zDim; z++) {
                        sum += drawMe.GetCurr(x, y, z);
                    }
                    SetColor(x, y, ColorFn.DoubleToInt(sum/drawMe.zDim));
                }
            }
        }
    }
    public void DrawGridDiffYZ(GridDiff3 drawMe, DoubleToIntFunction ColorFn){
        if(active) {
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int x = 0; x < drawMe.xDim; x++) {
                        sum += drawMe.GetCurr(x, y, z);
                    }
                    SetColor(y, z, ColorFn.DoubleToInt(sum/drawMe.zDim));
                }
            }
        }
    }
    public void DrawGridDiffXZ(GridDiff3 drawMe, DoubleToIntFunction ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int y = 0; y < drawMe.yDim; y++) {
                        sum += drawMe.GetCurr(x, y, z);
                    }
                    SetColor(x, z, ColorFn.DoubleToInt(sum / drawMe.xDim));
                }
            }
        }
    }
    public void DrawGridDiffXY(GridDiff3 drawMe, double min, double max){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    double sum = 0;
                    for (int z = 0; z < drawMe.zDim; z++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(x, y, sum / drawMe.zDim);
                }
            }
        }
    }
    public void DrawGridDiffYZ(GridDiff3 drawMe, double min, double max){
        if(active) {
            double range = max - min;
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int x = 0; x < drawMe.xDim; x++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(y, z, sum / drawMe.xDim);
                }
            }
        }
    }
    public void DrawGridDiffXZ(GridDiff3 drawMe, double min, double max){
        if(active) {
            double range = max - min;
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    double sum = 0;
                    for (int y = 0; y < drawMe.yDim; y++) {
                        sum += (drawMe.GetCurr(x, y, z) - min) / range;
                    }
                    SetColorHeat(x, z, sum / drawMe.xDim);
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid2<Q>>void DrawAgents(T drawMe, AgentToColorInt<Q> ColorFn, double bkR, double bkG, double bkB){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    Q a = drawMe.GetAgent(x, y);
                    if (a != null) {
                        SetColor(x, y, ColorFn.AgentToColor(a));
                    } else {
                        SetColor(x, y, bkR, bkG, bkB);
                    }
                }
            }
        }
    }
//    public <Q extends AgentBaseSpatial,T extends Grid2<Q>>void DrawAgentDensity(T drawMe, int maxDensity,String colorOrder){
//        if(active) {
//            for (int x = 0; x < drawMe.xDim; x++) {
//                for (int y = 0; y < drawMe.yDim; y++) {
//                    //SetColorHeatBound(x,y,drawMe.PopAt(x,y)*1.0/maxDensity,colorOrder);
//                }
//            }
//        }
//        throw new IllegalStateException("this function should be deleted");
//    }
//    public <Q extends AgentBaseSpatial,T extends Grid2<Q>>void DrawAgentDensity(T drawMe, int maxDensity){
//        if(active) {
//            for (int x = 0; x < drawMe.xDim; x++) {
//                for (int y = 0; y < drawMe.yDim; y++) {
//                    //SetColorHeatBound(x,y,drawMe.PopAt(x,y)*1.0/maxDensity);
//                }
//            }
//        }
//        throw new IllegalStateException("this function should be deleted");
//    }
    public <Q extends AgentBaseSpatial,T extends Grid2<Q>>void DrawAgents(T drawMe, AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    Q a = drawMe.GetAgent(x, y);
                    if (a != null) {
                        SetColor(x, y, ColorFn.AgentToColor(a));
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsXY(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    for (int z = 0; z <= drawMe.zDim; z++) {
                        if (z == drawMe.zDim) {
                            SetColor(x, y, bkColor);
                        } else {
                            Q a = drawMe.GetAgent(x, y, z);
                            if (a != null) {
                                SetColor(x, y, ColorFn.AgentToColor(a));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsXZ(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int y = 0; y <= drawMe.yDim; y++) {
                        if (y == drawMe.yDim) {
                            SetColor(x, z, bkColor);
                        } else {
                            Q a = drawMe.GetAgent(x, y, z);
                            if (a != null) {
                                SetColor(x, z, ColorFn.AgentToColor(a));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsYZ(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int x = 0; x <= drawMe.xDim; x++) {
                        if (x == drawMe.yDim) {
                            SetColor(y, z, bkColor);
                        } else {
                            Q a = drawMe.GetAgent(x, y, z);
                            if (a != null) {
                                SetColor(y, z, ColorFn.AgentToColor(a));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsXY(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    for (int z = drawMe.zDim-1; z >= -1; z++) {
                        if (z == -1) {
                            SetColor(x, y, bkColor);
                        } else {
                            Q a = drawMe.GetAgent(x, y, z);
                            if (a != null) {
                                SetColor(x, y, ColorFn.AgentToColor(a));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsXZ(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int y = drawMe.yDim-1; y >= -1; y++) {
                        if (y == -1) {
                            SetColor(x, z, bkColor);
                        } else {
                            Q a = drawMe.GetAgent(x, y, z);
                            if (a != null) {
                                SetColor(x, z, ColorFn.AgentToColor(a));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsYZ(T drawMe,AgentToColorInt<Q> ColorFn,double bkR,double bkG,double bkB){
        if(active) {
            int bkColor = Utils.ColorInt(bkR, bkG, bkB);
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int x = drawMe.xDim-1; x >= 0; x++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(y, z, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsXY(T drawMe,AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    for (int z = 0; z < drawMe.zDim; z++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(x, y, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsXZ(T drawMe,AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int y = 0; y < drawMe.yDim; y++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(x, z, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawClosestAgentsYZ(T drawMe,AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int x = 0; x < drawMe.xDim; x++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(y, z, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsXY(T drawMe,AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int y = 0; y < drawMe.yDim; y++) {
                    for (int z = drawMe.zDim - 1; z >= 0; z++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(x, y, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsXZ(T drawMe,AgentToColorInt<Q> ColorFn){
        if(active) {
            for (int x = 0; x < drawMe.xDim; x++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int y = drawMe.yDim - 1; y >= 0; y++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(x, z, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawFurthestAgentsYZ(T drawMe,AgentToColorInt<Q> ColorFn) {
        if (active) {
            for (int y = 0; y < drawMe.yDim; y++) {
                for (int z = 0; z < drawMe.zDim; z++) {
                    for (int x = drawMe.xDim - 1; x >= 0; x++) {
                        Q a = drawMe.GetAgent(x, y, z);
                        if (a != null) {
                            SetColor(y, z, ColorFn.AgentToColor(a));
                            break;
                        }
                    }
                }
            }
        }
    }
//    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawAgentDensityXY(T drawMe,int maxDensity,String colorOrder){
//        int[]densities=new int[drawMe.xDim*drawMe.yDim];
//        for (Q agent : drawMe) {
//            int i=agent.Isq();
//        }
//        if(active) {
//            for (int x = 0; x < drawMe.xDim; x++) {
//                for (int y = 0; y < drawMe.yDim; y++) {
//                    int sum = 0;
//                    for (int z = 0; z < drawMe.zDim; z++) {
////                        sum += drawMe.PopAt(x, y, z);
//                    }
//                    SetColorHeatBound(x, y, sum * 1.0 / maxDensity,colorOrder);
//                }
//            }
//        }
//    }
//    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawAgentDensityXZ(T drawMe,int maxDensity,String colorOrder){
//        if(active) {
//            for (int x = 0; x < drawMe.xDim; x++) {
//                for (int z = 0; z < drawMe.zDim; z++) {
//                    int sum = 0;
//                    for (int y = 0; y < drawMe.yDim; y++) {
////                        sum += drawMe.PopAt(x, y, z);
//                    }
//                    SetColorHeatBound(x, z, sum * 1.0 / maxDensity,colorOrder);
//                }
//            }
//        }
//    }
//    public <Q extends AgentBaseSpatial,T extends Grid3<Q>>void DrawAgentDensityYZ(T drawMe,int maxDensity,String colorOrder){
//        if(active) {
//            for (int y = 0; y < drawMe.yDim; y++) {
//                for (int z = 0; z < drawMe.zDim; z++) {
//                    int sum = 0;
//                    for (int x = 0; x < drawMe.xDim; x++) {
////                        sum += drawMe.PopAt(x, y, z);
//                    }
//                    SetColorHeatBound(y, z, sum * 1.0 / maxDensity,colorOrder);
//                }
//            }
//        }
//    }

    /**
     * called by the GuiWindow to draw the vis
     */
    public void ToPNG(String path){
        if(active) {
            File out = new File(path);
            try {
                ImageIO.write(buff, "png", out);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
    public void ToJPG(String path){
        if(active) {
            File out = new File(path);
            try {
                ImageIO.write(buff, "jpg", out);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
    public void ToGIF(String path) {
        if (active) {
            File out = new File(path);
            try {
                ImageIO.write(buff, "gif", out);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
}
