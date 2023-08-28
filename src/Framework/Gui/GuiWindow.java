package Framework.Gui;

import Framework.Misc.GuiCloseAction;
import Framework.Misc.GuiComp;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;

/**
 * a JFrame wrapper that holds menu and gui items
 * the window that the other gui items sit in
 */
public class GuiWindow{
    boolean alive;
    final boolean active;
    final boolean main;
    final private GridBagConstraints gbc;
    public final JFrame frame;
    public final JPanel panel;
    final ArrayList<GuiComp> comps;
    final ArrayList<Integer> compCoords;
    final ArrayList<Component> subComps;
    final ArrayList<Integer> subCompCoords;
    final ArrayList<Integer> subCompSizes;
    final int[] locs;
    final GuiCloseAction closeAction;

    /**
     * @param title the title that will appear at the top of the window
     * @param main whether the program should terminate on closing the window
     * @param closeAction function that will run when the window is closed
     * @param active
     */
    public GuiWindow(String title, boolean main, GuiCloseAction closeAction, boolean active){
        this.active=active;
        this.main=main;
        if(active) {
            this.frame = new JFrame();
            this.comps = new ArrayList<>();
            this.compCoords = new ArrayList<>();
            this.subComps = new ArrayList<>();
            this.subCompCoords = new ArrayList<>();
            this.subCompSizes = new ArrayList<>();
            this.frame.setResizable(false);//fixes window size
            this.frame.setLocationRelativeTo(null);//puts window in middle of screen
            this.closeAction = closeAction;
            if (main) {
                this.frame.addWindowListener(new WindowAdapter() {
                    @Override
                    public void windowClosing(WindowEvent e) {
                        if (closeAction != null) {
                            closeAction.Action(e);
                        }
                        frame.setVisible(false);
                        frame.dispose();
                        System.exit(0);
                    }
                });
            } else {
                this.frame.addWindowListener(new WindowAdapter() {
                    @Override
                    public void windowClosing(WindowEvent e) {
                        if (alive) {
                            alive = false;
                            if (closeAction != null) {
                                closeAction.Action(e);
                            }
                            frame.setVisible(false);
                            frame.dispose();
                        }
                    }
                });
            }
            this.frame.setTitle(title);
            panel = new JPanel();
            gbc = new GridBagConstraints();
            panel.setLayout(new GridBagLayout());
            this.frame.add(panel);
            this.locs = new int[1000 * 1000];
            for (int i = 0; i < 1000 * 1000; i++) {
                locs[i] = -1;
            }
        }
        else{
            this.frame=null;
            this.comps=null;
            this.compCoords=null;
            this.subComps=null;
            this.subCompCoords=null;
            this.panel=null;
            this.locs=null;
            this.subCompSizes=null;
            this.closeAction=null;
            this.gbc=null;
        }
    }
    public GuiWindow(String title, boolean main, GuiCloseAction closeAction){
        this.active=true;
        this.main=main;
        if(active) {
            this.frame = new JFrame();
            this.comps = new ArrayList<>();
            this.compCoords = new ArrayList<>();
            this.subComps = new ArrayList<>();
            this.subCompCoords = new ArrayList<>();
            this.subCompSizes = new ArrayList<>();
            this.frame.setResizable(false);//fixes window size
            this.frame.setLocationRelativeTo(null);//puts window in middle of screen
            this.closeAction = closeAction;
            if (main) {
                this.frame.addWindowListener(new WindowAdapter() {
                    @Override
                    public void windowClosing(WindowEvent e) {
                        if (closeAction != null) {
                            closeAction.Action(e);
                        }
                        frame.setVisible(false);
                        frame.dispose();
                        System.exit(0);
                    }
                });
            } else {
                this.frame.addWindowListener(new WindowAdapter() {
                    @Override
                    public void windowClosing(WindowEvent e) {
                        if (alive) {
                            alive = false;
                            if (closeAction != null) {
                                closeAction.Action(e);
                            }
                            frame.setVisible(false);
                            frame.dispose();
                        }
                    }
                });
            }
            this.frame.setTitle(title);
            panel = new JPanel();
            gbc = new GridBagConstraints();
            panel.setLayout(new GridBagLayout());
            this.frame.add(panel);
            this.locs = new int[1000 * 1000];
            for (int i = 0; i < 1000 * 1000; i++) {
                locs[i] = -1;
            }
        }
        else{
            this.frame=null;
            this.comps=null;
            this.compCoords=null;
            this.subComps=null;
            this.subCompCoords=null;
            this.panel=null;
            this.locs=null;
            this.subCompSizes=null;
            this.closeAction=null;
            this.gbc=null;
        }
    }
    /**
     * @param title the title that will appear at the top of the window
     * @param main whether the program should terminate on closing the window
     */
    public GuiWindow(String title, boolean main,boolean active) {
        this.active=active;
        this.main=main;
        if(active) {
            this.frame = new JFrame();
            this.comps = new ArrayList<>();
            this.compCoords = new ArrayList<>();
            this.subComps = new ArrayList<>();
            this.subCompCoords = new ArrayList<>();
            this.subCompSizes = new ArrayList<>();
            this.closeAction = null;
            this.frame.setResizable(false);//fixes window size
            this.frame.setLocationRelativeTo(null);//puts window in middle of screen
            if (main) {
                this.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            } else {
                this.frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            }
            this.frame.setTitle(title);
            panel = new JPanel();
            gbc = new GridBagConstraints();
            panel.setLayout(new GridBagLayout());
            this.frame.add(panel);
            this.locs = new int[1000 * 1000];
            for (int i = 0; i < 1000 * 1000; i++) {
                locs[i] = -1;
            }
        }
        else{
            this.frame=null;
            this.comps=null;
            this.compCoords=null;
            this.subComps=null;
            this.subCompCoords=null;
            this.panel=null;
            this.locs=null;
            this.subCompSizes=null;
            this.closeAction=null;
            this.gbc=null;
        }
    }
    public GuiWindow(String title, boolean main) {
        this.active=true;
        this.main=main;
        if(active) {
            this.frame = new JFrame();
            this.comps = new ArrayList<>();
            this.compCoords = new ArrayList<>();
            this.subComps = new ArrayList<>();
            this.subCompCoords = new ArrayList<>();
            this.subCompSizes = new ArrayList<>();
            this.closeAction = null;
            this.frame.setResizable(false);//fixes window size
            this.frame.setLocationRelativeTo(null);//puts window in middle of screen
            if (main) {
                this.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            } else {
                this.frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            }
            this.frame.setTitle(title);
            panel = new JPanel();
            gbc = new GridBagConstraints();
            panel.setLayout(new GridBagLayout());
            this.frame.add(panel);
            this.locs = new int[1000 * 1000];
            for (int i = 0; i < 1000 * 1000; i++) {
                locs[i] = -1;
            }
        }
        else{
            this.frame=null;
            this.comps=null;
            this.compCoords=null;
            this.subComps=null;
            this.subCompCoords=null;
            this.panel=null;
            this.locs=null;
            this.subCompSizes=null;
            this.closeAction=null;
            this.gbc=null;
        }
    }

    /**
     * sets the gui background color
     * @param backgroundColor default color of any empty space on the gui
     */
    public void SetColor(Color backgroundColor){
        if(active) {
            this.panel.setOpaque(true);
            this.panel.setBackground(backgroundColor);
        }
    }

    /**
     * Disables or enables all interacton with the GuiWindow
     * @param onOff whether to enable or disable the gui
     */
    public void GreyOut(boolean onOff){
        if(active) {
            this.frame.setEnabled(!onOff);
            this.panel.setEnabled(!onOff);
            for (GuiComp gc : comps) {
                subComps.clear();
                subCompCoords.clear();
                subCompSizes.clear();
                gc.GetComps(subComps, subCompCoords, subCompSizes);
                for (Component sc : subComps) {
                    sc.setEnabled(!onOff);
                }
            }
            frame.repaint();
        }
    }

    void PlaceComponent(GuiComp comp,int x,int y,int w,int h){
        int iComp=comps.size();
        comps.add(comp);
        compCoords.add(x);
        compCoords.add(y);
        for(int i=x;i<x+w;i++){
            for(int j=y;j<y+h;j++){
                locs[CtoI(i,j)]=iComp;
            }
        }
    }

    /**
     * starts the gui thread and reveals the gui to the user
     */
    public void RunGui(){
        if(alive){
            throw new IllegalStateException("RunGui has already been called on this window!");
        }
        if(active) {
            alive = true;
            for (int i = 0; i < comps.size(); i++) {
                int compX = compCoords.get(i * 2);
                int compY = compCoords.get(i * 2 + 1);
                GuiComp comp = comps.get(i);
                subComps.clear();
                subCompCoords.clear();
                subCompSizes.clear();
                comp.GetComps(subComps, subCompCoords, subCompSizes);
                for (int j = 0; j < subComps.size(); j++) {
                    Component subComp = subComps.get(j);
                    int subX = subCompCoords.get(j * 2);
                    int subY = subCompCoords.get(j * 2 + 1);
                    int subW = subCompSizes.get(j * 2);
                    int subH = subCompSizes.get(j * 2 + 1);
                    AddComponent(subComp, compX + subX, compY + subY, subW, subH);
                }
            }
            this.frame.pack();
            this.frame.setVisible(true);
            panel.setVisible(true);
        }
    }

    /**
     * ingore
     */
    void AddComponent(Component comp, int x, int y, int w, int h) {
        gbc.gridx = x;
        gbc.gridy = y;
        gbc.gridwidth = w;
        gbc.gridheight = h;
        panel.add(comp,gbc);
        comp.setVisible(true);
//        if(obj instanceof VisualizerOld){
//            ((VisualizerOld)obj).createBufferStrategy(3);
//        }
        if(comp instanceof PaintPanel){
            ((PaintPanel)comp).g=(Graphics2D)(((PaintPanel)comp).getGraphics());
        }
    }
    int CtoI(int x,int y){ return x*1000+y; }
    int CtoLocVal(int x,int y){
        if(x<0||x>=1000||y<0||y>=1000){ throw new IllegalArgumentException("going for coord outside GuiWindow range"); }
        return this.locs[x*1000+y];
    }

    /**
     * adds component by dropping down into the specified column
     * think connect 4 or tetris
     * @param comp component to be added
     * @param col column to drop from. the left end of the component will occupy this column
     */
    public void AddCol(GuiComp comp, int col){
        if(!active){
            comp.SetActive(false);
        }
        if(active) {
            if (comp.IsActive()) {
                int w = comp.compX();
                int h = comp.compY();
                int found = 0;
                for (int y = 999; y >= 0; y--) {
                    for (int x = col; x < col + w; x++) {
                        if (CtoLocVal(x, y) != -1) {
                            found = y + 1;
                            break;
                        }
                    }
                    if (found != 0) {
                        break;
                    }
                }
                PlaceComponent(comp, col, found, w, h);
            }
        }
    }

    /**
     * destroys the gui with the window closing event
     */
    public void Dispose(){
        if(active) {
            frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
        }
    }
}

