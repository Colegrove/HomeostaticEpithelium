package Framework.Gui;

import Framework.Misc.*;

import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import static Framework.Tools.Utils.*;


/**
 * a menu item that takes double input
 */
public class DoubleParam extends JFormattedTextField implements Framework.Misc.MenuItem,GuiComp {
    static final DecimalFormat doubleFormat=new DecimalFormat("#.##########################");
    String labelText;
    JLabel label;
    double min;
    double max;
    int compX;
    int compY;

    /**
     * @param mySet the ParamSet that can be queried to get the value of this combo box
     * @param label the label of the menuDouble
     * @param initVal the starting setting of the menuDouble
     */
    public DoubleParam(ParamSet mySet, String label, double initVal, double min, double max){
        super(doubleFormat);
        this.setColumns(10);
        this.min=min;
        this.max=max;
        this.labelText=label;
        this.label=new JLabel(labelText);
        this.compX=1;
        this.compY=2;
        mySet.AddGuiMenuItem(this,Double.toString(initVal));
    }

    /**
     * @param mySet the ParamSet that can be queried to get the value of this combo box
     * @param label the label of the menuDouble
     * @param initVal the starting value of the menuDouble
     * @param nCols the number of digits that will fit on the display
     * @param compX the width on the gui GridBagLayout
     * @param compY the height on the gui GridBagLayout
     */
    public DoubleParam(ParamSet mySet, String label, double initVal, double min, double max, int nCols, int compX, int compY){
        super(doubleFormat);
        this.setColumns(nCols);
        this.min=min;
        this.max=max;
        this.labelText=label;
        this.label=new JLabel(labelText);
        this.compX=compX;
        this.compY=compY;
        mySet.AddGuiMenuItem(this,Double.toString(initVal));
    }

    /**
     * sets the foreground and background of the DoubleParam
     * @param foregroundColor color of the text if null the GuiWindow color will be used
     * @param backgroundColor color of the background, if null the GuiWindow color will be used
     */
    public DoubleParam SetColor(Color foregroundColor, Color backgroundColor){
        if(backgroundColor!=null){
            setOpaque(true);
            setBackground(backgroundColor);
            label.setOpaque(true);
            label.setBackground(backgroundColor);
        }
        if(foregroundColor !=null) {
            setForeground(foregroundColor);
            label.setForeground(foregroundColor);
            setCaretColor(foregroundColor);
        }
        return this;
    }

    /**
     * ignore
     */
    @Override
    public int TypeID() {
        return 2;
    }

    /**
     * ignore
     */
    @Override
    public void Set(String val) {
        this.setText(Double.toString(BoundVal(Double.parseDouble(val),min,max)));
    }

    /**
     * ignore
     */
    @Override
    public String Get() {
        String val=Double.toString(BoundVal(Double.parseDouble(this.getText()),min,max));
        this.Set(val);
        return val;
    }

    /**
     * ignore
     */
    @Override
    public String GetLabel() {
        return labelText;
    }

    /**
     * ignore
     */
    @Override
    public int NEntries() {
        return 2;
    }

    /**
     * ignore
     */
    @Override
    public <T extends Component> T GetEntry(int iEntry) {
        switch(iEntry){
            case 0: return (T)label;
            case 1: return (T)this;
            default: throw new IllegalArgumentException(iEntry+" does not match to an item!");
        }
    }


    /**
     * ignore
     */
    @Override
    public int compX() {
        return compX;
    }

    /**
     * ignore
     */
    @Override
    public int compY() {
        return compY;
    }

    @Override
    public boolean IsActive() {
        return true;
    }

    @Override
    public void SetActive(boolean isActive) {

    }

    /**
     * ignore
     */
    @Override
    public void GetComps(ArrayList<Component> putHere, ArrayList<Integer> coordsHere, ArrayList<Integer> compSizesHere){
        int labelEnd=compY/2;
        putHere.add(this.label);
        coordsHere.add(0);
        coordsHere.add(0);
        compSizesHere.add(compX);
        compSizesHere.add(labelEnd);
        putHere.add(this);
        coordsHere.add(0);
        coordsHere.add(labelEnd);
        compSizesHere.add(compX);
        compSizesHere.add(compY-labelEnd);
    }
}
