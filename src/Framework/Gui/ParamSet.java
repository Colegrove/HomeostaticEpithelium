package Framework.Gui;

import Framework.Misc.MenuItem;

import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.util.*;

/*PLAN
Get Base Version Working
Refine Syntax and Methods
*/


/**
 * the ParamSet class is used to keep track of all of the items in a menu
 * it can be queried to get the current values of the menu
 * values and labels can also be loaded as arrays of strings, allowing for the menuset to be used without a gui
 */
public class ParamSet {
    GuiWindow win;
    final int typeID = 1;//integer
    boolean column;
    final ArrayList<MenuItem> items;
    HashMap<String,Integer>labels;
    ArrayList<String>vals;

    public ParamSet() {
        this.items = new ArrayList<MenuItem>();
        this.labels=new HashMap<>();
        this.vals=new ArrayList<>();
    }

    /**
     * @param labelArr array of all label names
     * @param vals array of values, must be in the same order as the labels array
     */
    public ParamSet(String[] labelArr, String[] vals){
        labels=new HashMap<>();
        for(int i=0;i<labelArr.length;i++){
            labels.put(labelArr[i],i);
        }
        this.vals=new ArrayList<String>(Arrays.asList(vals));
        this.items = new ArrayList<MenuItem>();
    }

    /**
     * gets the boolean value from the set with the specified label
     */
    public boolean GetBool(String label){
        String val=vals.get(labels.get(label)).toLowerCase();
        boolean setVal=false;
        if(val.equals("true")||val.equals("t")){
            setVal=true;
        }else if(val.equals("false")||val.equals("f")) {
            setVal = false;
        }
        return setVal;
    }
    /**
     * gets the integer value from the set with the specified label
     */
    public int GetInt(String label){
        return (int)Double.parseDouble(vals.get(labels.get(label)));
    }
    /**
     * gets the double value from the set with the specified label
     */
    public double GetDouble(String label){
        return Double.parseDouble(vals.get(labels.get(label)));
    }
    /**
     * gets the string value from the set with the specified label
     */
    public String GetStr(String label){
        return vals.get(labels.get(label));
    }

    /**
     * sets all values with the vals array provided
     */
    public void SetVals(String[] vals){
        if(vals.length==labels.size()) {
            this.vals = new ArrayList<String>(Arrays.asList(vals));
        }
        else{
            throw new IllegalArgumentException("Values array has wrong size!");
        }
    }

    /**
     * sets the value at the specified label to the specified value
     */
    public void Set(String label,String value){
        vals.set(labels.get(label),value);
    }

    /**
     * sets all labels with the labels array provided
     */
    public void SetLabels(String[] labels){
        for(int i=0;i<labels.length;i++){
            this.labels.put(labels[i],i);
            if(vals.size()<=i){
                vals.add("");
            }
        }
    }

    /**
     * sets both the vals and labels with the arrays provided
     */
    public void SetValsAndLabels(String[]Labels,String[] vals) {
        for(int i=0;i<Labels.length;i++){
            Set(Labels[i],vals[i]);
        }
    }

    /**
     * adds a new gui menu item to the set
     */
    MenuItem AddGuiMenuItem(MenuItem addMe, String initValue) {
        String name = addMe.GetLabel();
        labels.put(addMe.GetLabel(),items.size());
        items.add(addMe);
        addMe.Set(initValue);
        vals.add(initValue);
        //((java.awt.Component)addMe).addPropertyChangeListener(new PropertyChangeListener() {
        ((java.awt.Component)addMe).addFocusListener(new FocusListener() {
            @Override
            public void focusGained(FocusEvent e) {
            }

            @Override
            public void focusLost(FocusEvent e) {
                Set(addMe.GetLabel(),addMe.Get());
            }
        });
        return addMe;
    }

    /**
     * set all gui values with the vals array provided
     */
    public void SetGuiAll(String[] vals) {
        for (int i = 0; i < items.size(); i++) {
            items.get(i).Set(vals[i]);
        }
        SetVals(vals.clone());
    }

    /**
     * sets the gui component with the provided label to val
     */
    public void SetGui(String label, String val){
        items.get(labels.get(label)).Set(val);
        Set(label,val);
    }

    /**
     * returns the values in the set as an array of strings
     */
    public String[] ValueStrings() {
        String[] elems = new String[vals.size()];
        for (int i = 0; i < vals.size(); i++) {
            elems[i] = vals.get(i);
        }
        return elems;
    }

    /**
     * returns the labels in the set as an array of strings
     */
    public String[] LabelStrings() {
        String[] newLabels = new String[labels.size()];
        for (Map.Entry<String,Integer> entry:labels.entrySet()) {
            newLabels[entry.getValue()] = entry.getKey();
        }
        return newLabels;
    }
}
