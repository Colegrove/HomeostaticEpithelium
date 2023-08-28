package Framework.Tools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by bravorr on 10/24/16.
 */
public class FileIO {
    public final String fileName;
    public final char ReadWriteAppend;
    public final boolean binary;
    BufferedReader Reader;
    BufferedWriter Writer;
    DataOutputStream WriterBin;
    DataInputStream ReaderBin;

    /**
     * @param fileName name of the file to read from or write to
     * @param mode should be either "r":read, "w":write, "rb":readBinary, "wb":writeBinary
     */
    public FileIO(String fileName, String mode) {
        char[] modeChars=mode.toCharArray();
        if(modeChars.length>2||(modeChars.length>1&&modeChars[1]!='b')||(modeChars[0]!='w'&&modeChars[0]!='r'&&modeChars[0]!='a')){
            throw new IllegalArgumentException("inccorect mode argument! mode should be 'r' for read, 'w' for write, or 'a' for append, followed by optional 'b' for binary");
        }
        this.fileName=fileName;
        this.ReadWriteAppend=modeChars[0];
        this.binary= modeChars.length > 1;
        boolean appendOut=false;
        if(ReadWriteAppend=='a'){
            appendOut=true;
        }

        try {
            if (ReadWriteAppend == 'r') {
                if (this.binary) {
                    ReaderBin = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));
                } else {
                    Reader = new BufferedReader(new FileReader(fileName));
                }
            } else if (ReadWriteAppend == 'w' || ReadWriteAppend == 'a') {
                if (this.binary) {
                    this.WriterBin = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fileName, appendOut)));
                } else {
                    Writer = new BufferedWriter(new FileWriter(fileName, appendOut));
                }
            }
            else{
                throw new IllegalArgumentException("rwa character must be one of r(read) w(write) or a(append)");
            }
        }catch(IOException e){
            e.printStackTrace();
        }
    }

    //READ FUNCTIONS

    /**
     * requires read mode ("r")
     * pulls a line from the file, splits it by the delimiter, and returns an array of line segments
     * @param delimiter the delimiter used to divide the line
     */
    public String[] ReadLineDelimit(String delimiter){
        String[] read=null;
        try{
            String line=Reader.readLine();
            if(line!=null){
               read=line.split(delimiter);
            }
        }
        catch(IOException e){
            e.printStackTrace();
        }
        return read;
    }

    /**
     * requires read mode ("r")
     * pulls a line from the file, splits it by the delimeter and returns an array of integers from that line
     * @param delimiter the delimiter used to divide the line
     */
    public int[]ReadLineIntDelimit(String delimiter){
        String[] raw=ReadLineDelimit(delimiter);
        if(raw==null) {return null;}
        int[] ret=new int[raw.length];
        for (int i = 0; i < raw.length; i++) {
            ret[i] = Integer.parseInt(raw[i]);
        }
        return ret;
    }
    /**
     * requires read mode ("r")
     * pulls a line from the file, splits it by the delimeter and returns an array of floats from that line
     * @param delimiter the delimiter used to divide the line
     */
    public float[]ReadLineFloatDelimit(String delimiter){
        String[] raw=ReadLineDelimit(delimiter);
        float[] ret=new float[raw.length];
        for(int i=0;i<ret.length;i++){
            ret[i]=Float.parseFloat(raw[i]);
        }
        return ret;
    }
    /**
     * requires read mode ("r")
     * pulls a line from the file, splits it by the delimeter and returns an array of doubles from that line
     * @param delimiter the delimiter used to divide the line
     */
    public double[]ReadLineDoubleDelimit(String delimiter){
        String[] raw=ReadLineDelimit(delimiter);
        if(raw == null) {return null;}
        double[] ret=new double[raw.length];
        for (int i = 0; i < raw.length; i++) {
            ret[i] = Double.parseDouble(raw[i]);
        }
        return ret;
    }
    /**
     * requires read mode ("r")
     * pulls a line from the file, splits it by the delimeter and returns an array of longs from that line
     * @param delimiter the delimiter used to divide the line
     */
    public long[]ReadLineLongDelimit(String delimiter){
        String[] raw=ReadLineDelimit(delimiter);
        if(raw==null){return null;}
        long[] ret=new long[raw.length];
        for(int i=0;i<raw.length;i++){
            ret[i]=Long.parseLong(raw[i]);
        }
        return ret;
    }

    /**
     * requires read mode ("r")
     * pulls all lines from the file, splits them by the delimiter, and returns an arraylist, where each entry is an array of line segments from that line
     * @param delimiter the delimiter used to divide the lines
     */
    public ArrayList<String[]> ReadDelimit(String delimiter){
        ArrayList<String[]>lines=new ArrayList<>();
        String[] read=ReadLineDelimit(delimiter);
        while(read!=null){
            lines.add(read);
            read=ReadLineDelimit(delimiter);
        }
        return lines;
    }
    /**
     * requires read mode ("r")
     * pulls all lines from the file, splits them by the delimiter, and returns an arraylist, where each entry is an array of doubles from that line
     * @param delimiter the delimiter used to divide the lines
     */
    public ArrayList<double[]> ReadDoubleDelimit(String delimiter){
        ArrayList<double[]>lines=new ArrayList<>();
        double[] read=ReadLineDoubleDelimit(delimiter);
        while(read!=null){
            lines.add(read);
            read=ReadLineDoubleDelimit(delimiter);
        }
        return lines;
    }
    /**
     * requires read mode ("r")
     * pulls all lines from the file, splits them by the delimiter, and returns an arraylist, where each entry is an array of ints from that line
     * @param delimiter the delimiter used to divide the lines
     */
    public ArrayList<int[]> ReadIntDelimit(String delimiter){
        ArrayList<int[]>lines=new ArrayList<>();
        int[] read=ReadLineIntDelimit(delimiter);
        while(read!=null){
            lines.add(read);
            read=ReadLineIntDelimit(delimiter);
        }
        return lines;
    }
    /**
     * requires read mode ("r")
     * pulls all lines from the file, splits them by the delimiter, and returns an arraylist, where each entry is an array of longs from that line
     * @param delimiter the delimiter used to divide the lines
     */
    public ArrayList<long[]> ReadLongDelimit(String delimiter){
        ArrayList<long[]>lines=new ArrayList<>();
        long[] read=ReadLineLongDelimit(delimiter);
        while(read!=null){
            lines.add(read);
            read=ReadLineLongDelimit(delimiter);
        }
        return lines;
    }

    /**
     * requires read mode ("r")
     * returns one line from the file as a string
     */
    public String ReadLine(){
        String read="";
        try {
            read=Reader.readLine();
        }catch(IOException e){
            e.printStackTrace();
        }
        return read;
    }
    /**
     * requires read mode ("r")
     * returns an array list of all lines from the file as strings
     */
    public ArrayList<String> Read() {
        ArrayList<String>ret=new ArrayList<>();
        String line=ReadLine();
        while(line!=null){
            ret.add(line);
            line=ReadLine();
        }
        return ret;
    }
    //WRITE FUNCTIONS

    /**
     * requires write mode ("w")
     * @param line writes the line to the file
     */
    public void Write(String line) {
        try {
            this.Writer.write(line);
        }catch(IOException e){
            e.printStackTrace();
        }
    }
    /**
     * requires write mode ("w")
     * the objects in the list are written to the output file using the toString() method
     * @param data a list of objects to be written
     * @param delimiter the delimiter used to separate each object
     */
    public void WriteDelimit(List<Object> data, String delimiter){
        if(data.size()>0) {
            for (int i = 0; i < data.size() - 1; i++) {
                Write(data.get(i).toString() + delimiter);
            }
            Write(data.get(data.size() - 1).toString());
        }
    }

    public void WriteStrings(List<String> data, String delimiter){
        if(data.size()>0) {
            for (int i = 0; i < data.size() - 1; i++) {
                Write(data.get(i) + delimiter);
            }
            Write(data.get(data.size() - 1));
        }
    }
    /**
     * requires write mode ("w")
     * the objects in the array are written to the output file using the toString() method
     * @param data an array of objects to be written
     * @param delimiter the delimiter used to separate each object
     */
    public void WriteDelimit(Object[] data,String delimiter){
        if(data.length>0) {
            for (int i = 0; i < data.length - 1; i++) {
                Write(data[i].toString() + delimiter);
            }
            Write(data[data.length - 1].toString());
        }
    }

    /**
     * requires write mode ("w")
     * writes the array of ints to the file, separated by the delimiter
     */
    public void WriteDelimit(int[] data,String delimiter){
        if(data.length>0) {
            for (int i = 0; i < data.length - 1; i++) {
                Write(String.valueOf(data[i]) + delimiter);
            }
            Write(String.valueOf(data[data.length - 1]));
        }
    }
    /**
     * requires write mode ("w")
     * writes the array of longs to the file, separated by the delimiter
     */
    public void WriteDelimit(long[] data,String delimiter){
        if(data.length>0) {
            for (int i = 0; i < data.length - 1; i++) {
                Write(String.valueOf(data[i]) + delimiter);
            }
            Write(String.valueOf(data[data.length - 1]));
        }
    }
    /**
     * requires write mode ("w")
     * writes the array of floats to the file, separated by the delimiter
     */
    public void WriteDelimit(float[] data,String delimiter){
        if(data.length>0) {
            for (int i = 0; i < data.length - 1; i++) {
                Write(String.valueOf(data[i]) + delimiter);
            }
            Write(String.valueOf(data[data.length - 1]));
        }
    }
    /**
     * requires write mode ("w")
     * writes the array of doubles to the file, separated by the delimiter
     */
    public void WriteDelimit(double[] data,String delimiter){
        if(data.length>0) {
            for (int i = 0; i < data.length - 1; i++) {
                Write(String.valueOf(data[i]) + delimiter);
            }
            Write(String.valueOf(data[data.length - 1]));
        }
    }

    //BIN WRITE FUNCTIONS

    /**
     * requires writeBinary model ("wb")
     * writes a single double to the binary file
     */
    public void WriteBinDouble(double writeMe) {
        try{
            WriterBin.writeDouble(writeMe);
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    /**
     * requires writeBinary model ("wb")
     * writes a single float to the binary file
     */
    public void WriteBinFloat(float writeMe) {
        try{
        WriterBin.writeFloat(writeMe);
    }
        catch(IOException e){
        e.printStackTrace();
    }
    }
    /**
     * requires writeBinary model ("wb")
     * writes a single int to the binary file
     */
    public void WriteBinInt(int writeMe) {
        try{
        WriterBin.writeInt(writeMe);
    }
        catch(IOException e){
                e.printStackTrace();
                }
    }
    /**
     * requires writeBinary model ("wb")
     * writes a single long to the binary file
     */
    public void WriteBinLong(long writeMe) {
        try{
        WriterBin.writeLong(writeMe);
        }
        catch(IOException e){
        e.printStackTrace();
        }
    }
    /**
     * requires writeBinary model ("wb")
     * writes a string to the binary file
     */
    public void WriteBinString(String writeMe){
        try{
        WriterBin.writeChars(writeMe);
        }
        catch(IOException e){
        e.printStackTrace();
        }
    }
    /**
     * requires writeBinary model ("wb")
     * writes a single boolean to the binary file
     */
    public void WriteBinBoolean(boolean writeMe){
        try{
        WriterBin.writeBoolean(writeMe);
        }
        catch(IOException e){
        e.printStackTrace();
        }
    }

    //BIN READ FUNCTIONS

    /**
     * requires readBinary mode ("rb")
     * reads a single double from the binary file
     */
    public double ReadBinDouble() {
        double ret=0;
        try{
        ret=ReaderBin.readDouble();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }
    /**
     * requires readBinary mode ("rb")
     * reads a single float from the binary file
     */
    public float ReadBinFloat() {
        float ret=0;
        try{
        ret=ReaderBin.readFloat();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }
    /**
     * requires readBinary mode ("rb")
     * reads a single int from the binary file
     */
    public int ReadBinInt() {
        int ret=0;
        try{
        ret=ReaderBin.readInt();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }
    /**
     * requires readBinary mode ("rb")
     * reads a single long from the binary file
     */
    public long ReadBinLong() {
        long ret=0;
        try{
        ret=ReaderBin.readLong();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }
    /**
     * requires readBinary mode ("rb")
     * reads a string from the binary file
     */
    public String ReadBinString(){
        String ret="";
        try{
        ret=ReaderBin.readUTF();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }
    /**
     * requires readBinary mode ("rb")
     * reads a boolean from the binary file
     */
    public boolean ReadBinBoolean(){
        boolean ret=false;
        try{
        ret=ReaderBin.readBoolean();
        }
        catch(IOException e){
        e.printStackTrace();
        }
        return ret;
    }

    /**
     * returns the length of the file, or 0 if the file does not exist
     */
    public double length(){
        File file=new File(fileName);
        return file.exists()?file.length():0;
    }

    //CLOSE

    /**
     * used to close the file. if forgotten, write calls may never finish
     */
    public void Close() {
        try {
            if(ReadWriteAppend=='r') {
                if (binary) {
                    ReaderBin.close();
                    ReaderBin=null;
                } else {
                    Reader.close();
                    Reader=null;
                }
            }
            else if (ReadWriteAppend == 'w'||ReadWriteAppend=='a') {
                if (binary) {
                    WriterBin.close();
                    WriterBin=null;
                } else {
                    Writer.close();
                    Writer=null;
                }
            }
        }catch(IOException e){
            e.printStackTrace();
        }
    }

}
