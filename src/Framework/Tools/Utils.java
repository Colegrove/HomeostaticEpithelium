package Framework.Tools;

import Framework.Misc.SweepRunFunction;
import Framework.Misc.SweepRun;

import java.awt.*;
import java.lang.reflect.Array;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 *
 * A collection of helpful static utility functions
 * recommended import: import static Utils.*
 * Created by rafael on 10/11/16.
 */
public final class Utils {

    static double DOUBLE_EPSILON=2.22E-16;

    public static int ColorInt(double r,double g,double b){
        int ri=(int)(r*255);
        int gi=(int)(g*255);
        int bi=(int)(b*255);
        return (((ri&0x0ff)<<16)|((gi&0x0ff)<<8)|(bi&0x0ff))-16777216;

    }
    public static void HeatMapping(double val, float[]out){
        out[0]=(float) Math.min(1, val * 4);
        out[1]=val > 0.25?(float) Math.min(1, (val - 0.25) * 2):0;
        out[2]=val > 0.75?(float) Math.min(1, (val - 0.75) * 4):0;
    }
    public static void HeatMapping(double val, double[]out){
        out[0]=(float) Math.min(1, val * 4);
        out[1]=val > 0.25? Math.min(1, (val - 0.25) * 2):0;
        out[2]=val > 0.75? Math.min(1, (val - 0.75) * 4):0;
    }
    public static void HSBtoRGB(double hue, double saturation, double brightness,double[] retRGB) {
        float scale = 255f;
        int rgb = Color.HSBtoRGB((float)hue, (float)saturation, (float)brightness);
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >> 8) & 0xFF;
        int b = rgb & 0xFF;
        retRGB[0]=r/scale;
        retRGB[1]=g/scale;
        retRGB[2]=b/scale;
    }
    //Added from github version Darryl3DVersion 25Aug23HLC
    public static float[] GetHSBtoRGB(float hue, float saturation, float brightness) {
        float scale = 255f;
        int rgb = Color.HSBtoRGB(hue, saturation, brightness);
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >> 8) & 0xFF;
        int b = rgb & 0xFF;
        float[] RGBreturn = {r/scale,g/scale,b/scale};
        return RGBreturn;
    }
    public static double SumArray(double[] arr){
        double sum=0;
        for (double val : arr) {
            sum+=val;
        }
        return sum;
    }
    public static int SumArray(int[] arr){
        int sum=0;
        for (int val : arr) {
            sum+=val;
        }
        return sum;
    }
    /**
     * prints an array to a string
     *
     * @param arr   array to be printed
     * @param delim the delimiter used to separate entries
     * @param <T>   the type of the data entries in the array
     */
    public static <T> String ArrToString(T[] arr, String delim) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i] + delim);
        }
        return sb.toString();
    }

    /**
     * prints an array
     *
     * @param arr   array to be printed
     * @param delim the delimiter used to separate entries
     */
    public static String ArrToString(double[] arr, String delim) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i] + delim);
        }
        return sb.toString();
    }
    public static String ArrToString(float[] arr, String delim) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i] + delim);
        }
        return sb.toString();
    }

    /**
     * prints an array
     *
     * @param arr   array to be printed
     * @param delim the delimiter used to separate entries
     */
    public static String ArrToString(int[] arr, String delim) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i] + delim);
        }
        return sb.toString();
    }
    public static String ArrToString(long[] arr, String delim) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i] + delim);
        }
        return sb.toString();
    }

    /**
     * Returns a new array that is the first array with the second concatenated to the end of it
     *
     * @param <T> the type of the input and output arrays
     */
    public static <T> T[] Concat(T[] first, T[] second) {
        int firstLen = first.length;
        int secondLen = second.length;
        T[] ret = (T[]) Array.newInstance(first.getClass().getComponentType(), firstLen + secondLen);
        System.arraycopy(first, 0, ret, 0, firstLen);
        System.arraycopy(second, 0, ret, firstLen, firstLen + secondLen);
        return ret;
    }

    /**
     * Returns a new array that is the first array with the appendMe object appended to the end of it
     *
     * @param <T> the type of the inputs and output array
     */
    public static <T> T[] Append(T[] arr, T appendMe) {
        int firstLen = arr.length;
        T[] ret = (T[]) Array.newInstance(arr.getClass().getComponentType(), firstLen + 1);
        System.arraycopy(arr, 0, ret, 0, firstLen);
        ret[firstLen] = appendMe;
        return ret;
    }

    public static int[] IndicesArray(int nEntries){
        int indices[] = new int[nEntries];
        for (int i = 0; i < nEntries; i++) {
            indices[i] = i;
        }
        return indices;
    }
    /**
     * Generates a list of sequential integers starting at 0, and then shuffles them
     *
     * @param nEntries    the length of the list to be generated
     * @param CountRandom the number of elements that should be shuffled
     * @param rn          the random number generator to be used
     * @return returns the array of indices after being shuffled
     */

    public static int[] RandomIndices(int nEntries, int CountRandom, Random rn) {
        int indices[] = new int[nEntries];
        for (int i = 0; i < nEntries; i++) {
            indices[i] = i;
        }
        Shuffle(indices, indices.length, CountRandom, rn);
        return indices;
    }

    /**
     * Fills out with random doubles between min and max inclusive
     *
     * @param out the array the random doubles should be written to. the length of the input array defines the number of doubles to be generated
     * @param rn  the random number generator to be used
     */
    public static void RandomDS(double[] out, double min, double max, Random rn) {
        for (int i = 0; i < out.length; i++) {
            out[i] = rn.nextDouble() * (max - min) + min;
        }
    }
    /**
     * Fills out with random integers between min (inclusive) and max (exclusive)
     *
     * @param out the array the random doubles should be written to. the length of the input array defines the number of doubles to be generated
     * @param rn  the random number generator to be used
     */
    public static void RandomIS(int[] out, int min, int max, Random rn) {
        for (int i = 0; i < out.length; i++) {
            out[i] = rn.nextInt(max - min) + min;
        }
    }


    /**
     * returns the mean value of the provided array
     */
    static public double Mean(double[] a) {
        double tot = 0;
        for (int i = 0; i < a.length; i++) {
            tot += a[i];
        }
        return tot / a.length;
    }



    /**
     * Samples a gaussian with the provided mean and standard deviation
     *
     * @param mean   the mean of the gaussian
     * @param stdDev the standard deviation of the gaussian
     * @param rn     the random number generator that will be used
     * @return a single value sampled from the defined Gaussian
     */
    public static double Gaussian(double mean, double stdDev, Random rn) {
        return rn.nextGaussian() * stdDev + mean;
    }

    /**
     * Returns the coordinates defining the Von Neumann neighborhood centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim...]
     */
    public static int[] VonNeumannHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 1, 0, -1, 0, 0, 1, 0, -1};
        } else {
            return new int[]{1, 0, -1, 0, 0, 1, 0, -1};
        }
    }

    /**
     * Returns the coordinates defining the Moore neighborhood centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] MooreHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 1, 1, 1, 0, 1, -1, 0, -1, -1, -1, -1, 0, -1, 1, 0, 1};
        } else {
            return new int[]{1, 1, 1, 0, 1, -1, 0, -1, -1, -1, -1, 0, -1, 1, 0, 1};
        }
    }

    /**
     * Returns the coordinates defining the Hexagonal neighborhood for even yDim coordinates centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] HexEvenHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 1, 1, 1, 0, 1, -1, 0, -1, -1, 0, 0, 1};
        } else {
            return new int[]{1, 1, 1, 0, 1, -1, 0, -1, -1, 0, 0, 1};
        }
    }

    /**
     * Returns the coordinates defining the Hexagonal neighborhood for odd yDim coordinates centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] HexOddHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 1, 0, 0, -1, -1, -1, -1, 0, -1, 1, 0, 1};
        } else {
            return new int[]{1, 0, 0, -1, -1, -1, -1, 0, -1, 1, 0, 1};
        }
    }

    /**
     * Returns the coordinates defining the Triangular neighborhood for even xDim, even yDim or oddx, odd yDim. centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] TriangleSameParityHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, -1, 0, 1, 0, 0, 1};
        } else {
            return new int[]{-1, 0, 1, 0, 0, 1};
        }
    }

    /**
     * Returns the coordinates defining the Triangular neighborhood for even xDim, odd yDim or oddx, even yDim. centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] TriangleDiffParityHood(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, -1, 0, 1, 0, 0, -1};
        } else {
            return new int[]{-1, 0, 1, 0, 0, -1};
        }
    }


    /**
     * Returns the coordinates defining the Von Neumann neighborhood centered on (0,0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0,0)
     * @return coordinates returned as an array of the form [xDim,yDim,z,xDim,yDim,z,...]
     */
    public static int[] VonNeumannHood3D(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
        } else {
            return new int[]{1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
        }
    }

    public int[] MooreHood3d(boolean includeOrigin) {
        if (includeOrigin) {
            return new int[]{0, 0, 0,
                    0, 0, 1,
                    0, 0, -1,
                    1, 0, 0,
                    1, 0, 1,
                    1, 0, -1,
                    1, 1, 0,
                    1, 1, 1,
                    1, 1, -1,
                    0, 1, 0,
                    0, 1, 1,
                    0, 1, -1,
                    -1, 0, 0,
                    -1, 0, 1,
                    -1, 0, -1,
                    -1, 1, 0,
                    -1, 1, 1,
                    -1, 1, -1,
                    -1, -1, 0,
                    -1, -1, 1,
                    -1, -1, -1,
                    0, -1, 0,
                    0, -1, 1,
                    0, -1, -1,
                    1, -1, 0,
                    1, -1, 1,
                    1, -1, -1,
            };
        } else {
            return new int[]{
                    0, 0, 1,
                    0, 0, -1,
                    1, 0, 0,
                    1, 0, 1,
                    1, 0, -1,
                    1, 1, 0,
                    1, 1, 1,
                    1, 1, -1,
                    0, 1, 0,
                    0, 1, 1,
                    0, 1, -1,
                    -1, 0, 0,
                    -1, 0, 1,
                    -1, 0, -1,
                    -1, 1, 0,
                    -1, 1, 1,
                    -1, 1, -1,
                    -1, -1, 0,
                    -1, -1, 1,
                    -1, -1, -1,
                    0, -1, 0,
                    0, -1, 1,
                    0, -1, -1,
                    1, -1, 0,
                    1, -1, 1,
                    1, -1, -1,
            };
        }
    }

    /**
     * Returns an array of all squares touching a line between the positions provided
     *
     * @param x0 the xDim coordinate of the starting position
     * @param y0 the yDim coordinate of the starting position
     * @param x1 the xDim coordinate of the ending position
     * @param y1 the yDim coordinate of the ending position
     * @return coordinates return as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int[] SquaresAlongLine(double x0, double y0, double x1, double y1) {
        double dx = Math.abs(x1 - x0);
        double dy = Math.abs(y1 - y0);

        int x = (int) (Math.floor(x0));
        int y = (int) (Math.floor(y0));

        int n = 1;
        int x_inc, y_inc;
        double error;

        if (dx == 0) {
            x_inc = 0;
            error = Double.MAX_VALUE;
        } else if (x1 > x0) {
            x_inc = 1;
            n += (int) (Math.floor(x1)) - x;
            error = (Math.floor(x0) + 1 - x0) * dy;
        } else {
            x_inc = -1;
            n += x - (int) (Math.floor(x1));
            error = (x0 - Math.floor(x0)) * dy;
        }

        if (dy == 0) {
            y_inc = 0;
            error -= Double.MAX_VALUE;
        } else if (y1 > y0) {
            y_inc = 1;
            n += (int) (Math.floor(y1)) - y;
            error -= (Math.floor(y0) + 1 - y0) * dx;
        } else {
            y_inc = -1;
            n += y - (int) (Math.floor(y1));
            error -= (y0 - Math.floor(y0)) * dx;
        }

        int[] writeHere = new int[n ^ 2];
        int Count = 0;
        for (; n > 0; --n) {
            writeHere[Count * 2] = (int) Math.floor(x);
            writeHere[Count * 2] = (int) Math.floor(y);
            Count++;

            if (error > 0) {
                y += y_inc;
                error -= dx;
            } else {
                x += x_inc;
                error += dy;
            }
        }
        return writeHere;
    }
    /**
     * Returns an array of all squares touching a line between the positions provided
     *
     * @param x0 the xDim coordinate of the starting position
     * @param y0 the yDim coordinate of the starting position
     * @param x1 the xDim coordinate of the ending position
     * @param y1 the yDim coordinate of the ending position
     * @return coordinates return as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    public static int SquaresAlongLine(double x0, double y0, double x1, double y1,int[] returnCoords) {
        double dx = Math.abs(x1 - x0);
        double dy = Math.abs(y1 - y0);

        int x = (int) (Math.floor(x0));
        int y = (int) (Math.floor(y0));

        int n = 1;
        int x_inc, y_inc;
        double error;

        if (dx == 0) {
            x_inc = 0;
            error = Double.MAX_VALUE;
        } else if (x1 > x0) {
            x_inc = 1;
            n += (int) (Math.floor(x1)) - x;
            error = (Math.floor(x0) + 1 - x0) * dy;
        } else {
            x_inc = -1;
            n += x - (int) (Math.floor(x1));
            error = (x0 - Math.floor(x0)) * dy;
        }

        if (dy == 0) {
            y_inc = 0;
            error -= Double.MAX_VALUE;
        } else if (y1 > y0) {
            y_inc = 1;
            n += (int) (Math.floor(y1)) - y;
            error -= (Math.floor(y0) + 1 - y0) * dx;
        } else {
            y_inc = -1;
            n += y - (int) (Math.floor(y1));
            error -= (y0 - Math.floor(y0)) * dx;
        }

        int Count = 0;
        for (; n > 0; --n) {
            returnCoords[Count * 2] = (int) Math.floor(x);
            returnCoords[Count * 2] = (int) Math.floor(y);
            Count++;

            if (error > 0) {
                y += y_inc;
                error -= dx;
            } else {
                x += x_inc;
                error += dy;
            }
        }
        return Count;
    }

    /**
     * Returns the coordinates of all squares whose centers lie within a circle of the provided radius, centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @param radius        the radius of the circle
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    static public int[] CircleCentered(boolean includeOrigin, double radius) {
        double distSq = radius * radius;
        int min = (int) Math.floor(-radius);
        int max = (int) Math.ceil(radius);
        int[] retLong = new int[((max + 1 - min) * (max + 1 - min)) * 2];
        int ct = 0;
        if (includeOrigin) {
            ct++;
            retLong[0] = 0;
            retLong[1] = 0;
        }
        for (int x = min; x <= max; x++) {
            for (int y = min; y <= max; y++) {
                if (Utils.DistSq2D(0, 0, x, y) <= distSq) {
                    if (x == 0 && y == 0) {
                        continue;
                    }
                    retLong[ct * 2] = x;
                    retLong[ct * 2 + 1] = y;
                    ct++;
                }
            }
        }
        int[] ret = new int[ct * 2];
        System.arraycopy(retLong, 0, ret, 0, ret.length);
        return ret;
    }

    static public int CircleCentered(boolean includeOrigin, double radius,int[] returnCoords) {
        double distSq = radius * radius;
        int min = (int) Math.floor(-radius);
        int max = (int) Math.ceil(radius);
        int ct = 0;
        if (includeOrigin) {
            ct++;
            returnCoords[0] = 0;
            returnCoords[1] = 0;
        }
        for (int x = min; x <= max; x++) {
            for (int y = min; y <= max; y++) {
                if (Utils.DistSq2D(0, 0, x, y) <= distSq) {
                    if (x == 0 && y == 0) {
                        continue;
                    }
                    returnCoords[ct * 2] = x;
                    returnCoords[ct * 2 + 1] = y;
                    ct++;
                }
            }
        }
        return ct;
    }
    /**
     * Returns the coordinates of all squares whose centers lie within a rectangle of the provided radius, centered on (0,0)
     *
     * @param includeOrigin defines whether to include the origin (0,0)
     * @param radX          the radius of the rectangle in the xDim direction
     * @param radY          the radius of the rectangle in the yDim direction
     * @return coordinates returned as an array of the form [xDim,yDim,xDim,yDim,...]
     */
    static public int[] RectCentered(boolean includeOrigin, int radX, int radY) {
        //returns a square with a center location at 0,0
        int[] dataIn;
        int nCoord;
        if (includeOrigin) {
            dataIn = new int[(radX * 2 + 1) * (radY * 2 + 1) * 2];
            dataIn[0] = 0;
            dataIn[1] = 0;
            nCoord = 1;
        } else {
            dataIn = new int[(radX * 2 + 1) * (radY * 2 + 1) * 2 - 1];
            nCoord = 0;
        }
        for (int x = -radX; x <= radX; x++) {
            for (int y = -radY; y <= radY; y++) {
                if (x == 0 && y == 0) {
                    continue;
                }
                dataIn[nCoord * 2] = x;
                dataIn[nCoord * 2 + 1] = y;
                nCoord++;
            }
        }
        return dataIn;
    }
    static public int RectCentered(boolean includeOrigin, int radX, int radY,int[]returnCoords) {
        //returns a square with a center location at 0,0
        int nCoord;
        if (includeOrigin) {
            returnCoords[0] = 0;
            returnCoords[1] = 0;
            nCoord = 1;
        } else {
            returnCoords = new int[(radX * 2 + 1) * (radY * 2 + 1) * 2 - 1];
            nCoord = 0;
        }
        for (int x = -radX; x <= radX; x++) {
            for (int y = -radY; y <= radY; y++) {
                if (x == 0 && y == 0) {
                    continue;
                }
                returnCoords[nCoord * 2] = x;
                returnCoords[nCoord * 2 + 1] = y;
                nCoord++;
            }
        }
        return nCoord;
    }


    /**
     * an efficient implementation of the multinomial function
     */
    public static void Multinomial(double[] probabilities,int n,Binomial bn,int[] ret){
        double pSum=1;
        for (int i = 0; i < probabilities.length; i++) {
            int ni=bn.generateBinomial(n,probabilities[i]/pSum);
            ret[i]=ni;
            n-=ni;
            pSum-=probabilities[i];
        }
        if(Math.abs(pSum-1)>DOUBLE_EPSILON){
            throw new IllegalArgumentException("Multinomial probabilities array must sum to 1");
        }
    }
    public static void Multinomial(double[] probabilities,int n,Binomial bn,int[] ret,int nProbs){
        double pSum=1;
        for (int i = 0; i < nProbs; i++) {
            int ni=bn.generateBinomial(n,probabilities[i]/pSum);
            ret[i]=ni;
            n-=ni;
            pSum-=probabilities[i];
        }
        if(Math.abs(pSum-1)>DOUBLE_EPSILON){
            throw new IllegalArgumentException("Multinomial probabilities array must sum to 1");
        }
    }
    public static void Multinomial(double[] probabilities,long n,Binomial bn,long[] ret){
        double pSum=1;
        for (int i = 0; i < probabilities.length; i++) {
            long ni=bn.generateBinomial(n,probabilities[i]/pSum);
            ret[i]=ni;
            n-=ni;
            pSum-=probabilities[i];
        }
        if(Math.abs(pSum-1)>DOUBLE_EPSILON){
            throw new IllegalArgumentException("Multinomial probabilities array must sum to 1");
        }
    }
    public static void Multinomial(double[] probabilities,long n,Binomial bn,long[] ret,int nProbs){
        double pSum=1;
        for (int i = 0; i < nProbs; i++) {
            long ni=bn.generateBinomial(n,probabilities[i]/pSum);
            ret[i]=ni;
            n-=ni;
            pSum-=probabilities[i];
        }
        if(Math.abs(pSum-1)>DOUBLE_EPSILON){
            throw new IllegalArgumentException("Multinomial probabilities array must sum to 1");
        }
    }
    /**
     * returns the number of heads from nTrials coin flips, where successProb is the probability of heads
     */
    public static int BinomialNaive(double successProb, int nTrials, Random rn) {
        int ret = 0;
        for (int iTrial = 0; iTrial < nTrials; iTrial++) {
            if (rn.nextDouble() < successProb) {
                ret++;
            }
        }
        return ret;
    }

    /**
     *gets a random point on the surface of a sphere centered at 0,0,0, with the provided radius. the x,y,z coords are put in the double[] ret
     */
    public static void RandomPointOnSphereEdge(double radius, Random rn, double[] ret){
        double x=Gaussian(0,radius,rn);
        double y=Gaussian(0,radius,rn);
        double z=Gaussian(0,radius,rn);
        double norm=Norm(x,y,z);
        ret[0]=(x*radius)/norm;
        ret[1]=(y*radius)/norm;
        ret[2]=(z*radius)/norm;
    }

    /**
     *gets a random point on the surface of a circle centered at 0,0, with the provided radius. the x,y coords are put in the double[] ret
     */
    public static void RandomPointOnCircleEdge(double radius, Random rn, double[] ret){
        double x=Gaussian(0,radius,rn);
        double y=Gaussian(0,radius,rn);
        double norm=Norm(x,y);
        ret[0]=(x*radius)/norm;
        ret[1]=(y*radius)/norm;
    }

    public static void RandomPointInCircle(double radius,Random rn,double[] ret){
        double r=Math.sqrt(rn.nextDouble())*radius;
        double a=rn.nextDouble()*Math.PI*2;
        ret[0]=r*Math.cos(a);
        ret[1]=r*Math.sin(a);
    }



    //OTHER COORDINATE FUNCTIONS

    /**
     * Writes to the array argument all squares touching a line between the positions provided
     *
     * @param writeHere the array to write the coordinates into
     * @param x0        the xDim coordinate of the starting position
     * @param y0        the yDim coordinate of the starting position
     * @param x1        the xDim coordinate of the ending position
     * @param y1        the yDim coordinate of the ending position
     * @return the number of touching squares
     */
    public static int SquaresAlongLine(int[] writeHere, double x0, double y0, double x1, double y1) {
        double dx = Math.abs(x1 - x0);
        double dy = Math.abs(y1 - y0);

        int x = (int) (Math.floor(x0));
        int y = (int) (Math.floor(y0));

        int n = 1;
        int x_inc, y_inc;
        double error;

        if (dx == 0) {
            x_inc = 0;
            error = Double.MAX_VALUE;
        } else if (x1 > x0) {
            x_inc = 1;
            n += (int) (Math.floor(x1)) - x;
            error = (Math.floor(x0) + 1 - x0) * dy;
        } else {
            x_inc = -1;
            n += x - (int) (Math.floor(x1));
            error = (x0 - Math.floor(x0)) * dy;
        }

        if (dy == 0) {
            y_inc = 0;
            error -= Double.MAX_VALUE;
        } else if (y1 > y0) {
            y_inc = 1;
            n += (int) (Math.floor(y1)) - y;
            error -= (Math.floor(y0) + 1 - y0) * dx;
        } else {
            y_inc = -1;
            n += y - (int) (Math.floor(y1));
            error -= (y0 - Math.floor(y0)) * dx;
        }

        int Count = 0;
        for (; n > 0; --n) {
            writeHere[Count * 2] = (int) Math.floor(x);
            writeHere[Count * 2] = (int) Math.floor(y);
            Count++;

            if (error > 0) {
                y += y_inc;
                error -= dx;
            } else {
                x += x_inc;
                error += dy;
            }
        }
        return Count;
    }


    //MATH FUNCTIONS

    /**
     * Samples a discrete random variable from the probabilities provided
     *
     * @param probs an array of probabilities. should sum to 1
     * @param rn    the random number generator to be used
     * @return the index of the probability bin that was sampled
     */
    public static int RandomVariable(double[] probs, Random rn) {
        double rand = rn.nextDouble();
        for (int i = 0; i < probs.length; i++) {
            rand -= probs[i];
            if (rand <= 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * samples a random variable nSamples times, results are put in ret
     */
    public static void RandomVariableSample(double[] probs, Random rn, int[] ret, int nSamples) {
        for (int i = 0; i < nSamples; i++) {
            ret[i] = RandomVariable(probs, rn);
        }
    }

    /**
     * sets the values in the array such that they sum to 1
     *
     * @param vals an array of values
     */
    public static void SumTo1(double[] vals) {
        double tot = 0;
        for (int i = 0; i < vals.length; i++) {
            tot += vals[i];
        }
        for (int i = 0; i < vals.length; i++) {
            vals[i] = vals[i] / tot;
        }
    }

    /**
     * uses the Michaelis Menten equation to compute the reaction rate for a given substrate concentration
     *
     * @param conc         concentration of the reaction limiting substrate
     * @param maxRate      reaction rate given maximum concentration
     * @param halfRateConc substrate concentration at which the reaction rate is 1/2 the maximum
     * @return the reaction rate at the given substrate concentration
     */
    public static double MichaelisMenten(final double conc, final double maxRate, final double halfRateConc) {
        if (conc > 0) {
            return (maxRate * conc) / (halfRateConc + conc);
        }
        return 0;
    }

    public static double HillEqn(final double conc, final double dissociationRate, final double hillCoef) {
        if (conc > 0) {
            double adj_conc = Math.pow(conc, hillCoef);
            return adj_conc/(dissociationRate + adj_conc);
        }
        return 0;
    }

    /**
     * Shuffles an array of integers
     *
     * @param arr          array to be shuffled
     * @param lenToShuffle number of elements from array that shuffling can nextField
     * @param Count        number of elements that will be shuffled, should not exceed lenToShuffle
     * @param rn           the random number generator to be used
     */
    public static void Shuffle(int[] arr, int lenToShuffle, int Count, Random rn) {
        for (int i = 0; i < Count; i++) {
            int iSwap = rn.nextInt(lenToShuffle - i) + i;
            int swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }

    public static void Shuffle(int[] arr, Random rn) {
        for (int i = 0; i < arr.length; i++) {
            int iSwap = rn.nextInt(arr.length - i) + i;
            int swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }
    /**
     * Shuffles an array of doubles
     *
     * @param arr          array to be shuffled
     * @param lenToShuffle number of elements from array that shuffling can nextField
     * @param Count        number of elements that will be shuffled, should not exceed lenToShuffle
     * @param rn           the random number generator to be used
     */
    public static void Shuffle(double[] arr, int lenToShuffle, int Count, Random rn) {
        for (int i = 0; i < Count; i++) {
            int iSwap = rn.nextInt(lenToShuffle - i) + i;
            double swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }

    public static void Shuffle(double[] arr, Random rn) {
        for (int i = 0; i < arr.length; i++) {
            int iSwap = rn.nextInt(arr.length - i) + i;
            double swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }
    /**
     * Shuffles an array of objects
     *
     * @param arr          array to be shuffled
     * @param lenToShuffle number of elements from array that shuffling can nextField
     * @param Count        number of elements that will be shuffled, should not exceed lenToShuffle
     * @param rn           the random number generator to be used
     */
    public static void Shuffle(Object[] arr, int lenToShuffle, int Count, Random rn) {
        for (int i = 0; i < Count; i++) {
            int iSwap = rn.nextInt(lenToShuffle - i) + i;
            Object swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }

    public static void Shuffle(Object[] arr, Random rn) {
        for (int i = 0; i < arr.length; i++) {
            int iSwap = rn.nextInt(arr.length - i) + i;
            Object swap = arr[iSwap];
            arr[iSwap] = arr[i];
            arr[i] = swap;
        }
    }

    /**
     * returns the distance squared between the two position provided in 2D
     *
     * @param x1 the xDim coordinate of the first position
     * @param y1 the yDim coordinate of the first position
     * @param x2 the xDim coordinate of the second position
     * @param y2 the yDim coordinate of the second position
     * @return the distance squared between the first and second position
     */
    public static double DistSq2D(double x1, double y1, double x2, double y2) {
        double xDist = x2 - x1, yDist = y2 - y1;
        return xDist * xDist + yDist * yDist;
    }
    public static double DistSq2D(double x1, double y1, double x2, double y2, double xDim, double yDim, boolean wrapX, boolean wrapY) {
        double xDist,yDist;
        if(wrapX){
            xDist= DistWrap(x1,x2,xDim);
        }
        else {
            xDist = x2 - x1;
        }
        if(wrapY){
            yDist= DistWrap(y1,y2,yDim);
        }
        else{
            yDist = y2 - y1;
        }
        return xDist * xDist + yDist * yDist;
    }

    /**
     * @param x1 the xDim coordinate of the first position
     * @param y1 the yDim coordinate of the first position
     * @param x2 the xDim coordinate of the second position
     * @param y2 the yDim coordinate of the second position
     * @return the distance squared between the first and second position
     */
    public static double DistSq3D(double x1,double y1,double z1, double x2, double y2, double z2) {
        double xDist = x2 - x1, yDist = y2 - y1,zDist=z2-z1;
        return xDist * xDist + yDist * yDist + zDist*zDist;
    }
    public static double DistSq3D(double x1, double y1, double z1, double x2, double y2, double z2, int xDim, int yDim, int zDim, boolean wrapX, boolean wrapY, boolean wrapZ) {
        double xDist,yDist,zDist;
        if(wrapX){
            xDist= DistWrap(x1,x2,xDim);
        }
        else {
            xDist = x2 - x1;
        }
        if(wrapY){
            yDist= DistWrap(y1,y2,yDim);
        }
        else{
            yDist = y2 - y1;
        }
        if(wrapZ){
            zDist= DistWrap(z1,z2,zDim);
        }
        else{
            zDist = z2 - z1;
        }
        return xDist * xDist + yDist * yDist + zDist * zDist;
    }

    /**
     * returns the distance squared between the two position provided in any number of dimensions
     *
     * @param p1 the coordinates of the first position
     * @param p2 the coordinates of the second position
     * @return the distance squared between the first and second position
     */
    public static double DistSqND(double[] p1, double[] p2) {
        double sum = 0;
        for (int i = 0; i < p1.length; i++) {
            double diff = p1[i] - p2[i];
            sum += diff * diff;
        }
        return sum;
    }
    public static double Norm(double v1,double v2){
        return Math.sqrt((v1*v1)+(v2*v2));
    }
    public static double Norm(double v1,double v2,double v3){
        return Math.sqrt((v1*v1)+(v2*v2)+(v3*v3));
    }
    public static double Norm(double v1,double v2,double v3,double v4){
        return Math.sqrt((v1*v1)+(v2*v2)+(v3*v3)+(v4*v4));
    }
    public static double Norm(double... vals){
        double tot=0;
        for (double val : vals) {
            tot+=val*val;
        }
        return Math.sqrt(tot);
    }
    public static double NormSq(double v1,double v2){
        return (v1*v1)+(v2*v2);
    }
    public static double NormSq(double v1,double v2,double v3){
        return (v1*v1)+(v2*v2)+(v3*v3);
    }
    public static double NormSq(double v1,double v2,double v3,double v4){
        return (v1*v1)+(v2*v2)+(v3*v3)+(v4*v4);
    }

    /**
     * returns the evenly spaced coordinates along the edge of a circle in 2D, centered on (0,0)
     */
    public static float[] GenCirclePoints(final float rad,final int nCorners){
        float[]ret=new float[nCorners*2+2];
        double step=(Math.PI*2)/(nCorners);
        double pos=0;
        for (int i = 0; i < nCorners+1; i++) {
            ret[i*2]=(float)Math.cos(pos)*rad;
            ret[i*2+1]=(float)Math.sin(pos)*rad;
            pos+=step;
        }
        return ret;
    }

//    public static <T extends AgentPT2,G extends Grid2<T>> void AgentsInRad(G searchMe, final ArrayList<T> putHere, final double x, final double y, final double rad, boolean wrapX, boolean wrapY){
//        putHere.clear();
//        int nAgents;
//        for (int xSq = (int)Math.floor(x-rad); xSq <(int)Math.ceil(x+rad) ; xSq++) {
//            for (int ySq = (int)Math.floor(y-rad); ySq <(int)Math.ceil(y+rad) ; ySq++) {
//                int retX=xSq; int retY=ySq;
//                boolean inX=Utils.InDim(searchMe.xDim,retX);
//                boolean inY=Utils.InDim(searchMe.yDim,retY);
//                if((!wrapX&&!inX)||(!wrapY&&!inY)){
//                    continue;
//                }
//                if(wrapX&&!inX){
//                    retX=Utils.ModWrap(retX,searchMe.xDim);
//                }
//                if(wrapY&&!inY){
//                    retY=Utils.ModWrap(retY,searchMe.yDim);
//                }
//                searchMe.GetAgents(putHere,searchMe.I(retX,retY));
//            }
//        }
//    }
//
    public static double DistWrap(double p1, double p2, double dim){
        if(Math.abs(p2-p1)>dim/2){
            if(p1>p2){
                p2=p2+dim;
            }
            else{
                p2=p2-dim;
            }
        }
        return p2-p1;
    }

//    public static <T extends AgentPhys2,Q extends AgentPhys2,G extends Grid2<Q>>double CollisionSum2D(T agent,G searchMe, final ArrayList<Q> putAgentsHere,RadToForceMap ForceFun,double searchRad,boolean wrapX,boolean wrapY){
//        double ret=0;
//        putAgentsHere.clear();
//        AgentsInRad(searchMe,putAgentsHere,agent.Xpt(),agent.Ypt(),searchRad,wrapX,wrapY);
//        for (Q a : putAgentsHere) {
//            if(a!=agent){
//                double xComp=wrapX?DistWrap(agent.Xpt(), a.Xpt(), searchMe.xDim):a.Xpt()-agent.Xpt();
//                double yComp=wrapY?DistWrap(agent.Ypt(),a.Ypt(),searchMe.yDim):a.Ypt()-agent.Ypt();
//                double dist=Math.sqrt(xComp*xComp+yComp*yComp)-(agent.radius+a.radius);
//                double force=ForceFun.DistToForce(dist);
//                agent.AddForce(xComp,yComp,force);
//                ret+=force;
//            }
//        }
//        return ret;
//    }

    /**
     * returns the original value bounded by min and max inclusive
     */
    public static double BoundVal(double val, double min, double max) {
        return val < min ? min : (val > max ? max : val);
    }

    /**
     * returns the original value bounded by min and max inclusive
     */
    public static int BoundValI(int val, int min, int max) {
        return val < min ? min : (val > max ? max : val);
    }

    /**
     * returns the original value bounded by min and max inclusive
     */
    public static float BoundValF(float val, double min, double max) {
        return (float) (val < min ? min : (val > max ? max : val));
    }

    /**
     * returns where the value is from min to max as a number from 0 to 1
     */
    public static double ScaleVal(double val, double min, double max) {
        return (val - min) / (max - min);
    }

    /**
     * returns value with wraparound between 0 and max
     */
    public static int ModWrap(int val, int max) {
        return val < 0 ? max + val % max : val % max;
    }
    /**
     * returns value with wraparound between 0 and max
     */
    public static double ModWrap(double val, double max) {
        return val < 0 ? max + val % max : val % max;
    }

    /**
     * converts proton concentration to pH
     */
    public static double ConvertHToPh(double h) {
        return -Math.log10(h) + 3;
    }

    /**
     * converts pH to proton concentration
     */
    public static double ConvertPhToH(double ph) {
        return Math.pow(10, 3.0 - ph);
    }

    /**
     * adjusts probability that an event will occur in 1 unit of time to the probability that the event will occur in timeFraction duration
     *
     * @param prob         probability that an event occurs in 1 unit of time
     * @param timeFraction duration to over which event may occur
     * @return the probability that the event will occur in timeFraction
     */
    public static double ProbScale(double prob, double timeFraction) {
        return 1.0f - (Math.pow(1.0 - prob, timeFraction));

    }

    //LIST FUNCTIONS
    public static double LogDist(double min,double max,Random rn){
        if(min<=0||max<=0){
            System.err.println("Error, LogDist contains range value < 0!");
        }
        double logMin=Math.log(min);
        double logMax=Math.log(max);
        double inVal=rn.nextDouble()*(logMax-logMin)+logMin;
        return Math.exp(inVal);
    }


    /**
     * runs a fully connected neural network layer
     *
     * @param neurons      array of all neurons in the network
     * @param weights      array of weights in the fully connected layer
     * @param iFromStart   index of the start of the input to the fully connected layer
     * @param iFromEnd     index of the end of the input to the fully connected layer
     * @param iToStart     index of the start of the output of the fully connected layer
     * @param iToEnd       index of the end of the output of the fully connected layer
     * @param iWeightStart index of the start of the weights for this layer out of the weights array
     */
    public static void NNfullyConnectedLayer(double[] neurons, double[] weights, int iFromStart, int iFromEnd, int iToStart, int iToEnd, int iWeightStart) {
        int iWeight = iWeightStart;
        for (int iFrom = iFromStart; iFrom < iFromEnd; iFrom++) {
            for (int iTo = iToStart; iTo < iToEnd; iTo++) {
                neurons[iTo] += neurons[iFrom] * weights[iWeight];
                iWeight++;
            }
        }
    }

    /**
     * set all neurons between iStart and iEnd with the given value
     */
    public static void NNset(double[] neurons, int iStart, int iEnd, double val) {
        Arrays.fill(neurons, iStart, iEnd, val);
    }

    //UTILITIES

    /**
     * returns a timestamp of the form "yyyy_MM_dd_HH_mm_ss" as a string
     */
    static public String TimeStamp() {
        return new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss").format(new Date());

    }

    /**
     * gets the current working directory as a string
     */
    static public String PWD() {
        return Paths.get("").toAbsolutePath().toString() + "/";
        //return System.getProperty("user.dir");
    }

    /**
     * gets information about the memory usage and max memory allocated for the program
     */
    public static String MemoryUsageStr() {
        int mb = 1024 * 1024;
        Runtime rt = Runtime.getRuntime();
        return "Used Memory: " + (rt.totalMemory() - rt.freeMemory()) / mb + " mb, Free Momory: " + rt.freeMemory() / mb + " mb, Total Memory:" + rt.totalMemory() / mb + " mb, Max  Memory: " + rt.maxMemory() / mb + " mb";
    }
    /**
     * Runs quicksort on an object that implements Sortable
     *
     * @param sortMe          the object to be sorted
     * @param greatestToLeast if true, sorting will be form greatest to least, otherwise will be least to greatest
     */
    public static <T extends Sortable> void QuickSort(T sortMe, boolean greatestToLeast) {
        SortHelper(sortMe, 0, sortMe.Length() - 1, greatestToLeast);
    }

    static <T extends Sortable> void SortHelper(T sortMe, int lo, int hi, boolean greatestToLeast) {
        if (lo < hi) {
            int p = Partition(sortMe, lo, hi, greatestToLeast);
            SortHelper(sortMe, lo, p - 1, greatestToLeast);
            SortHelper(sortMe, p + 1, hi, greatestToLeast);
        }
    }

    static <T extends Sortable> int Partition(T sortMe, int lo, int hi, boolean greatestToLeast) {
        if (greatestToLeast) {
            for (int j = lo; j < hi; j++) {
                if (sortMe.Compare(hi, j) <= 0) {
                    sortMe.Swap(lo, j);
                    lo++;
                }
            }
            sortMe.Swap(lo, hi);
            return lo;
        } else {
            for (int j = lo; j < hi; j++) {
                if (sortMe.Compare(hi, j) >= 0) {
                    sortMe.Swap(lo, j);
                    lo++;
                }
            }
            sortMe.Swap(lo, hi);
            return lo;
        }
    }

    /**
     * returns whether the input value is between 0 and the dimension value
     */
    public static boolean InDim(int Dim, int Val) {
        return Val >= 0 && Val < Dim;
    }
    /**
     * returns whether the input value is between 0 and the dimension value
     */
    public static boolean InDim(double Dim, double Val) {
        return Val >= 0 && Val < Dim;
    }
    public double GradientX(double[] vals,int xDim,int yDim,int centerX,int centerY,boolean boundaryCond,double boundaryValue,boolean wrapX){
        double xP1,xM1;
        if(InDim(xDim,centerX+1)){
            xP1=vals[(centerX+1)*yDim+centerY];
        }
        else if(boundaryCond){
            xP1=boundaryValue;
        }
        else if(wrapX){
            xP1=vals[(0)*yDim+centerY];
        }
        else{
            xP1=vals[centerX*yDim+centerY];
        }
        if(InDim(xDim,centerX-1)){
            xM1=vals[(centerX-1)*yDim+centerY];
        }
        else if(boundaryCond){
            xM1=boundaryValue;
        }
        else if(wrapX){
            xM1=vals[(xDim-1)*yDim+centerY];
        }
        else{
            xM1=vals[centerX*yDim+centerY];
        }
        return xP1-xM1;
    }
    public double GradientY(double[] vals,int xDim,int yDim,int centerX,int centerY,boolean boundaryCond,double boundaryValue,boolean wrapY){
        double yP1,yM1;
        if(InDim(yDim,centerY+1)){
            yP1=vals[centerX*yDim+(centerY+1)];
        }
        else if(boundaryCond){
            yP1=boundaryValue;
        }
        else if(wrapY){
            yP1=vals[centerX*yDim+(0)];
        }
        else{
            yP1=vals[centerX*yDim+centerY];
        }
        if(InDim(yDim,centerY-1)){
            yM1=vals[centerX*yDim+(centerY-1)];
        }
        else if(boundaryCond){
            yM1=boundaryValue;
        }
        else if(wrapY){
            yM1=vals[centerX*yDim+(yDim-1)];
        }
        else{
            yM1=vals[centerX*yDim+centerY];
        }
        return yP1-yM1;
    }


    public static <T> ArrayList<T> ParallelSweep(int nRuns, int nThreads, SweepRunFunction<T> RunFn) {
        ArrayList<T> runOuts = new ArrayList<>(nRuns);
        ArrayList<SweepRun<T>> runners = new ArrayList<>(nRuns);
        for (int i = 0; i < nRuns; i++) {
            runOuts.add(null);
            runners.add(new SweepRun<T>(RunFn, runOuts, i));
        }
        ExecutorService exec = Executors.newFixedThreadPool(nThreads);
        for (SweepRun<T> run : runners) {
            exec.execute(run);
        }
        exec.shutdown();
        while (!exec.isTerminated()) ;
        return runOuts;
    }

    /**
     * runs the finite differences equation in 2 dimensions
     *
     * @param inGrid        an array of values holding the starting state of the diffusible
     * @param outGrid       an array into which the result of diffusion will be written
     * @param xDim          xDim dimenison of the inGrid and outGrid
     * @param yDim          yDim dimension of the inGrid and outGrid
     * @param diffRate      diffusion rate for the diffusion equaition
     * @param boundaryCond  defines whether a constant boundary condition value should diffuse in from the boundaries
     * @param boundaryValue only impacts diffusion if boundaryCond is true, sets the boundary condition value
     * @param wrapX         whether to wrap around diffusion over the left and right boundaries
     */
    public static void Diffusion(double[] inGrid, double[] outGrid, int xDim, int yDim, double diffRate, boolean boundaryCond, double boundaryValue, final boolean wrapX, final boolean wrapY) {
        //This code is ugly and repetitive to improve performance by getting around bounds checking
        int x, y;
        //first we do the corners
        if (boundaryCond) {
            outGrid[0] = inGrid[0] + diffRate * (-inGrid[0] * 4 + inGrid[1] + inGrid[yDim] + 2 * boundaryValue);
            outGrid[(xDim - 1) * yDim] = inGrid[(xDim - 1) * yDim] + diffRate * (-inGrid[(xDim - 1) * yDim] * 4 + inGrid[(xDim - 2) * yDim] + inGrid[(xDim - 1) * yDim + 1] + 2 * boundaryValue);
            outGrid[(xDim - 1) * yDim + yDim - 1] = inGrid[(xDim - 1) * yDim + yDim - 1] + diffRate * (-inGrid[(xDim - 1) * yDim + yDim - 1] * 4 + inGrid[(xDim - 2) * yDim + yDim - 1] + inGrid[(xDim - 1) * yDim + yDim - 2] + 2 * boundaryValue);
            outGrid[yDim - 1] = inGrid[yDim - 1] + diffRate * (-inGrid[yDim - 1] * 4 + inGrid[yDim + yDim - 1] + inGrid[yDim - 2] + 2 * boundaryValue);
        } else {
            outGrid[0] = inGrid[0] + diffRate * (-inGrid[0] * 4 + (wrapX?inGrid[1*yDim+0]+inGrid[(xDim-1)*yDim]:inGrid[(1*yDim)+0]*2) + (wrapY?inGrid[0*yDim+yDim]+inGrid[0*yDim+1]:inGrid[0*yDim+1]*2));
            outGrid[(xDim - 1) * yDim] = inGrid[(xDim - 1) * yDim] + diffRate * (-inGrid[(xDim - 1) * yDim] * 4 + (wrapX?inGrid[(xDim - 2) * yDim]+inGrid[0]:2*inGrid[(xDim - 2)*yDim]) + (wrapY?inGrid[(xDim - 1) * yDim + 1]+inGrid[xDim*yDim-1]:2*inGrid[(xDim - 1) * yDim + 1]));
            outGrid[(xDim - 1) * yDim + yDim - 1] = inGrid[(xDim - 1) * yDim + yDim - 1] + diffRate * (-inGrid[(xDim - 1) * yDim + yDim - 1] * 4 + (wrapX?inGrid[(xDim - 2) * yDim + yDim - 1]+inGrid[yDim-1]:2*inGrid[(xDim - 2) * yDim + yDim - 1]) + (wrapY?inGrid[(xDim-1)*yDim]+inGrid[(xDim - 1) * yDim + yDim - 2]:2*inGrid[(xDim - 1) * yDim + yDim - 2]));
            outGrid[yDim - 1] = inGrid[yDim - 1] + diffRate * (-inGrid[yDim - 1] * 4 + (wrapX?inGrid[yDim + yDim - 1]+inGrid[(xDim-1)*yDim+yDim-1]:2*inGrid[yDim + yDim - 1]) + (wrapY?inGrid[0]+inGrid[yDim - 2]:2*inGrid[yDim - 2]));
        }
        //then we do the sides
        if (boundaryCond) {
            x = 0;
            for (y = 1; y < yDim - 1; y++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[(x + 1) * yDim + y] + inGrid[x * yDim + y + 1] + inGrid[x * yDim + y - 1] + boundaryValue);
            }
            x = xDim - 1;
            for (y = 1; y < yDim - 1; y++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[(x - 1) * yDim + y] + inGrid[x * yDim + y + 1] + inGrid[x * yDim + y - 1] + boundaryValue);
            }
            y = 0;
            for (x = 1; x < xDim - 1; x++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[x * yDim + y + 1] + inGrid[(x + 1) * yDim + y] + inGrid[(x - 1) * yDim + y] + boundaryValue);
            }
            y = yDim - 1;
            for (x = 1; x < xDim - 1; x++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[x * yDim + y - 1] + inGrid[(x + 1) * yDim + y] + inGrid[(x - 1) * yDim + y] + boundaryValue);
            }
        } else{
            x = 0;
            for (y = 1; y < yDim - 1; y++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[(x + 1) * yDim + y] + inGrid[x * yDim + y + 1] + inGrid[x * yDim + y - 1] + (wrapX?inGrid[(xDim - 1) * yDim + y]:inGrid[(x + 1) * yDim + y]));
            }
            x = xDim - 1;
            for (y = 1; y < yDim - 1; y++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[(x - 1) * yDim + y] + inGrid[x * yDim + y + 1] + inGrid[x * yDim + y - 1] + (wrapX?inGrid[0 * yDim + y]:inGrid[(x-1)*yDim+y]));
            }
            y = 0;
            for (x = 1; x < xDim - 1; x++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[x * yDim + y + 1] + inGrid[(x + 1) * yDim + y] + inGrid[(x - 1) * yDim + y]+(wrapY?inGrid[x*yDim+yDim-1]:inGrid[x*yDim+1]));
            }
            y = yDim - 1;
            for (x = 1; x < xDim - 1; x++) {
                outGrid[x * yDim + y] = inGrid[x * yDim + y] + diffRate * (-inGrid[x * yDim + y] * 4 + inGrid[x * yDim + y - 1] + inGrid[(x + 1) * yDim + y] + inGrid[(x - 1) * yDim + y]+(wrapY?inGrid[x*yDim]:inGrid[x*yDim+y-2]));
            }
        }
        //then we do the middle
        for (x = 1; x < xDim - 1; x++) {
            for (y = 1; y < yDim - 1; y++) {
                int i = x * yDim + y;
                outGrid[i] = inGrid[i] + diffRate * (-inGrid[i] * 4 + inGrid[(x + 1) * yDim + y] + inGrid[(x - 1) * yDim + y] + inGrid[x * yDim + y + 1] + inGrid[x * yDim + y - 1]);
            }
        }
    }

    //boolean in3(int compX,int compY,int zDim,int xDim,int yDim,int z){
    //    return xDim>=0&&xDim<compX&&yDim>=0&&yDim<compY&&z>=0&&z<zDim;
    //}


    /**
     * runs the diffusion equation in 3 dimensions
     *
     * @param inGrid        an array of values holding the starting state of the diffusible
     * @param outGrid       an array into which the result of diffusion will be written
     * @param xDim          xDim dimenison of the inGrid and outGrid
     * @param yDim          yDim dimension of the inGrid and outGrid
     * @param zDim          z dimension of the inGrid and outGrid
     * @param diffRate      diffusion rate for the diffusion equaition
     * @param boundaryCond  defines whether a constant boundary condition value should diffuse in from the boundaries
     * @param boundaryValue only impacts diffusion if boundaryCond is true, sets the boundary condition value
     * @param wrapX        whether to wrap around diffusion over the left and right and front and back boundaries
     * @param wrapY        whether to wrap around diffusion over the left and right and front and back boundaries
     * @param wrapZ        whether to wrap around diffusion over the left and right and front and back boundaries
     */
    public static void Diffusion3(final double[] inGrid, final double[] outGrid, final int xDim, final int yDim, final int zDim, final double
            diffRate, final boolean boundaryCond, final double boundaryValue, final boolean wrapX,final boolean wrapY,final boolean wrapZ) {
        int x, y, z;
        double valSum;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {
                for (z = 0; z < zDim; z++) {
                    //6 squares to check
                    valSum = 0;
                    if (InDim(xDim, x + 1)) {
                        valSum += inGrid[(x + 1) * yDim * zDim + (y) * zDim + (z)];
                    }  else if (boundaryCond) {
                        valSum += boundaryValue;
                    }else if (wrapX) {
                        valSum += inGrid[(0) * yDim * zDim + (y) * zDim + (z)];
                    } else{
                        valSum+=inGrid[(x-1)*yDim*zDim+(y)*zDim+(z)];
                    }

                    if (InDim(xDim, x - 1)) {
                        valSum += inGrid[(x - 1) * yDim * zDim + (y) * zDim + (z)];
                    }  else if (boundaryCond) {
                        valSum += boundaryValue;
                    }else if (wrapX) {
                        valSum += inGrid[(xDim - 1) * yDim * zDim + (y) * zDim + (z)];
                    } else{
                        valSum += inGrid[(x + 1) * yDim * zDim + (y) * zDim + (z)];
                    }

                    if (InDim(yDim, y + 1)) {
                        valSum += inGrid[(x) * yDim * zDim + (y + 1) * zDim + (z)];
                    } else if (boundaryCond) {
                        valSum += boundaryValue;
                    } else if(wrapY){
                        valSum += inGrid[(x) * yDim * zDim + (0) * zDim + (z)];
                    } else{
                        valSum += inGrid[(x) * yDim * zDim + (y-1) * zDim + (z)];
                    }

                    if (InDim(yDim, y - 1)) {
                        valSum += inGrid[(x) * yDim * zDim + (y - 1) * zDim + (z)];
                    } else if (boundaryCond) {
                        valSum += boundaryValue;
                    } else if(wrapY){
                        valSum += inGrid[(x) * yDim * zDim + (yDim-1) * zDim + (z)];
                    } else{
                        valSum += inGrid[(x) * yDim * zDim + (y+1) * zDim + (z)];
                    }

                    if (InDim(zDim, z + 1)) {
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (z+1)];
                    } else if (boundaryCond) {
                        valSum += boundaryValue;
                    } else if(wrapZ){
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (0)];
                    } else{
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (z-1)];
                    }

                    if (InDim(zDim, z - 1)) {
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (z-1)];
                    } else if (boundaryCond) {
                        valSum += boundaryValue;
                    } else if(wrapZ){
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (zDim-1)];
                    } else{
                        valSum += inGrid[(x) * yDim * zDim + (y) * zDim + (z+1)];
                    }
                    int i = x * yDim * zDim + y * zDim + z;
                    outGrid[i] = inGrid[i] + diffRate * (-inGrid[i] * 6 + valSum);
                }
            }
        }
    }
    public static void Diffusion2inhomogeneous(final double[] inGrid, final double[] outGrid,final double[] diffRates, final int xDim, final int yDim, final boolean boundaryCond, final double boundaryValue, final boolean wrapX,final boolean wrapY) {
        int x, y;
        double valSum;
        double rateSum;
        double currRate;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {

                int i = x * yDim + y;
                //4 squares to check
                valSum = 0;
                rateSum = 0;
                double diffRate=diffRates[x*yDim+y];

                if (InDim(xDim, x + 1)) {
                    currRate=diffRate+diffRates[(x+1)*yDim+(y)];
                    valSum += inGrid[(x + 1) * yDim + (y)]*currRate;
                    rateSum+=currRate;
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                    rateSum++;
                } else if (wrapX) {
                    currRate=diffRate+diffRates[(0) * yDim + (y)];
                    valSum += inGrid[(0) * yDim + (y)]*currRate;
                    rateSum+=currRate;
                } else {
                    currRate=diffRate+diffRates[(x - 1) * yDim + (y)];
                    valSum += inGrid[(x - 1) * yDim + (y)]*currRate;
                    rateSum += currRate;
                }

                if (InDim(xDim, x - 1)) {
                    currRate=diffRate+diffRates[(x - 1) * yDim + (y)];
                    valSum += inGrid[(x - 1) * yDim + (y)]*currRate;
                    rateSum+=currRate;
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                    rateSum++;
                } else if (wrapX) {
                    currRate=diffRate+diffRates[(xDim - 1) * yDim + (y)];
                    valSum += inGrid[(xDim - 1) * yDim + (y)]*currRate;
                    rateSum+=currRate;
                } else {
                    currRate=diffRate+diffRates[(x + 1) * yDim + (y)];
                    valSum += inGrid[(x + 1) * yDim + (y)]*currRate;
                    rateSum+=currRate;
                }

                if (InDim(yDim, y + 1)) {
                    currRate=diffRate+diffRates[(x) * yDim + (y + 1)];
                    valSum += inGrid[(x) * yDim + (y + 1)]*currRate;
                    rateSum+=currRate;
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                    rateSum++;
                } else if (wrapY) {
                    currRate=diffRate+diffRates[(x) * yDim + (0)];
                    valSum += inGrid[(x) * yDim + (0)]*currRate;
                    rateSum+=currRate;
                } else {
                    currRate=diffRate+diffRates[(x) * yDim + (y - 1)];
                    valSum += inGrid[(x) * yDim + (y - 1)]*currRate;
                    rateSum+=currRate;
                }

                if (InDim(yDim, y - 1)) {
                    currRate=diffRate+diffRates[(x) * yDim + (y - 1)];
                    valSum += inGrid[(x) * yDim + (y - 1)]*currRate;
                    rateSum+=currRate;
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                    rateSum++;
                } else if (wrapY) {
                    currRate=diffRate+diffRates[(x) * yDim + (yDim - 1)];
                    valSum += inGrid[(x) * yDim + (yDim - 1)]*currRate;
                    rateSum+=currRate;
                } else {
                    currRate=diffRate+diffRates[(x) * yDim + (y + 1)];
                    valSum += inGrid[(x) * yDim + (y + 1)]*currRate;
                    rateSum+=currRate;
                }
                outGrid[i] = inGrid[i] +  (-inGrid[i] * rateSum + valSum)/2;
            }
        }
    }
    public static void Diffusion2(final double[] inGrid, final double[] outGrid, final int xDim, final int yDim, final double diffRate, final boolean boundaryCond, final double boundaryValue, final boolean wrapX,final boolean wrapY) {
        int x, y;
        double valSum;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {
                //4 squares to check
                valSum = 0;
                if (InDim(xDim, x + 1)) {
                    valSum += inGrid[(x + 1) * yDim + (y)];
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                } else if (wrapX) {
                    valSum += inGrid[(0) * yDim + (y)];
                } else {
                    valSum += inGrid[(x - 1) * yDim + (y)];
                }

                if (InDim(xDim, x - 1)) {
                    valSum += inGrid[(x - 1) * yDim + (y)];
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                } else if (wrapX) {
                    valSum += inGrid[(xDim - 1) * yDim + (y)];
                } else {
                    valSum += inGrid[(x + 1) * yDim + (y)];
                }

                if (InDim(yDim, y + 1)) {
                    valSum += inGrid[(x) * yDim + (y + 1)];
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                } else if (wrapY) {
                    valSum += inGrid[(x) * yDim + (0)];
                } else {
                    valSum += inGrid[(x) * yDim + (y - 1)];
                }

                if (InDim(yDim, y - 1)) {
                    valSum += inGrid[(x) * yDim + (y - 1)];
                } else if (boundaryCond) {
                    valSum += boundaryValue;
                } else if (wrapY) {
                    valSum += inGrid[(x) * yDim + (yDim - 1)];
                } else {
                    valSum += inGrid[(x) * yDim + (y + 1)];
                }
                int i = x * yDim + y;
                outGrid[i] = inGrid[i] + diffRate * (-inGrid[i] * 4 + valSum);
            }
        }
    }
    public static void Advection2(final double[] inGrid, final double[] outGrid, final int xDim, final int yDim, final double[] xVels, final double[] yVels, final boolean boundaryCond, final double boundaryValue) {
        //mostly works as intended, value can be deleted if velocity direction changes
        int x, y;
        double xComp;
        double yComp;
        double currVal;
        double xVel;
        double yVel;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {
                //2 squares to check
                xComp = 0;
                yComp = 0;
                int i = x * yDim + y;
                xVel=xVels[i];
                yVel=yVels[i];

                currVal=inGrid[i];
                if(xVel!=0) {
                    if(xVel<0) {
                        if (InDim(xDim, x + 1)) {
                            xComp = inGrid[(x + 1) * yDim + (y)];
                        } else if (boundaryCond) {
                            xComp = boundaryValue;
                        } else {
                            xComp = inGrid[(0) * yDim + (y)];
                        }
                        xComp=xVel*(currVal-xComp);
                    }
                    else {
                        if (InDim(xDim, x - 1)) {
                            xComp = inGrid[(x - 1) * yDim + (y)];
                        } else if (boundaryCond) {
                            xComp = boundaryValue;
                        } else {
                            xComp = inGrid[(xDim - 1) * yDim + (y)];
                        }
                        xComp=xVel*(xComp-currVal);
                    }
                }
                if(yVel!=0) {
                    if (yVel < 0) {
                        if (InDim(yDim, y + 1)) {
                            yComp += inGrid[(x) * yDim + (y + 1)];
                        } else if (boundaryCond) {
                            yComp += boundaryValue;
                        } else {
                            yComp += inGrid[(x) * yDim + (0)];
                        }
                        yComp=yVel*(currVal-yComp);
                    } else {
                        if (InDim(yDim, y - 1)) {
                            yComp += inGrid[(x) * yDim + (y - 1)];
                        } else if (boundaryCond) {
                            yComp += boundaryValue;
                        } else {
                            yComp += inGrid[(x) * yDim + (yDim - 1)];
                        }
                        yComp=yVel*(yComp-currVal);
                    }
                }
                outGrid[i] = inGrid[i] + xComp + yComp;
            }
        }
    }
    public static void ConservativeTransportStep2(final double[] inGrid,final double[] midGrid, final double[] outGrid, final int xDim, final int yDim,final double[] xVels,final double[] yVels, final boolean boundaryCond, final double boundaryValue){
        int x,y;
        double vXp1,vXm1,vYp1,vYm1,dVx,dVy;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {
                int i=x*yDim+y;
                if(InDim(xDim,x+1)){
                    vXp1=xVels[(x+1)*yDim+y];
                }
                else if(boundaryCond){
                    vXp1=boundaryValue;
                }
                else{
                    vXp1=xVels[(0)*yDim+y];
                }
                if(InDim(xDim,x-1)){
                    vXm1=xVels[(x-1)*yDim+y];
                }
                else if(boundaryCond){
                    vXm1=boundaryValue;
                }
                else{
                    vXm1=xVels[(xDim-1)*yDim+y];
                }

                if(InDim(yDim,y+1)){
                    vYp1=yVels[x*yDim+(y+1)];
                }
                else if(boundaryCond){
                    vYp1=boundaryValue;
                }
                else{
                    vYp1=yVels[x*yDim+(0)];
                }
                if(InDim(yDim,y-1)){
                    vYm1=yVels[x*yDim+(y-1)];
                }
                else if(boundaryCond){
                    vYm1=boundaryValue;
                }
                else{
                    vYm1=yVels[x*yDim+(yDim-1)];
                }
                dVx=vXp1-vXm1;
                dVy=vYp1-vYm1;
                outGrid[i]=midGrid[i]-(dVx/2+dVy/2)*inGrid[i];
            }
        }
    }
    public static void Advection2(final double[] inGrid, final double[] outGrid, final int xDim, final int yDim, final double xVel, final double yVel, final boolean boundaryCond, final double boundaryValue) {
        int x, y;
        double xComp;
        double yComp;
        double currVal;
        for (x = 0; x < xDim; x++) {
            for (y = 0; y < yDim; y++) {
                //2 squares to check
                xComp = 0;
                yComp = 0;
                int i = x * yDim + y;
                currVal=inGrid[i];
                if(xVel!=0) {
                    if(xVel<0) {
                        if (InDim(xDim, x + 1)) {
                            xComp = inGrid[(x + 1) * yDim + (y)];
                        } else if (boundaryCond) {
                            xComp = boundaryValue;
                        } else {
                            xComp = inGrid[(0) * yDim + (y)];
                        }
                        xComp=xVel*(currVal-xComp);
                    }
                    else {
                        if (InDim(xDim, x - 1)) {
                            xComp = inGrid[(x - 1) * yDim + (y)];
                        } else if (boundaryCond) {
                            xComp = boundaryValue;
                        } else{
                            xComp = inGrid[(xDim - 1) * yDim + (y)];
                        }
                        xComp=xVel*(xComp-currVal);
                    }
                }
                if(yVel!=0) {
                    if (yVel < 0) {
                        if (InDim(yDim, y + 1)) {
                            yComp += inGrid[(x) * yDim + (y + 1)];
                        } else if (boundaryCond) {
                            yComp += boundaryValue;
                        } else{
                            yComp += inGrid[(x) * yDim + (0)];
                        }
                        yComp=yVel*(currVal-yComp);
                    } else {
                        if (InDim(yDim, y - 1)) {
                            yComp += inGrid[(x) * yDim + (y - 1)];
                        } else if (boundaryCond) {
                            yComp += boundaryValue;
                        } else {
                            yComp += inGrid[(x) * yDim + (yDim - 1)];
                        }
                        yComp=yVel*(yComp-currVal);
                    }
                }
                outGrid[i] = inGrid[i] + xComp + yComp;
            }
        }
    }
public static void TDMAx(final double[] in, final double[] out,final double[] scratch,final int xDim, final int yDim, final int iRow, final double diffRate) {
    final double ac=diffRate;
    final double b=-2*(diffRate+1);
    final double db=2*(diffRate-1);
    final double dac=-diffRate;

    final int len =  xDim;
    final int max = yDim;

    //Doing the 0 entries
    scratch[0] = (2.0*ac)/b;
    double above = iRow == max - 1 ? in[(0)*yDim+ (iRow - 1)]: in[(0)*yDim+(iRow + 1)];
    double below = iRow == 0 ? above : in[(0 *yDim)+(iRow - 1)];
    double middle=in[(0*yDim)+iRow];
    double di=db*middle + above*dac + below*dac;
    scratch[len] = (di) / b;

    //Doing the forward passes
    for (int i = 1; i < len-1; i++) {
        scratch[i] = ac / (b - ac*scratch[i - 1]);
    }
    for (int i = 1; i < len; i++) {
        above = iRow == max - 1 ? in[(i)*yDim +(iRow - 1)] : in[(i)*yDim+ (iRow + 1)];
        below = iRow == 0 ? above : in[(i)*yDim+(iRow - 1)];
        middle=in[(i)*yDim+(iRow)];
        di=(db*middle + above*dac + below*dac);
        if(i<len-1){scratch[len + i] = (di - ac*scratch[len + i - 1]) / (b - ac*scratch[i - 1]);}
        else{scratch[len + i] = (di - 2*ac*scratch[len + i - 1]) / (b - 2*ac*scratch[i - 1]);}
    }

    //backward pass, do not touch!
    out[(len - 1)*yDim+(iRow)] = scratch[len * 2 - 1];
    for (int i = len - 2; i >= 0; i--) {
        out[(i)*yDim+(iRow)] = scratch[len + i] - scratch[i] * out[(i + 1)*yDim+(iRow)];
    }
}
    public static void TDMAy(final double[] in, final double[] out,final double[] scratch,final int xDim, final int yDim, final int iRow, final double diffRate) {
        final double ac=diffRate;
        final double b=-2*(diffRate+1);
        final double db=2*(diffRate-1);
        final double dac=-diffRate;


        final int len = yDim;
        final int max = xDim;

        //Doing the 0 entries
        scratch[0] = (2.0*ac)/b;
        double above = iRow == max - 1 ? in[(iRow - 1)*yDim+(0)] : in[(iRow + 1)*yDim+(0)];
        double below = iRow == 0 ? above : in[(iRow - 1)* yDim+(0)];
        double middle=in[(iRow)*yDim+(0)];
        double di=db*middle + above*dac + below*dac;
        scratch[len] = (di) / b;

        //Doing the forward passes
        for (int i = 1; i < len-1; i++) {
            scratch[i] = ac / (b - ac*scratch[i - 1]);
        }
        for (int i = 1; i < len; i++) {
            above = iRow == max - 1 ? in[(iRow - 1)*yDim+(i)] : in[(iRow + 1)*yDim+(i)];
            below = iRow == 0 ? above : in[(iRow - 1)*yDim+(i)];
            middle=in[(iRow)*yDim+(i)];
            di=(db*middle + above*dac + below*dac);
            if(i<len-1){scratch[len + i] = (di - ac*scratch[len + i - 1]) / (b - ac*scratch[i - 1]);}
            else{scratch[len + i] = (di - 2*ac*scratch[len + i - 1]) / (b - 2*ac*scratch[i - 1]);}
        }

        //backward pass, do not touch!
        out[(iRow)*yDim+(len - 1)] = scratch[len * 2 - 1];
        for (int i = len - 2; i >= 0; i--) {
            out[(iRow)*yDim+(i)] = scratch[len + i] - scratch[i] * out[(iRow)*yDim+(i + 1)];
        }
    }
    //only changes from 2D to 3D are 4x diffRate for db, front and back added to d side of equation, and array accesses are different
    public static void TDMA3x(final double[] in, final double[] out,final double[] scratch,final int xDim, final int yDim,final int zDim, final int iRow1,final int iRow2, final double diffRate) {
        //doing:x, iRow1:y, iRow2:z
        final double ac=diffRate;
        final double b=-2*(2*diffRate+1);
        final double db=2*(2*diffRate-1);
        final double dac=-diffRate;

        final int len =  xDim;
        final int max1 = yDim;
        final int max2 = zDim;

        //Doing the 0 entries
        scratch[0] = (2.0*ac)/b;
        double above = iRow1 == max1 - 1 ? in[(0)*yDim*zDim+ (iRow1 - 1)*yDim+iRow2]: in[(0)*yDim*zDim+(iRow1 + 1)*yDim+iRow2];
        double below = iRow1 == 0 ? above : in[(0)*yDim*zDim+(iRow1 - 1)*yDim+iRow2];
        double front = iRow2 == max2 - 1 ? in[(0)*yDim*zDim+ (iRow1)*yDim+iRow2-1]: in[(0)*yDim*zDim+(iRow1)*yDim+iRow2+1];
        double back = iRow1 == 0 ? above : in[(0)*yDim*zDim+(iRow1)*yDim+iRow2-1];
        double middle=in[(0)*yDim*zDim+iRow1*yDim+iRow2];
        double di=db*middle + above*dac + below*dac + front*dac + back*dac;
        scratch[len] = (di) / b;

        //Doing the forward passes
        for (int i = 1; i < len-1; i++) {
            scratch[i] = ac / (b - ac*scratch[i - 1]);
        }
        for (int i = 1; i < len; i++) {

            above = iRow1 == max1 - 1 ? in[(i)*yDim*zDim+ (iRow1 - 1)*yDim+iRow2]: in[(i)*yDim*zDim+(iRow1 + 1)*yDim+iRow2];
            below = iRow1 == 0 ? above : in[(i)*yDim*zDim+(iRow1 - 1)*yDim+iRow2];
            front = iRow2 == max2 - 1 ? in[(i)*yDim*zDim+ (iRow1)*yDim+iRow2-1]: in[(i)*yDim*zDim+(iRow1)*yDim+iRow2+1];
            back = iRow1 == 0 ? above : in[(i)*yDim*zDim+(iRow1)*yDim+iRow2-1];
            middle=in[(i)*yDim*zDim+iRow1*yDim+iRow2];
            di=db*middle + above*dac + below*dac + front*dac + back*dac;
            if(i<len-1){scratch[len + i] = (di - ac*scratch[len + i - 1]) / (b - ac*scratch[i - 1]);}
            else{scratch[len + i] = (di - 2*ac*scratch[len + i - 1]) / (b - 2*ac*scratch[i - 1]);}
        }

        //backward pass, do not touch!
        out[(len - 1)*yDim*zDim+(iRow1)*yDim+iRow2] = scratch[len * 2 - 1];
        for (int i = len - 2; i >= 0; i--) {
            out[(i)*yDim*zDim+(iRow1)*yDim+iRow2] = scratch[len + i] - scratch[i] * out[(i + 1)*yDim*zDim+(iRow1)*yDim+iRow2];
        }
    }
    //only changes from 2D to 3D are 4x diffRate for db, front and back added to d side of equation, and array accesses are different
    public static void TDMA3y(final double[] in, final double[] out,final double[] scratch,final int xDim, final int yDim,final int zDim, final int iRow1,final int iRow2, final double diffRate) {
        //doing:y, iRow1:x, iRow2:z
        final double ac=diffRate;
        final double b=-2*(2*diffRate+1);
        final double db=2*(2*diffRate-1);
        final double dac=-diffRate;

        final int max1 = xDim;
        final int len =  yDim;
        final int max2 = zDim;

        //Doing the 0 entries
        scratch[0] = (2.0*ac)/b;
        double above = iRow1 == max1 - 1 ? in[(iRow1 - 1)*yDim*zDim+ (0)*yDim+iRow2]: in[(iRow1+1)*yDim*zDim+(0)*yDim+iRow2];
        double below = iRow1 == 0 ? above : in[(iRow1 - 1)*yDim*zDim+(0)*yDim+iRow2];
        double front = iRow2 == max2 - 1 ? in[(iRow1)*yDim*zDim+ (0)*yDim+iRow2-1]: in[(iRow1)*yDim*zDim+(0)*yDim+iRow2+1];
        double back = iRow1 == 0 ? above : in[(iRow1)*yDim*zDim+(0)*yDim+iRow2-1];
        double middle=in[(iRow1)*yDim*zDim+0*yDim+iRow2];
        double di=db*middle + above*dac + below*dac + front*dac + back*dac;
        scratch[len] = (di) / b;

        //Doing the forward passes
        for (int i = 1; i < len-1; i++) {
            scratch[i] = ac / (b - ac*scratch[i - 1]);
        }
        for (int i = 1; i < len; i++) {

            above = iRow1 == max1 - 1 ? in[(iRow1 - 1)*yDim*zDim+ (i)*yDim+iRow2]: in[(iRow1+1)*yDim*zDim+(i)*yDim+iRow2];
            below = iRow1 == 0 ? above : in[(iRow1 - 1)*yDim*zDim+(i)*yDim+iRow2];
            front = iRow2 == max2 - 1 ? in[(iRow1)*yDim*zDim+ (i)*yDim+iRow2-1]: in[(iRow1)*yDim*zDim+(i)*yDim+iRow2+1];
            back = iRow1 == 0 ? above : in[(iRow1)*yDim*zDim+(i)*yDim+iRow2-1];
            middle=in[(iRow1)*yDim*zDim+i*yDim+iRow2];
            di=db*middle + above*dac + below*dac + front*dac + back*dac;
            if(i<len-1){scratch[len + i] = (di - ac*scratch[len + i - 1]) / (b - ac*scratch[i - 1]);}
            else{scratch[len + i] = (di - 2*ac*scratch[len + i - 1]) / (b - 2*ac*scratch[i - 1]);}
        }

        //backward pass, do not touch!
        out[(iRow1)*yDim*zDim+(len - 1)*yDim+iRow2] = scratch[len * 2 - 1];
        for (int i = len - 2; i >= 0; i--) {
            out[(iRow1)*yDim*zDim+(i)*yDim+iRow2] = scratch[len + i] - scratch[i] * out[(iRow1)*yDim*zDim+(i + 1)*yDim+iRow2];
        }
    }

    //only changes from 2D to 3D are 4x diffRate for db, front and back added to d side of equation, and array accesses are different
    public static void TDMA3z(final double[] in, final double[] out,final double[] scratch,final int xDim, final int yDim,final int zDim, final int iRow1,final int iRow2, final double diffRate) {
        //doing:z, iRow1:x, iRow2:y
        final double ac=diffRate;
        final double b=-2*(2*diffRate+1);
        final double db=2*(2*diffRate-1);
        final double dac=-diffRate;

        final int max1 = xDim;
        final int max2 = yDim;
        final int len =  zDim;

        //Doing the 0 entries
        scratch[0] = (2.0*ac)/b;
        double above = iRow1 == max1 - 1 ? in[(iRow1 - 1)*yDim*zDim+ (iRow2)*yDim+0]: in[(iRow1+1)*yDim*zDim+(iRow2)*yDim+0];
        double below = iRow1 == 0 ? above : in[(iRow1 - 1)*yDim*zDim+(iRow2)*yDim+0];
        double front = iRow2 == max2 - 1 ? in[(iRow1)*yDim*zDim+ (iRow2-1)*yDim+0]: in[(iRow1)*yDim*zDim+(iRow2+1)*yDim+0];
        double back = iRow1 == 0 ? above : in[(iRow1)*yDim*zDim+(iRow2-1)*yDim+0];
        double middle=in[(iRow1)*yDim*zDim+iRow2*yDim+0];
        double di=db*middle + above*dac + below*dac + front*dac + back*dac;
        scratch[len] = (di) / b;

        //Doing the forward passes
        for (int i = 1; i < len-1; i++) {
            scratch[i] = ac / (b - ac*scratch[i - 1]);
        }
        for (int i = 1; i < len; i++) {

            above = iRow1 == max1 - 1 ? in[(iRow1 - 1)*yDim*zDim+ (iRow2)*yDim+i]: in[(iRow1+1)*yDim*zDim+(iRow2)*yDim+i];
            below = iRow1 == 0 ? above : in[(iRow1 - 1)*yDim*zDim+(iRow2)*yDim+i];
            front = iRow2 == max2 - 1 ? in[(iRow1)*yDim*zDim+ (iRow2-1)*yDim+i]: in[(iRow1)*yDim*zDim+(iRow2+1)*yDim+i];
            back = iRow1 == 0 ? above : in[(iRow1)*yDim*zDim+(iRow2-1)*yDim+i];
            middle=in[(iRow1)*yDim*zDim+iRow2*yDim+i];
            di=db*middle + above*dac + below*dac + front*dac + back*dac;
            if(i<len-1){scratch[len + i] = (di - ac*scratch[len + i - 1]) / (b - ac*scratch[i - 1]);}
            else{scratch[len + i] = (di - 2*ac*scratch[len + i - 1]) / (b - 2*ac*scratch[i - 1]);}
        }

        //backward pass, do not touch!
        out[(iRow1)*yDim*zDim+(iRow2)*yDim+len - 1] = scratch[len * 2 - 1];
        for (int i = len - 2; i >= 0; i--) {
            out[(iRow1)*yDim*zDim+(iRow2)*yDim+i] = scratch[len + i] - scratch[i] * out[(iRow1)*yDim*zDim+(iRow2)*yDim+i + 1];
        }
    }
    public static void DiffusionADI3(int iAxis,double[]inGrid,double[]outGrid,final double[] scratch,final int xDim,final int yDim,final int zDim,final double diffRate){
        switch (iAxis){
            case 0://x axis case
                for (int y = 0; y < yDim; y++) {
                    for (int z = 0; z < zDim; z++) {
                        TDMA3x(inGrid,outGrid,scratch,xDim,yDim,zDim,y,z,diffRate);
                    }
                }
                break;
            case 1://y axis case
                for (int x = 0; x < yDim; x++) {
                    for (int z = 0; z < zDim; z++) {
                        TDMA3y(inGrid,outGrid,scratch,xDim,yDim,zDim,x,z,diffRate);
                    }
                }
                break;
            case 2://z axis case
                for (int x = 0; x < xDim; x++) {
                    for (int y = 0; y < zDim; y++) {
                        TDMA3z(inGrid,outGrid,scratch,xDim,yDim,zDim,x,y,diffRate);
                    }
                }
                break;
            default:throw new IllegalArgumentException("iAxis variable must be one of 0(x),1(y),2(z)");
        }
    }

    public static void DiffusionADI2(boolean xAxis,final double[]inGrid, final double[]outGrid,final double[]scratch,final int xDim,final int yDim,final double diffRate) {
        int len=xAxis?yDim:xDim;
        if(xAxis){
            for (int i = 0; i < len; i++) {
                TDMAx(inGrid, outGrid, scratch, xDim, yDim, i, diffRate);
            }
        }
        else {
            for (int i = 0; i < len; i++) {
                TDMAy(inGrid, outGrid, scratch, xDim, yDim, i, diffRate);
            }
        }
    }
}



