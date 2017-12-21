package mapreduce.cfsdp;


import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;

/**
 * Created by Chengheng on 2016/10/21.
 */
public class Grid {
    static final int leftShif=10;
    static final double[] xparts = new double[]{-20, -15,-10,-5,0,5,10, 20, 25,30};
    static final double[] yparts = new double[]{-10,-5, 0, 5,10,15,20,25};
    static final double xBorderWidth = 1;
    static final double yBorderWidth = 1;

    public static void main(String[] args) {
        for (int i = 0; i < 100000; i++) {
            List<Long> l=getSendTo(36,15.2);
            for (Long integer : l) {
                System.out.println(integer/16+" "+integer%1000);
            }
        }
    }

    static private int getSection(double[] a, double key) {
        int low = 0;
        int high = a.length - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            double midVal = a[mid];

            if (midVal < key) low = mid + 1;  // Neither val is NaN, thisVal is smaller
            else if (midVal > key) high = mid - 1; // Neither val is NaN, thisVal is larger
            else {
                long midBits = Double.doubleToLongBits(midVal);
                long keyBits = Double.doubleToLongBits(key);
                if (midBits == keyBits)     // Values are equal
                    return mid;             // Key found
                else if (midBits < keyBits) // (-0.0, 0.0) or (!NaN, NaN)
                    low = mid + 1;
                else                        // (0.0, -0.0) or (NaN, !NaN)
                    high = mid - 1;
            }
        }
        return low;  // key not found.
    }


    public static List<Long> getSendTo(double x, double y) {
        List<Long> sendTo = Lists.newArrayList();
        int gx = getSection(xparts, x);
        int gy = getSection(yparts, y);
        sendTo.add(((long)gx << leftShif) + gy);
        double xleft = gx > 0 ? xparts[gx - 1] : Double.MIN_VALUE;
        double xright = gx < xparts.length ? xparts[gx ] : Double.MAX_VALUE;
        double ydown = gy > 0 ? yparts[gy - 1] : Double.MIN_VALUE;
        double yup = gy < yparts.length ? yparts[gy] : Double.MAX_VALUE;


        if (x - xleft >= xBorderWidth && xright - x >= xBorderWidth) {//1
            if (y - ydown >= yBorderWidth && yup - y >= yBorderWidth) {
                //1a donothing
            } else if (y - ydown < yBorderWidth) {
                //1b
                if (gy > 0) {
                    sendTo.add(((long)gx << leftShif) + gy - 1);
                }
            } else {
                //1c
                if (gy < yparts.length) {
                    sendTo.add(((long)gx << leftShif) + gy + 1);
                }
            }
        } else if (x - xleft < xBorderWidth) {//2
            if (y - ydown >= yBorderWidth && yup - y >= yBorderWidth) {
                //2a
                if (gx > 0) {
                    sendTo.add((((long)gx - 1) << leftShif) + gy);
                }
            } else if (y - ydown < yBorderWidth) {
                //2b
                if (gx > 0) {
                    sendTo.add((((long)gx - 1) << leftShif) + gy);
                }
                if (gy > 0) {
                    sendTo.add(((long)gx << leftShif) + gy - 1);
                }
                if (gx > 0 && gy > 0) {
                    sendTo.add((((long)gx - 1) << leftShif) + gy - 1);
                }
            } else {
                //2c
                if (gx > 0) {
                    sendTo.add((((long)gx - 1) << leftShif) + gy);
                }
                if (gy < yparts.length) {
                    sendTo.add(((long)gx << leftShif) + gy + 1);
                }
                if (gx > 0 && gy < yparts.length) {
                    sendTo.add((((long)gx - 1) << leftShif) + gy + 1);
                }

            }
        } else {
            //3
            if (y - ydown >= yBorderWidth && yup - y >= yBorderWidth) {
                //3a
                if (gx < xparts.length) {
                    sendTo.add((((long)gx + 1) << leftShif) + gy);
                }
            } else if (y - ydown < yBorderWidth) {
                //3b
                if (gx < xparts.length) {
                    sendTo.add((((long)gx + 1) << leftShif) + gy);
                }
                if (gy > 0) {
                    sendTo.add(((long)gx << leftShif) + gy - 1);
                }
                if (gx < xparts.length && gy > 0) {
                    sendTo.add((((long)gx + 1) << leftShif) + gy - 1);
                }
            } else {
                //3c
                if (gx < xparts.length) {
                    sendTo.add((((long)gx + 1) << leftShif) + gy);
                }
                if (gy < yparts.length) {
                    sendTo.add(((long)gx << leftShif) + gy + 1);
                }
                if (gx < xparts.length && gy < yparts.length) {
                    sendTo.add((((long)gx + 1) << leftShif) + gy + 1);
                }
            }
        }

        return sendTo;
    }
}
