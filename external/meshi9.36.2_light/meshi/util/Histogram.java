package meshi.util;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 12/01/12
 * Time: 22:31
 * To change this template use File | Settings | File Templates.
 */
public class Histogram {
    private int[] hist;
    private double min,max,bin;
    private int nBins;
    private int nObservations, tooLow, tooHigh;

    public Histogram(double min, double max, double bin) {
        double n = 1+(max-min)/bin;
        nBins = (int) Math.round(n);
        if (nBins != n) throw new RuntimeException("bad min max bin "+min+" "+max+" "+bin);
        hist = new int[nBins];
        this.max = max;
        this.min = min;
        tooHigh = 0;
        tooLow = 0;
        this.bin = bin;
    }
    public void add(double addMe) {
        if (addMe < min) tooLow++;
        else if (addMe > max) tooHigh++;
        else {
            nObservations++;
            int iBin = (int) Math.floor((addMe-min)/bin);
            if (iBin >= hist.length)
                throw new RuntimeException("hist.length = "+hist.length+" iBin = "+iBin+"\n"+
                                           "addMe = "+addMe+" min = "+min+" bin = "+bin);
            hist[iBin]++;
        }
    }
     public double median() {
            double sum = 0;
            for (int i = 0; i < hist.length; i++) {
                  sum += 1.0*hist[i]/nObservations;
                  if (sum >= 0.5) return min+bin*i+bin/2;
            }
            throw new RuntimeException("This is weird");
     }

    public int tooLow() {return tooLow;}
    public int tooHigh() {return tooHigh;}
    public int nObservations() {return nObservations;}
    public void print() {
        for (int i = 0; i < hist.length; i++) {
            System.out.println(""+min+" - "+(i+0.5)*bin+" "+hist[i]);
        }
    }
    public double[] fractionalHistogram() {
        double[] out = new double[hist.length];
        for (int i = 0; i < hist.length; i++)
            out[i] = (hist[i]*1.0)/nObservations;
        return out;
    }
    public int getBin(int iBin) {
        return hist[iBin];
    }
}

