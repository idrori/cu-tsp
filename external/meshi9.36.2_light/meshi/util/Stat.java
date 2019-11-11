package meshi.util;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 31/12/11
 * Time: 23:25
 * To change this template use File | Settings | File Templates.
 */
public class Stat implements Comparable{
    private double sum,sum2;
    private int numberOfObservations;
    private Histogram histogram;
    private double min,max;
    private double score;
    public Stat() {
        reset();
    }
    public Stat(double minHistogram, double maxHistogram, double binHistogram) {
        this();
        histogram = new Histogram(minHistogram, maxHistogram, binHistogram);
    }

    public void reset() {
        sum = 0;
        sum2 = 0;
        numberOfObservations = 0;
        histogram = null;
        min = Double.MAX_VALUE;
        max = Double.MIN_VALUE;
        score = 0;
    }


    public double getScore() {return score;}
    public void setScore(double score) {this.score = score;}
    public int compareTo(Object o) {
        Stat other = (Stat)o;
        double diff = score - other.score;
        if (diff > 0) return 1;
        if (diff == 0) return 0;
        if (diff < 0) return -1;
        throw new RuntimeException("This is weird "+score+" "+other.score+" "+diff+"\n"+this+"\n"+other);
    }
    public void add(double addMe) {
        sum +=  addMe;
        sum2 += addMe*addMe;
        if (histogram != null) histogram.add(addMe);
        numberOfObservations++;
        if (addMe < min ) min = addMe;
        if (addMe > max ) max = addMe;
    }
    public double getMax() {return max;}
    public double getMin() {return min;}
    public double getMean() {
        return sum/numberOfObservations;
    }
    public Double getMedian() {
        if (histogram == null) return null;

        return new Double(histogram.median());
    }
    public double getStd() {
        if (numberOfObservations == 0) return -1;
        double mn = getMean();
        return Math.sqrt(sum2/numberOfObservations-mn*mn);
    }
    public double getRms() {
       if (numberOfObservations == 0) return -1;
        return Math.sqrt(sum2/numberOfObservations);
    }
    public    int numberOfObservations() {return numberOfObservations; }
    public double zScore(double d) {
        return (d-getMean())/getStd();
    }
    public Histogram getHistogram() {return histogram;}
}

