package meshi.molecularElements;


public abstract class Conservation {
    public final double conservationScore;
    public Conservation(double conservationScore) {
        this.conservationScore = conservationScore;
    }
    public abstract double conservationWeight();
}
