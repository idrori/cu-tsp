package meshi.geometry;

import meshi.molecularElements.atoms.Atom;

public class AlphaCoordinates {
    // Consider 4 consecutive CA atoms 0-3
    private double alpha, angle0, angle1;



    boolean valid;

    public AlphaCoordinates(Atom atom0, Atom atom1, Atom atom2, Atom atom3, DistanceMatrix distanceMatrix) {
        Angle angle0 = new Angle(atom0, atom1, atom2, distanceMatrix);
        this.angle0 = angle0.angle();
        Angle angle1 = new Angle(atom1, atom2, atom3, distanceMatrix);
        this.angle1 = angle1.angle();
        alpha = (new Torsion(angle0, angle1, distanceMatrix)).torsion();
        valid = true;
    }
    public String toString() {
        if (valid)
            return "["+Angle.rad2deg(alpha)+", "+Angle.rad2deg(angle0)+", "+Angle.rad2deg(angle1)+"]";
        return "AlphaCoordinates invalid";
    }


    public AlphaCoordinates() {
        valid = false;
    }

    public double getAngle0() {
        return angle0;
    }

    public double getAngle1() {
        return angle1;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setValid(boolean valid) {
        this.valid = valid;
    }

    public boolean isValid() {
        return valid;
    }
}
