package meshi.molecularElements.atoms;
import meshi.parameters.*;
/**
 * A shortcut into the atom's guts.
 **/
public class AtomCore {
    protected    AtomStatus status;
    public final AtomStatus status() {return status;}

    protected double[] x,y,z;
    public final void addForce(double fx, double fy, double fz) {
        x[1] += fx;
        y[1] += fy;
        z[1] += fz;
    }
    public final double x() {
        return x[0];
    }
    public final double y() {return y[0];}
    public final double z() {return z[0];}

    protected    AtomType type;
    protected    AtomType samudralaType;
    public final AtomType type() {return type;}

    public final Atom     atom;
    public final int      number;
    protected AtomCore(Atom atom, AtomType type, AtomStatus status, int number, double x, double y, double z) {
            this.atom = atom;
        this.type   = type;
        this.samudralaType   = type;
        this.status = status;
        this.number = number;
        this.x = new double[2];
        this.y = new double[2];
        this.z = new double[2];
        this.x[0]   = x;
        this.y[0]   = y;
        this.z[0]   = z;
        this.x[1]   = 0;
        this.y[1]   = 0;
        this.z[1]   = 0;
    }

    public String toString() {
        return "AtomCore ("+number+" "+status+" "+type+" "+x[0]+" "+x[1]+" ; "+y[0]+" "+y[1]+" ; "+z[0]+" "+z[1]+") ";
    }
}
  
