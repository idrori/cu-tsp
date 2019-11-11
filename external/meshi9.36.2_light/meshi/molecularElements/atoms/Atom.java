/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

import meshi.PDB.PdbLine;
import meshi.geometry.Coordinates;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.parameters.AtomType;
import meshi.parameters.BBatom;
import meshi.util.*;
import meshi.util.string.StringList;
import meshi.util.string.StringParser;

import java.util.HashMap;

/**
 * A generic atom.
 * <br><hr>
 * <b>A possible pitfall</b> <br>
 * A {@link meshi.geometry.Coordinates Coordinates} object
 * ({@link meshi.molecularElements.atoms.Atom#core#coordinates coordinates})is associated with an atom
 * and handle all geometric issues(position, distance, forces etc.).
 * The "geometric" methods ({@link meshi.molecularElements.atoms.Atom#x x()},
 * {@link meshi.molecularElements.atoms.Atom#y y()},
 * {@link meshi.molecularElements.atoms.Atom#distanceFrom distanceFrom(Atom)}, etc.)
 * are used to encapsulate this implementation. The {@link meshi.geometry.Coordinates Coordinates}
 * object itself though, is publicly accessible. Removing the overhead of function calls accelerates
 * {@link meshi.geometry.DistanceMatrix DistanceMatrix} updating which is the bottleneck of
 * energy based calculations. <b> It may though, open the way to spectacular bugs!!!</b>.
 *
 * @see AtomList
 */
public class Atom implements Comparable, Attributable {
    //------------------------------------------------------------------------------
    //---------------------- class and object variables ----------------------------
    //------------------------------------------------------------------------------

    public final AtomCore core;

    /**
     * The getAtoms name. In proteins a unique specification of it's position within the residue.
     */
    public final String name;

    /**
     * A unique identifier of the atom in the MolecularSystem.
     */
    public final int ID;
    public final MolecularSystem molecularSystem;


    /**
     * A unique identifier of the atom in the molecular system.
     * Among other things serves as an index into a {@link meshi.geometry.DistanceMatrix}.
     */
    public final int number() {
        return core.number;
    }


    private double energy = 0;

    /**
     * Atom's Residue.
     */
    protected Residue residue;

    private Double temperatureFactor = null;

    /**
     * The atom type.
     */
    public AtomType type() {
        return core.type;
    }

    public void setType(AtomType newType) {
        if (((type().backboneN()) & (newType == AtomType.TRN)) ||
                ((type().backboneC()) & (newType == AtomType.TRC)) ||
                ((type().backboneO()) & (newType == AtomType.TRO)))
            core.type = newType;
        else throw new RuntimeException("Cannot setResidue atom " + this + " of type " + core.type + " to " + newType);
    }


    /**
     * Atoms bonded to this atom.
     */
    private final int ATOM_BONDED_CAPACITY = 4;
    protected AtomList bonded;


    /**
     * Number of getAtoms in all the molecular systems.
     */
    private static int numberOfAtoms = 0;


    private PdbLine pdbLine;

    /* private double reliability = 0;
        public double reliability() {return reliability;}
        public void setReliability(double reliability) {this.reliability = reliability;}
    */   
     public AtomType samudralaType() {return core.samudralaType;}

    //------------------------------------------------------------------------------
    //----------------------------- constructors -----------------------------------
    //------------------------------------------------------------------------------

    /**
     * A generic atom with unspecified coordinates.
     * Mainly for use by other more specific constructors.
     */
    public Atom(String name, Residue residue, AtomType type, Coordinates coordinates, double occupency, Double temperatureFactor,MolecularSystem molecularSystem) {
        AtomStatus status;
        this.molecularSystem = molecularSystem;
        bonded = new AtomList(ATOM_BONDED_CAPACITY,molecularSystem);
        if (coordinates.nowhere()) status = AtomStatus.NOWHERE;
        else status = AtomStatus.NORMAL;
        core = molecularSystem.createAtomCore(this, type, status, coordinates.x(), coordinates.y(), coordinates.z());

        this.name = name;
        ID = numberOfAtoms();
        numberOfAtoms++;
        this.residue = residue;
        this.temperatureFactor = temperatureFactor;
        if (residue != null) {
            pdbLine = new PdbLine(number(), name,
                    "", // default alternative location
                    residue.name,
                    residue.chain(), residue.number(),
                    coordinates.x(), coordinates.y(), coordinates.z(),
                    occupency, temperatureFactor.doubleValue()); // default occupency
        } else pdbLine = new PdbLine(number(), name,
                "", // default alternative location
                "UNK",
                " ", -99,
                coordinates.x(), coordinates.y(), coordinates.z(),
                temperatureFactor.doubleValue(), 1); // default occupency and temperature factor
    }



        public Atom(PdbLine line,MolecularSystem molecularSystem) {
        this(line.name(),
                new Residue(new ResidueIdentifier(line.chain(),
                        line.residueNumber(),-1),  //An atom can not know the chain number
                        line.residueName()), // residue needs to be assigned
                        line.type(), new Coordinates(line.x(),
                        line.y(), line.z()), line.occupancy(), new Double(line.temperatureFactor()),molecularSystem);
    }


    //------------------------------------------------------------------------------
    //------------------------------- methods --------------------------------------
    //------------------------------------------------------------------------------

    public boolean nowhere() {
        return (core.status == AtomStatus.NOWHERE);
    }

    public boolean active() {
        return core.status.active();
    }

    public PdbLine pdbLine() {
        return pdbLine;
    }

    /**
     * Returns the Atom's chain.
     *
     * @see meshi.molecularElements.atoms.Atom#chain
     */
    public String chain() {
        return residue.chain();
    }

    /**
     * Set the atom residue.
     *
     * @see meshi.molecularElements.atoms.Atom#residue
     */
    public void setResidue(Residue residue) {
        this.residue = residue;
        pdbLine = new PdbLine(number(), name,
                "", // default alternative location
                residue.name,
                residue.chain(), residue.number(),
                core.x[0], core.y[0], core.z[0],
                1, 1); // default occupency and temperature factor

    }

    /**
     * Returns the atom's residue.
     *
     * @see meshi.molecularElements.atoms.Atom#residue
     */
    public final Residue residue() {
        return residue;
    }


    /**
     * Returns true if the protein's name is a substring of the String parameter seperated by spaces.
     * Examples:
     * <ol><li>
     * String s = "CA CB";            <br>
     * Atom a;                          <br>
     * a.name = "CB";                   <br>
     * System.out.println(a.nameIs(s)); <br>
     * <b> prints true </b>             <br>
     * </li>
     * <li>
     * String s = "CACB";             <br>
     * Atom a;                          <br>
     * a.name = "CB";                   <br>
     * System.out.println(a.nameIs(s)); <br>
     * <b> prints false </b>            <br>
     * </li>
     * </ol>
     *
     * @see meshi.molecularElements.atoms.Atom#name
     */
    public boolean nameIs(String keys) {
        return nameIs(StringParser.standard(keys));
    }

    /**
     * Returns true if the protein's name is included in the StringList parameter.
     */
    public boolean nameIs(StringList keys) {
        //Iterator iter = keys.iterator();
        //String string;
        //while ((string = (String) iter.next()) != null)
        for (String string : keys)
            if (string.equals(name)) return true;
        return false;
    }

    /**
     * Returns the X coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#x
     */
    public final double x() {
        if (nowhere()) throw new RuntimeException("This Atom:\n" +
                this + "\n" +
                "is NOWHERE and thus, its position is meaningless.");
        return core.x[0];
    }

    public final double[] X() {
        return core.x;
    }

    /**
     * Returns the Y coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#Y
     */
    public final double y() {
        return core.y[0];
    }

    public final double[] Y() {
        return core.y;
    }

    /**
     * Returns the Z coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#Z
     */
    public final double z() {
        return core.z[0];
    }

    public final double[] Z() {
        return core.z;
    }


    /**
     * Adds the parameter to the X coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#addToX
     */
    public void addToX(double addMe) {
        if (nowhere()) throw new RuntimeException("This Atom:\n" +
                this + "\n" +
                "is NOWHERE and thus, it makes no sence to move it relative to the\n" +
                "current position.");
    //    if (addMe > 1000) addMe = 1000;
        core.x[0] += addMe;
    }

    /**
     * Adds the parameter to the Y coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#addToY
     */
    public void addToY(double addMe) {
        if (addMe > 1000) addMe = 1000;
        core.y[0] += addMe;
    }

    /**
     * Adds the parameter to the Y coordinate of the atom.
     *
     * @see meshi.geometry.Coordinates#addToZ
     */
    public void addToZ(double addMe) {
        if (addMe > 1000) addMe = 1000;
        core.z[0] += addMe;
    }

    /**
     * Sets the X,Y,Z coordinate of the atom.
     */
    public void setXYZ(double newx, double newy, double newz) {
        //if (!nowhere()) throw new RuntimeException("Only the coordinates of nowhere getAtoms can be setResidue\n"+this+" "+"\n"+core+"\n");
        core.x[0] = newx;
        core.y[0] = newy;
        core.z[0] = newz;
        core.status = AtomStatus.NORMAL;
    }

    public void setXYZ(double newx, double newy, double newz, AtomStatus status) {
        //if (!nowhere()) throw new RuntimeException("Only the coordinates of nowhere getAtoms can be setResidue\n"+this+" "+"\n"+core+"\n");
         core.x[0] = newx;
        core.y[0] = newy;
        core.z[0] = newz;
        core.status = status;
    }

    public void setXYZ(Coordinates coor) {
        //if (!nowhere()) throw new RuntimeException("Only the coordinates of nowhere getAtoms can be setResidue\n"+this+" "+"\n"+core+"\n");
        core.x[0] = coor.x();
        core.y[0] = coor.y();
        core.z[0] = coor.z();
        if (coor.nowhere()) core.status = AtomStatus.NOWHERE;
        else core.status = AtomStatus.NORMAL;
    }

    public void resetCoordinates() {
        core.x[0] = Coordinates.NOWHERE_CONST;
        core.y[0] = Coordinates.NOWHERE_CONST;
        core.z[0] = Coordinates.NOWHERE_CONST;
        core.status = AtomStatus.NOWHERE;
    }

    /**
     * Returns the force operating isOn the atom in the X direction.
     *
     * @see meshi.geometry.Coordinates#fx
     */
    public double fx() {
        if (nowhere()) throw new RuntimeException("This Atom:\n" +
                this + "\n" +
                "is NOWHERE and thus, the force isOn it is meaningless.");
        return core.x[1];
    }

    /**
     * Returns the force operating isOn the atom in the Y direction.
     *
     * @see meshi.geometry.Coordinates#fy
     */
    public double fy() {
        return core.y[1];
    }

    /**
     * Returns the force operating isOn the atom in the Z direction.
     *
     * @see meshi.geometry.Coordinates#fz
     */
    public double fz() {
        return core.z[1];
    }

    /**
     * Sets the force operating isOn the atom in the X direction.
     *
     * @see meshi.geometry.Coordinates#setFx
     */
    public void setFx(double fx) {
        core.x[1] = fx;
    }

    /**
     * Sets the force operating isOn the atom in the Y direction.
     *
     * @see meshi.geometry.Coordinates#setFy
     */
    public void setFy(double fy) {
        core.y[1] = fy;
    }

    /**
     * Sets the force operating isOn the atom in the Z direction.
     *
     * @see meshi.geometry.Coordinates#setFz
     */
    public void setFz(double fz) {
        core.z[1] = fz;
    }

    /**
     * Adds the parameter to the force operating isOn the atom in the X direction.
     *
     * @see meshi.geometry.Coordinates#addToFx
     */
    public final void addToFx(double addMe) {
        if (nowhere()) throw new RuntimeException("This Atom:\n" +
                this + "\n" +
                "is NOWHERE it must not take part in energy calculations.\n" +
                "The forces isOn it. Which you are now trying to modify are \n" +
                "meaningless.");
        core.x[1] += addMe;
    }

    /**
     * Adds the parameter to the force operating isOn the atom in the Y direction.
     *
     * @see meshi.geometry.Coordinates#addToFy
     */
    public final void addToFy(double addMe) {
        core.y[1] += addMe;
    }

    /**
     * Adds the parameter to the force operating isOn the atom in the Z direction.
     *
     * @see meshi.geometry.Coordinates#addToFz
     */
    public final void addToFz(double addMe) {
        core.z[1] += addMe;
    }


    /**
     * The distance between this atom and the parameter.
     */
    public final double distanceFrom(Atom atom) {
        try {
            double dx = x() - atom.x();
            double dy = y() - atom.y();
            double dz = z() - atom.z();
            return Math.sqrt(dx * dx + dy * dy + dz * dz);
        } catch (RuntimeException ex) {
            throw new RuntimeException("Cannot measure the distance between\n" + this + "\nand\n" + atom + "\n" + ex);
        }
    }

    /**
     * Returns the atom as a PDB formatted String.
     */
    public String verbose(int level) {
        String frozenString;

        if (frozen()) frozenString = "frozen";
        else frozenString = "notFrozen";

        return toString() + " " + type() + " " + frozenString;
    }

    public String toString() {
        String residueName;
        String chain;
        int residueNumber;
        if (residue != null) {
            residueName = residue.name;
            chain = residue.chain();
            residueNumber = residue.number();
        } else {
            residueName = "UNK";
            chain = " ";
            residueNumber = -99;
        }
        return (new PdbLine(number(),
                name,
                "", // default alternative location
                residueName,
                chain, residueNumber,
                core.x[0], core.y[0], core.z[0],
                1, temperatureFactor) // default occupency and temperature factor
        ).toString();
    }

    /**
     * Move the atom to a random position within <b>radius</b> from (<b>centerx</b>,<b>centery</b>,<b>centerz</b>).
     */
    public void randomize(double radius,
                          double centerx, double centery, double centerz) {
        setXYZ(centerx + (MeshiProgram.randomNumberGenerator().nextDouble() - 0.5) * radius,
                centery + (MeshiProgram.randomNumberGenerator().nextDouble() - 0.5) * radius,
                centerz + (MeshiProgram.randomNumberGenerator().nextDouble() - 0.5) * radius);
    }

    public void randomize(double radius, Atom center) {
        randomize(radius, center.x(), center.y(), center.z());
    }

    /**
     * Returns atom's comment.
     */
    public String comment() {
        return "Atom";
    }

    /**
     * Returns the list of getAtoms bonded to this one.
     *
     * @see meshi.molecularElements.atoms.Atom#bonded
     */
    public AtomList bonded() {
        return bonded;
    }

    /**
     * Bonds the <b>other</b> atom  to this one.
     * Eeach atom is added to the bonded list of the other.
     */
    public AtomPair bond(Atom other) {
        bonded.add(other);
        other.bonded.add(this);
        return new AtomPair(other, this);
    }

    public final int residueNumber() {
        return residue.number();
    }

    public String name() {
        return name;
    }

    public static int numberOfAtoms() {
        return numberOfAtoms;
    }

    public double occupancy() {
        return pdbLine.occupancy();
    }

    public Double temperatureFactor() {
        return temperatureFactor;
    }

    public void setTemperatureFactor(Double tf) {
        if (tf > 99.99) throw new RuntimeException("This is weird. tf = " + tf);
        temperatureFactor = tf;
    }

    public String alternateLocation() {
        return pdbLine.alternateLocation();
    }

    public String residueName() {
        if (residue != null) return residue.name;
        return pdbLine.residueName();
    }


    public void freeze() {
        freeze(true);
    }
    public void freeze(boolean distanceMatrixCreateFlag) {
        if (core.status == AtomStatus.NOWHERE) {
            throw new RuntimeException("Cannot freeze atom of status " + core.status);
        }
        core.status = AtomStatus.FROZEN;
        if (distanceMatrixCreateFlag) molecularSystem.createDistanceMatrix("Atom "+this+" froze. You must instantiate a new TotalEnergy.");
    }

    public void hide() {
        if (core.status == AtomStatus.NOWHERE) {
            throw new RuntimeException("Cannot hide atom of status " + core.status);
        }
        core.status = AtomStatus.HIDDEN;
        molecularSystem.createDistanceMatrix("Atom "+this+" hided. You must instantiate a new TotalEnergy.");
   }

    public boolean frozen() {
        return (core.status.frozen());
    }

    public boolean normal() {
        return (core.status == AtomStatus.NORMAL);
    }

    /**
     * An insecure defrost  method. It allows the simulation to continue even after some getAtoms have melted.
     * If this happens during a minimization, the minimization would become severely instable.
     */
     public void defrost() {
         defrost(true);
     }
    public void defrost(boolean newDistanceMatrixFlag) {
        if (!frozen()) {
            throw new RuntimeException("Cannot defrost an atom of status " + core.status);
        }
        core.status = AtomStatus.NORMAL;
        core.x[1] = 0;
        core.y[1] = 0;
        core.z[1] = 0;
        if (newDistanceMatrixFlag) molecularSystem.createDistanceMatrix("Atom "+this+" melted. You must instantiate a new TotalEnergy.");
    }

    public int compareTo(Object obj) {
        Atom other = (Atom) obj;
        if (number() > other.number()) return 1;
        if (number() < other.number()) return -1;
        return 0;
    }


    public void resetEnergy() {
        energy = 0;
    }

    public void resetTemperatureFactor() {
        temperatureFactor = 0.0;
    }

    public void addEnergy(double add) {
        energy += add;
    }

    public double energy() {
        return energy;
    }

    public final boolean isBackbone() {
        return type().backbone();
    }

    public final boolean backboneH() {
        return type().backboneH();
    }

    public final boolean backboneN() {
        return type().backboneN();
    }

    public final boolean backboneCA() {
        return type().backboneCA();
    }

    public final boolean backboneC() {
        return type().backboneC();
    }

    public final boolean backboneO() {
        return type().backboneO();
    }

    public final boolean isCarbon() {
        return type().isCarbon();
    }

    public final boolean isOxygen() {
        return type().isOxygen();
    }

    public final boolean isNitrogen() {
        return type().isNitrogen();
    }

    public final boolean isSulfur() {
        return type().isSulfur();
    }

    public final boolean isHydrogen() {
        return type().isHydrogen();
    }

    public final BBatom bbAtom() {
        return type().bbAtom();
    }

    public void setStatus(AtomStatus status) {
        core.status = status;
    }

//     public Coordinates coordinates() {
// 	if (nowhere()) return new Coordinates();
// 	return new Coordinates(x(), y(), z());
//     }

    //----------------------------------- Attributes ---------------------------------------
    private HashMap attributes = new HashMap();

    public final void addAttribute(MeshiAttribute attribute) {
        attributes.put(attribute.key(), attribute);
    }

    public final MeshiAttribute getAttribute(int key) {
        return (MeshiAttribute) attributes.get(key);
    }

    // --------------------------------- for homology modeling ------------------------------


    public void setResidueNumber(int number) {
        throw new MeshiException("Chen and Oren: " +
                "Unsure what this method should do.");
    }

    public void emptyBonded() {
        bonded = new AtomList(ATOM_BONDED_CAPACITY,molecularSystem);
    }

    public void setNowhere() {
        core.status= AtomStatus.NOWHERE;

    }

}

