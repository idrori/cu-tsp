/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

import meshi.PDB.PdbLine;
import meshi.PDB.PdbLineATOM;
import meshi.PDB.PdbLineIterator;
import meshi.PDB.PdbReader;
import meshi.geometry.Coordinates;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.ResidueList;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueType;
import meshi.sequences.AtomAlignment;
import meshi.util.MeshiException;
import meshi.util.Rms;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

/**
 * A list of Atoms
 */
public class AtomList extends ArrayList<Atom> {
    private static int nUnknown = 0;
    private static final int maxUnknown = 1000;
    private static String[][] unknowns = new String[maxUnknown][2];
    private static String freezeComment = null;

    private MeshiLineReader sourceFile = null;
    public final MolecularSystem molecularSystem;

    public MolecularSystem molecularSystem() {
        return molecularSystem;
    }


    private String comment = null;

    public AtomList(MolecularSystem molecularSystem) {
        super();
        this.molecularSystem = molecularSystem;
    }

    public AtomList(Atom atom) {
        super();
        molecularSystem = atom.molecularSystem;
        add(atom);
    }

    public AtomList(int initialCapaciy, MolecularSystem molecularSystem) {
        super(initialCapaciy);
        this.molecularSystem = molecularSystem;
        this.print();
    }

    public AtomList(AtomList other) {
        molecularSystem = new MolecularSystem();
        for (Atom atom : other) {
            PdbLine pdbLine = atom.pdbLine();
            add(new Atom(pdbLine, molecularSystem));
        }
    }

    public void setNowhere(){
        for( Atom atom : this){
            atom.setNowhere();
        }
    }


    /**
     * An atom list based isOn a PDB file
     * In a multiple model file (such as NMR's) only model 1 is read.
     * Currently whenever there are two alternative atom positions present,
     * only position A is taken.
     */
    public AtomList(MeshiLineReader file, Filter lineFilter) {
        this(new PdbLineIterator((PdbReader) file), lineFilter);
    }

    public AtomList(PdbLineIterator lineIterator, Filter lineFilter) {
        PdbLine line;
        boolean contRead = true;
        molecularSystem = new MolecularSystem();
        sourceFile = lineIterator.sourceFile();
/*        if (!(file instanceof PdbReader))
            throw new MeshiException("Cannot construct AtomList from file: " + file + "\n" +
                    "Currently only PDB formated files are used.");   */
        while (lineIterator.hasNext() & contRead) {
            line = lineIterator.next();
            //if (line.isAModel() & (line.getModel() != 1))      ToDo Take care of NMR multiple models
            //contRead = false;
            if (lineFilter.accept(line)) {
                    if (line.type().equals(AtomType.XXX)) {
                        boolean found = false;
                        for (int i = 0; (i < nUnknown) & (!found); i++) {
                            if (unknowns[i][0].equals(line.residueName()) & unknowns[i][1].equals(line.name()))
                                found = true;
                        }
                        if (!found) {
                            if (nUnknown >= maxUnknown)
                                throw new RuntimeException("this is weird - too many unknown atom types");
                            Utils.println("Cannot identify Atom type for " + line);
                            unknowns[nUnknown][0] = line.residueName();
                            unknowns[nUnknown][1] = line.name();
                            nUnknown++;
                        }

                    } else {
                        try {
                            add(new Atom(line,molecularSystem));
                        }
                        catch (OutOfMemoryError er) {
                            System.out.println("************************\n" + "Cannot Create an atom from line: " + line + "\n" +
                                    "of file " + sourceFile.fileName() + "\n" +
                                    "This often happens when numerous protein structure files are read (say for some statistics and they all use the same Molecular System object.\n" +
                                    "As a result this MS growes even if previously used getAtoms are no longer needed.");
                            throw er;
                        }
                    }

            }
        }
        setComment(sourceFile.name());
    }

    /**
     * An atom list based isOn a PDB file
     */
    public AtomList(PdbReader file) {
        this(file, new PdbLineATOM());
    }

    public AtomList(String dir, String fileName) {
        this(new PdbReader(dir, fileName));
    }

    public AtomList(String fileName, Filter lineFilter) {
        this(new PdbReader(fileName), lineFilter);
    }

    public AtomList(String fileName) {
        this(fileName, new PdbLineATOM());
    }


    public AtomList(String fileName,boolean flag) throws IOException{
        super();
        molecularSystem = new MolecularSystem();
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String line;

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("ATOM")) {
                StringTokenizer tokenizer = new StringTokenizer(line);
                tokenizer.nextToken(); // ATOM
                tokenizer.nextToken(); // atom number

                String name = tokenizer.nextToken();
                String residueName = tokenizer.nextToken();
                int    residueNumber = (new Integer(tokenizer.nextToken())).intValue();
                String xString = tokenizer.nextToken();
                try {
                if (!xString.startsWith("-999.999")){
                    double x = (new Double(xString)).doubleValue();
                    double y = (new Double(tokenizer.nextToken())).doubleValue();
                    double z = (new Double(tokenizer.nextToken())).doubleValue();

                    ResidueType residueType = ResidueType.type(residueName);
                    AtomType    atomType    = AtomType.type(residueName, name);
                    ResidueIdentifier residueIdentifier = new ResidueIdentifier(residueNumber);
                    Residue residue = new Residue(residueIdentifier, residueName);
                    Coordinates coordinates = new Coordinates(x,y,z);
                    Atom atom = new Atom(name,residue,atomType,coordinates,1, 0.0,molecularSystem);
                    add(atom);
                }    }
                catch (java.lang.NumberFormatException ex) {
                    System.out.println(ex);
                }
            }
        }
    }


    /**
     * An atom list based on a residue list
     */
    public AtomList(ResidueList residues, MolecularSystem molecularSystem) {
        this(molecularSystem);
        for (Residue residue : residues) {
            if (residue.getAtoms() != null) {
                for (Atom atom : residue.getAtoms())
                    add(atom);
            }
            File tmp = residues.getSourceFile();
            if (tmp != null) {
                String tmpPath = tmp.getAbsolutePath();
                sourceFile = new MeshiLineReader(tmpPath);
            }
        }
    }

    /**
     * Extracts the list elements that pass the filter.
     */
    public AtomList filter(Filter filter) {
        AtomList out = new AtomList(molecularSystem);
        for (Atom atom : this)
            if (filter.accept(atom)) out.add(atom);
        return out;
    }


    public boolean add(Atom atom) {
            if (atom.molecularSystem != molecularSystem)
                throw new RuntimeException("An attempt to add atom:\n" +
                        atom + "\n" +
                        "of " + atom.molecularSystem + " to an atom list of " + molecularSystem);
        //if ((size() > 0) && (atom.type().equals(get(size()-1).type()))) return false;
        return (super.add(atom));
    }

    /**
     * RMS deviation between this list of atom and some other list.
     */
    public Rms rms(AtomList otherList) {
        return new Rms(new AtomAlignment(this, otherList));
    }


    public double getRms(AtomList otherList) {
        return rms(otherList).getRms();
    }

    /**
     * Returns the atom in position <b>i</b> of the list.
     */
    public Atom atomAt(int i) {
        return (Atom) get(i);
    }

    public static class IsAtom implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof Atom);
        }
    }

    public int whereIs(Atom atom) {
        for (int i = 0; i < size(); i++) {
            if (atomAt(i).equals(atom)) return i;
        }
        return -1;
    }

    /**
     * of gyration of the list getAtoms.
     * SQRT{ SIGMA_OVER_ALL_ATOMS([(x-X)^2+(y-Y)^2+(z-Z)^2]/N)}<br>
     * Where (x,y,z) are the atom coordinates and (X,Y,Z) are the
     * coordinates of the center of mass.
     */
    public double radius() {
        return Math.sqrt(radius2());
    }

    public double radius2() {
        double radius;
        double dx, dy, dz, d2;
        double cmx, cmy, cmz; // center of mass x, y and z
        radius = 0;
        cmx = cmy = cmz = 0.0;
        int nWithCoordinates = 0;
        for (Atom atom : this) {
            if ((!atom.nowhere()) &&
                 (atom.x() != -999.999)) {
                cmx += atom.x();
                cmy += atom.y();
                cmz += atom.z();
                nWithCoordinates++;
            }
        }
        cmx /= nWithCoordinates;
        cmy /= nWithCoordinates;
        cmz /= nWithCoordinates;
        for (Atom atom : this) {
            if ((!atom.nowhere()) &&
                 (atom.x() != -999.999)) {
                dx = cmx - atom.x();
                dy = cmy - atom.y();
                dz = cmz - atom.z();
                d2 = dx * dx + dy * dy + dz * dz;
                radius += d2;
            }
        }
        return radius / nWithCoordinates;
    }


    /**
     * Sets the parameter to be the residue of this atom.
     */
    public void setResidue(Residue residue) {
        for (Atom atom : this)
            atom.setResidue(residue);
    }

    public Atom[] toArrayOfAtoms() {
        return toArray(new Atom[size()]);
    }

    public MeshiLineReader sourceFile() {
        return sourceFile;
    }

    public boolean frozenAtomsExist() {
        for (Atom atom : this)
            if (atom.frozen()) return true;
        return false;
    }

    public void freeze(String comment) {
        for (Atom atom : this)
            if (atom.core.status == AtomStatus.NORMAL) atom.freeze(false);
        if (!comment.equals(freezeComment)) {
            molecularSystem.createDistanceMatrix(comment+" - Some getAtoms froze (1). You must instantiate a new TotalEnergy.");
            freezeComment = comment;
        }
        else molecularSystem.createDistanceMatrix();
    }

    public void freeze(Filter filter) {
        for (Atom atom : this) {
            if (filter.accept(atom) & (atom.core.status == AtomStatus.NORMAL)) atom.freeze(false);
        }
        molecularSystem.createDistanceMatrix("Some getAtoms froze (2). You must instantiate a new TotalEnergy.");
    }


    public void defrost() {
        for (Atom atom : this)
            if (atom.frozen()) {
                atom.defrost(false);
            }
        molecularSystem.createDistanceMatrix("Some getAtoms melted. You must instantiate a new TotalEnergy.");
    }

    /*public AtomList fullOccupancyFilter() {
        AtomList out = new AtomList();
        for (Atom atom : this)
            if (atom.occupancy() > 0.99) out.add(atom);
        return out;
    }*/

    public AtomList noOXTFilter() {
        AtomList out = new AtomList(molecularSystem);
        for (Atom atom : this)
            if (atom.name().trim().compareTo("OXT") != 0) out.add(atom);
        return out;
    }

    public AtomList backbone() {
        return filter(new BackboneFilter());
    }

    public AtomList sideChains() {
        return filter(new SideChainFilter());
    }

    public static class BackboneFilter implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            if (atom.name().equals("N")) return true;
            if (atom.name().equals("H")) return true;
            if (atom.name().equals("CA")) return true;
            if (atom.name().equals("CB")) return true;
            if (atom.name().equals("C")) return true;
            if (atom.name().equals("O")) return true;
            return false;
        }
    }

    public static class SideChainFilter implements Filter {
        private static BackboneFilter bbf = new BackboneFilter();
        public boolean accept(Object obj) {
            return (! bbf.accept(obj));
        }
    }

    public boolean sortable() {
        return false;
    }

    public String toString() {
        return "AtomList with " + size() + " getAtoms";
    }

    public AtomList CAFilter() {
        AtomList out = new AtomList(molecularSystem);
        out.sourceFile = sourceFile;
        for (Atom atom : this)
            if (atom.name().trim().compareTo("CA") == 0) out.add(atom);
        return out;
    }

    public AtomList CAnoImageFilter() {
        AtomList out = new AtomList(molecularSystem);
        out.sourceFile = sourceFile;
        for (Atom atom : this)
            if ((atom.name().trim().compareTo("CA") == 0) && (!atom.core.status().image())) out.add(atom);
        return out;
    }


    public AtomList frozenAtoms() {
        return filter(new IsFrozenFilter());
    }

    public AtomList defrostedAtoms() {
        return filter(new IsDefrostedFilter());
    }

    public static class IsDefrostedFilter implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            return !atom.frozen();
        }
    }

    public static class IsFrozenFilter implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            return atom.frozen();
        }
    }

    public static class NonHydrogen implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            return !(atom.isHydrogen());
        }
    }

    /**
     * Returns true if the parameter has the same length and the same composition in terms of atom names.
     * Current implementation, rather inefficient.
     */
    public boolean sameAtoms(AtomList other) {
        if (size() != other.size()) return false;

        for (int iAtom = 0; iAtom < size(); iAtom++) {
            boolean notFound = true;
            Atom atom = get(iAtom);
            for (int jAtom = 0; (jAtom < size()) & notFound; jAtom++) {
                Atom otherAtom = get(jAtom);
                if (atom.name().equals(otherAtom.name())) notFound = false;
            }
            if (notFound) return false;
        }
        return true;
    }


    /**
     * Asserts that <code>this</code> AtomList is comprised of
     * getAtoms belonging to a single residue.
     *
     * @throws meshi.util.MeshiException if the assertion fails.
     */
    public void assertSingleResidue() {
        if (!checkSingleResidue())
            throw new MeshiException("AtomList.assertSingleResidue -- " +
                    "It appears that this AtomList contains more than one Residue.");
    }

    /**
     * Checks if <code>this</code> AtomList is comprised of
     * getAtoms belonging to a single residue.
     */
    public boolean checkSingleResidue() {
        if (size() == 0) return true;
        Iterator atomIter = iterator();
        Atom atom = (Atom) atomIter.next();
        Residue residueOfFirstAtom = atom.residue();

        while (atomIter.hasNext()) {
            atom = (Atom) atomIter.next();
            if (!Utils.equals(atom.residue(), residueOfFirstAtom))
                return false;
        }

        return true;
    }


    public AtomList[] splitToChains() {
        AtomList leftToSplit = this;
        ArrayList<AtomList> atomLists = new ArrayList<AtomList>();
        String currentChain;

        while (!leftToSplit.isEmpty()) {
            currentChain = leftToSplit.get(0).chain();

            atomLists.add(leftToSplit.filter(new AcceptChainFilter(currentChain)));
            leftToSplit = leftToSplit.filter(new RejectChainFilter(currentChain));
        }

        AtomList[] out = new AtomList[atomLists.size()];
        for (int i = 0; i < out.length; i++)
            out[i] = atomLists.get(i);

        return out;
    }

    public AtomList getChainAtoms(String chainLetter) { //A == 1
        return this.filter(new AcceptChainFilter(chainLetter));
    }

    public AtomList getSourceAtoms() {
        return this.splitToChains()[0];
    }


    public static class AcceptChainFilter implements Filter {
        public final String CHAIN;

        public AcceptChainFilter(String chain) {
            if (chain == null || chain.length() != 1)
                throw new MeshiException
                        ("ChainFilter must receive a one-letter String");

            this.CHAIN = chain;
        }

        public boolean accept(Object obj) {
            return (obj instanceof Atom) &&
                    CHAIN.equals(((Atom) obj).chain());
        }
    }

    public static class RejectChainFilter implements Filter {
        public final String CHAIN;

        public RejectChainFilter(String chain) {
            if (chain == null || chain.length() != 1)
                throw new MeshiException
                        ("ChainFilter must receive a one-letter String");

            this.CHAIN = chain;
        }

        public boolean accept(Object obj) {
            return (obj instanceof Atom) &&
                    !CHAIN.equals(((Atom) obj).chain());
        }
    }

    public AtomList somewhere() {
        AtomList out = new AtomList(molecularSystem);
        for (Iterator atoms = iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (!atom.nowhere()) out.add(atom);
        }
        return out;
    }


    public void moveCMtoOrigin() {
        double cmx = 0.0;
        double cmy = 0.0;
        double cmz = 0.0;
        int counter = 0; // Number of getAtoms with coordinates (all but getAtoms that are nowhere).
        for (Atom atom : this) {
            if (!atom.nowhere()) {
                cmx += atom.x();
                cmy += atom.y();
                cmz += atom.z();
                counter++;
            }
        }
        cmx /= counter;
        cmy /= counter;
        cmz /= counter;
        for (Atom atom : this) {
            if (!atom.nowhere())
                atom.setXYZ(atom.x() - cmx,
                        atom.y() - cmy,
                        atom.z() - cmz);
        }
    }

    public boolean includesNowhereAtoms() {
        for (Iterator atoms = iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (atom.nowhere()) return true;
        }
        return false;
    }

    public void setComment(String s) {
        comment = s;
    }

    public void print() {
        for (Atom a : this) System.out.println(a);
    }

    public void print(PrintWriter writer) {
        print(writer, false);
    }

    public void print(PrintWriter writer, boolean nowhereFlag) {
        for (Atom a : this)
            if (nowhereFlag || !a.nowhere()) writer.println(a);
    }

    public String comment() {
        return comment;
    }

    public AtomList located() {
        AtomList out = new AtomList(molecularSystem);
        for (Atom atom :this)
            if (!atom.nowhere()) out.add(atom);
        return out;
    }
}


