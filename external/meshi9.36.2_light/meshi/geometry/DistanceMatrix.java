/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.util.Terminator;
import meshi.util.file.MeshiWriter;
import meshi.util.mathTools.Jama.Matrix;

import java.util.ArrayList;

/**
 * Where all the headache of {@link Distance distances} is handled.                   <br>
 * <p/>
 * Almost any measurable feature of a molecule is related to distances
 * between pairs of {@link meshi.molecularElements.atoms}.
 * Thus, The calculation of distances (and their inverse and derivatives)
 * are typically a computational bottleneck in structural
 * biology applications. In all applications that we are aware of (Please,
 * enlighten us if you know better) distance calculation is done as part
 * of the procedures that use it (say, as part of the van-der-Waals energy
 * calculation). As a result the distance between two atoms may be calculated
 * more then once. For example the distance between two atoms may be
 * calculated both during angle and torsion angle energies calculations.
 * In Meshi we tried to encapsulate all distance related issues in few
 * classes: {@link Distance Distance}, its subclasses and this one.
 * <p>
 * The motivation behind this class is twofold: first, to keep all the headache of
 * cutoff distances (see below) in a single place. Second, to arrange all the
 * {@link Distance
 * distances                                                                                     }
 * of a molecular system in a single data
 * structure so that no distance is calculated more than once in a single energy
 * evaluation.
 * </p><p>
 * <b>
 * Distance cutoff                                                                               </b><br>
 * Calculating and storing all the distances of a molecular system requires O(n^2) time
 * and storage, where <b>n</b> is the number of atoms. Heuristic algorithms reduces it to O(n),
 * under certain assumptions.
 * The heuristic algorithms (see references below) relays isOn three characteristics
 * of energy functions and energy based simulations:
 * <ol>
 * <li> Atoms have a characteristic minimal distance from other atoms. Thus, the distance
 * matrix is typically dominated by large distances.
 * <li> Energy functions typically decay with distance and become negligible when the
 * atoms are far apart. <br>
 * This implies that even if the energy function is formally defined over all distances,
 * distances above some threshold can be considered "infinite" and having no energy contribution.
 * <li> During most of an energy based simulation (e.g. MD, minimization, MC etc.) the atom
 * movements are rather slow, and the distance matrix changes very little between
 * consecutive simulation steps.<br>
 * This implies that most of the atom pairs that are "infinitely" distant in a given
 * simulation step will remain so in the next step.
 * </ol>
 * The current implementation requires O(n^2) storage and O(n) CPU time. Future
 * implementations are intended to be more efficient.
 * </p><p>
 * <p/>
 * The first step of the algorithm is to separate the setResidue of distances into two groups:<br>
 * <ol>
 * <li> bonded distances - these are distances between atoms that are not likely to be too
 * far apart at any time and thus, algorithms based isOn distance cutoffs are not applicable to them.
 * Typically these are atoms separated by up to three or four covalent bonds.
 * The number of these distances is O(n) so they do not add to the computational complexity.
 * <li> unbonded distances - these are distances between atoms that may or may not be close in space.
 * </ol>
 * <p/>
 * <p/>
 * by two predefined points: rMax & rmax+buffer<br>
 * <ol>
 * <li> 0 < distance <= rMax                                                                              <br>
 * It is assumed that all non zero interactions occure within this distance range.                   <br>
 * The algorithm garenties that if a pair of atoms have a distance within this range
 * it is included in the nonbonded list.</li>
 * <li> rMax < distance <= rMax + buffer                                                               <br>
 * Some of the atom pairs with a distance within this range are included in the nonbonded list.   </li>
 * <li> rMax < distance                                                                                <br>
 * These atom pairs are not included in the nonbonded list</li>
 * <p/>
 * The number of atom-pairs with inter-atomic distances in the first two regions is O(n). These
 * atom-pairs are stored in a list called the non-bonded-list. The only O(n^2) task is to test for
 * each of the O(n^2) atom-pairs currently in the third region whether it moved to the first two. This
 * is where our assumption that changes in the inter-atomic distances are slow enters. Consider a pair of
 * atoms A1 and A2 that were in the coordinates C1(S) and C2(S) at some step S of the simulation
 * such that the distance between them
 * D(C1(S),C2(S)) > rMax+buffer. D(C1(S+T),C2(S+T)) the distance between these atoms at some later
 * step S+T may be smaller or equal to rMax only if at least for one of the atoms
 * D(C(S),C(S+T)) > buffer/2. This condition is tested in O(n) and assuming slow atoms movements
 * should, for most atoms fail. Only when this condition is satisfied for one of the proteins
 * We should check all it's distances from other atoms again with again O(n) complexity.
 * <p/>
 * </p><p>
 */
public class DistanceMatrix extends DistanceLists implements Updateable {
    public final Terminator terminator;
    private static int counter = 0;
    public final int id;
    public static enum DistanceMatrixType {
        STANDARD(5.5, 2.0, DistanceMatrix.DEFAULT_EDGE(5.5, 2.0), 4),
        LONG_DISTANCES(10, 2.0, DistanceMatrix.DEFAULT_EDGE(10, 2.0), 4),
        SHORT_DISTANCES(2.5, 2.0, DistanceMatrix.DEFAULT_EDGE(2.5, 2.0), 4);
        public final double rMax;
        public final double buffer;
        public final double edge;
        public final int bondListDepth;

        private DistanceMatrixType(double rMax, double buffer,
                                   double edge, int bondListDepth){
            this.rMax          = rMax;
            this.buffer        = buffer;
            this.edge          = edge;
            this.bondListDepth = bondListDepth;
        }
    };


    /*------------------------------------------ object variables --------------------------------------------*/
    public static final double DEFAULT_RMAX = 5.5;
    public static final double DEFAULT_BUFFER = 1;
    public static final int DEFAULT_BONDED_LIST_DEPTH = 5;
    private Indicator indicatorToUpdateHB;
    public final DistanceMatrixType type;

    /**
     * The list of all atoms in the molecular system.
     */
    public final MolecularSystem molecularSystem;

    /**
     * Internal data structure.
     */
    protected Grid grid;

    /**
     * Maximal distance for nonzero interactions.
     */
    protected double rMax;
    protected double rMax2;
    protected double edge;

    public static final double DEFAULT_EDGE(double rMax, double buffer) {
        return rMax + 2.0 * buffer / 3.0;
    }

    /**
     * rMax+buffer
     */
    protected static double rMaxPlusBuffer;

    /**
     * (rMax+buffer)^2
     */
    protected static double rMaxPlusBuffer2;

    protected double buffer;

    /**
     * (buffer/3)^2
     */
    protected static double bufferOneThirdSqr;

    /**
     * Atom pairs with inter-atomic distances below rMax (and some of the pairs below rMax+buffer).
     */
    protected DistanceLists nonBondedList,newNonBondedList;


    /**
     * List of DistanceLists
     * Every DistanceLists contains distances needed for one EnergyTerm,
     * selected by its filter
     * For example:
     * - Distances of good hydrogen bonds candidate that were added in the current update opperation
     * - Applicable distances between nonBonded C-N candidate that were added in the current update opperation
     */
    protected ArrayList<DistanceLists> energyTermsDistanceLists;

    public ArrayList<DistanceLists> energyTermsDistanceLists() {
        return energyTermsDistanceLists;
    }

    /**
     * Atom pairs with inter-atomic distances that are always relatively small.
     */

    protected DistanceLists bondedList;

    private int newConstant = 0;
    private boolean nonBondedFlag = true;
    protected int bondedListDepth;
    protected int numberOfUpdates = 0;
    protected boolean debug = false;

    /**
     * Enter DistanceMatrix debug mode.
     */
    public void debugON() {
        debug = true;
    }

    /**
     * Exit DistanceMatrix debug mode.
     */
    public void DebugOFF() {
        debug = false;
    }

    protected final ArrayList<Distance> distancesToUpdate; //Atoms that need to be considered as moving (see MatrixRow)
    /*-------------------------------------------- constructors -----------------------------------------------*/

    public DistanceMatrix(MolecularSystem molecularSystem,DistanceMatrixType type) {
        super(molecularSystem.size());
        this.molecularSystem = molecularSystem;
        this.type = type;
        terminator = molecularSystem.terminator();
        buffer = type.buffer;
        rMax   = type.rMax;
        edge   = type.edge;
        bondedListDepth = type.bondListDepth;
        distancesToUpdate = new ArrayList<Distance>();
        reset();
        id = counter;
        counter++;
    }

    public DistanceMatrix(MolecularSystem molecularSystem) {
        this(molecularSystem, 5.5, 0.2, DistanceMatrix.DEFAULT_EDGE(5.5, 0.2), 4);
    }

    private DistanceMatrix(MolecularSystem molecularSystem, double rMax, double buffer, double edge,
                          int bondedListDepth) {
        this(molecularSystem,DistanceMatrixType.STANDARD);
        // Setting the contants used for calculating the reported distance.
        this.buffer = buffer;
        this.rMax = rMax;
        this.edge = edge;
        this.bondedListDepth = bondedListDepth;
        reset();
    }

    public void reset() {
        rMax2 = rMax * rMax;
        rMaxPlusBuffer = rMax + buffer;
        rMaxPlusBuffer2 = rMaxPlusBuffer * rMaxPlusBuffer;
        bufferOneThirdSqr = buffer * buffer / 9;
        indicatorToUpdateHB = new Indicator();

        // Initialize the distance matrix.
        clear();
        nonBondedList = new DistanceLists(molecularSystem.size());
        newNonBondedList = new DistanceLists(molecularSystem.size());
        bondedList = new DistanceLists(molecularSystem.size());
        for (AtomCore atom : molecularSystem) {
//             if (atom.status().active()) {
                DistanceList nonBondedRow = new DistanceList(DistanceLists.ROW_INITIAL_CAPACITY, atom);
                DistanceList newNonBondedRow = new DistanceList(1, atom);
                DistanceList bondedRow = new DistanceList(DistanceLists.ROW_INITIAL_CAPACITY, atom);
                add(new MatrixRow(atom, this, nonBondedRow, newNonBondedRow));
                nonBondedList.add(nonBondedRow);
                newNonBondedList.add(newNonBondedRow);
                bondedList.add(bondedRow);//The matrix row does not need to now about the bonded row. It does not update it.
//            }
        }
        /**
         * Every element of energyTermsDistanceLists
         * can be specified in Creator of EnergyTerm needed it (for example HydrogenBondsCreator)
         */
        energyTermsDistanceLists = new ArrayList<DistanceLists>();

        // Generate bonded list.
        getBondedList(molecularSystem, bondedListDepth);
        //bonded list.

        // We adjust the grid cell edge size so that the memory requirements of the
        // grid will not result in memory failure.
        try {
            while ((grid = new Grid(molecularSystem, edge, DEFAULT_EDGE(rMax, buffer))).status() != GridStatus.OK) {
                if (edge >= 100000) throw new RuntimeException("Very weird");
                edge *= 2;
            }
            update();
        }
        catch (UpdateableException e) {
            throw new RuntimeException(" Cannot create grid. Apparently the atoms have spreaded over" +
                    " a very large volume.");
        }
        if (MatrixRow.debug)
            MatrixRow.debugCounter = 0;
    }


    /*------------------------------------------ update --------------------------------------------*/

    /**
     * Updates the distance matrix.
     */
    public void update(int numberOfUpdates) throws UpdateableException {
        if (terminator.dead())
            throw new RuntimeException("Distance matrix \n" + this + "\n is dead. nDebug = "+Terminator.nDebug+"\n" + terminator.message());
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);

    }


    public void update() throws UpdateableException {
        if (terminator.dead())
            throw new RuntimeException("Distance matrix \n" + this + "\n is dead. nDebug = "+Terminator.nDebug+"\n" + terminator.message());

        for (DistanceList matrixRow : this) {
            ((MatrixRow) matrixRow).update();
            if (terminator.dead())
                throw new RuntimeException("In matrix row:"+matrixRow.atomOne+" Distance matrix \n" + this + "\n is dead. nDebug = "+Terminator.nDebug+"\n" + terminator.message());
        }

        GridStatus gridStatus = grid.build();
        if (gridStatus != GridStatus.OK) {
            System.out.println("Distance Matrix of type "+type+"\n"+"rMax = "+rMax+"\nedge = "+edge);
            throw new UpdateableException(gridStatus);
        }
        if (MatrixRow.debug)
            MatrixRow.debugCounter++;
        distancesToUpdate.clear(); //the rows will update it for the next round
        for (DistanceList row : this) {
            if (!row.atomOneNowhere) {
                GridCell gridCell = grid.getCell(row.atomOne);
                ((MatrixRow) row).addCell(gridCell);
            }
        }
    }



    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }


    /*------------------------------------------ other methods --------------------------------------------*/



    /**
     * Returns the non-bonded-list.
     */
    public DistanceLists nonBondedList() {
        return nonBondedList;
    }
    public DistanceLists newNonBondedList() {
        return newNonBondedList;
    }
    public DistanceLists bondedList() {
        return bondedList;
    }

//    public DistanceLists newNonBondedList() {return newNonBondedList;}

    /**
     * Returns the Distance object of the parameters.
     */
    public Distance distance(Atom atom1, Atom atom2) {
        try {
            return distance(atom1.number(), atom2.number());
        } catch (RuntimeException ex) {
            System.out.println(" A problem while finding the distance between:\n" +
                    "atom1 = " + atom1 + "\n"+atom1.core+"\nand\n" + "atom2 = " + atom2 + "\n" +atom2.core+"\n"+
                    "atom1.nowhere() = " + atom1.nowhere() + " ; " + "atom2.nowhere() = " + atom2.nowhere() + "\n");
            ex.printStackTrace();
            throw ex;
        }
    }

    public Distance distance(AtomPair atomPair) {
        return distance(atomPair.largeNumber(), atomPair.smallNumber());
    }

    /**
     * Returns the Distance object of the parameters.
     */
    public Distance distance(int atom1Number, int atom2Number) {
        if (terminator.dead()) throw new RuntimeException("This DistanceMatrix is dead. The terminator message is "+terminator.message());
//        try {
            if ((atom1Number < 0) ||
                    (atom2Number < 0) ||
                    (atom1Number >= size()) ||
                    (atom2Number >= size()) ||
                    (get(atom1Number) == null) ||
                    (get(atom2Number) == null)) //return null;
                    throw new RuntimeException("\n atom1Number " + atom1Number +"\n" +
                    "atom2Number " + atom2Number + "\n" +
                    "matrix.length " + size()+"\n");
/*
        } catch (RuntimeException ex) {
            System.out.println("atom1Number " + atom1Number + "\n" +
                    "atom2Number " + atom2Number + "\n" +
                    "matrix.length " + matrix.length);
            throw ex;
        }
*/
        Distance out;
            out = ((MatrixRow) search(atom1Number)).search(atom2Number);
       return out;
    }


    public double radius() {
        return Utils.radiusOfGyration(molecularSystem.toAtomList());
    }

    public String toString() {
        return ("DistanceMatrix "+id+" :\n" +
                "\t type \t"+type+
                "\t number of atoms \t" + molecularSystem.size() +
                "\t rMax \t" + rMax +
                "\t buffer\t" + buffer);
    }


    /**
     * Returns the bonded list
     */
    public int nonBondedListSize() {
        int size = 0;
        for (DistanceList distanceList : nonBondedList) {
            size += distanceList.size();
        }
        return size;
    }

    public void getBondedList(MolecularSystem molecularSystem, int depth) {
        AtomList bonded;
        Distance distance;
        double dx,dy,dz,d2;
        AtomCore atomOne;
        for (DistanceList  matrixRow : this) {
            atomOne = matrixRow.atomOne;
                bonded = getBonded(atomOne.atom, depth);
                for (Atom bondedAtom : bonded)   {
                    if ( bondedAtom.active()) {
                        dx = atomOne.x() - bondedAtom.x();
                        dy = atomOne.y() - bondedAtom.y();
                        dz = atomOne.z() - bondedAtom.z();
                        d2 = dx * dx + dy * dy + dz * dz;
                        if (atomOne.atom.frozen() & bondedAtom.frozen())
                            distance = new BondedFrozenDistance(atomOne.atom, bondedAtom,
                                    dx, dy, dz, Math.sqrt(d2));
                        else distance = new BondedDistance(atomOne.atom, bondedAtom, dx, dy, dz, Math.sqrt(d2));

                        bondedList.search(atomOne.number).add(distance);
                        matrixRow.add(distance);
                    }
                }
            }
    }



    public static AtomList getBonded(Atom atom, int depth) {
        AtomList out = new AtomList(atom.molecularSystem);
        getBonded(atom, depth, out, atom.number());
        return out;
    }

    public static void getBonded(Atom atom, int depth, AtomList out, int rootNumber) {
        if (depth == 0) return;
        for (Atom bondedAtom : atom.bonded()) {
            if ((rootNumber > bondedAtom.number()) &
                    (!out.contains(bondedAtom)))
                out.add(bondedAtom);
            getBonded(bondedAtom, depth - 1, out, rootNumber);
        }

    }



    public double rMax() {
        return rMax;
    }


    public double buffer() {
        return buffer;
    }


    public class Indicator {
        public Indicator() {
        }
    }


    public void writeNonBondedList(String s) {
        MeshiWriter writer=null;
        try {
            writer = new MeshiWriter("nonBondedList_"+s+".txt");
            if (nonBondedList == null)
                throw new RuntimeException(nonBondedList+" == null");
            for (DistanceList row : nonBondedList) {
                writer.println("--------------------------------------------------------");
                writer.println(" row "+row.atomOne);
                for (Distance distance :row) {
                    writer.println(distance);
                }
            }
            writer.close();
        } catch (Exception ex) {
            System.out.println("writeNonBondedList failed. "+ex+"\n"+ex.getMessage());
            ex.printStackTrace();
            if (writer != null) writer.close();
            throw new RuntimeException(ex);}
    }



}
