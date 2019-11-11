/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeSumma5PolarTypesEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandran;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandran;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.util.CommandList;
import meshi.util.Terminator;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.filters.Filter;
import meshi.util.formats.Fdouble;
import meshi.util.formats.Fint;
import meshi.util.formats.Format;
import meshi.util.info.ChainsInfo;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import programs.Optimize;

import java.util.*;

/**
 * A generic class for all meshi energy functions. Typically only a single such object exists, which
 * include several energy term objects (sub classes of AbstractEnergy).
 */
public class TotalEnergy {
    /**
     * Number of times the energy function was evaluated. This is a simple machine-independent
     * measure of convergence efficiency.
     */
    protected int numberOfEvaluations;
    protected int numberOfUpdates;
    protected int numberOfReports;
    protected AtomList atomList;
    protected double[][] coordinates;
    protected ArrayList<AbstractEnergy> energyTerms;
    protected ArrayList<Double> energyValues;
    protected double totalEnergy;
    protected Fdouble dformat = Fdouble.SHORT;
    protected Format sformat = Format.SHORT;
    protected Fint iformat = Fint.SHORT;
    protected final DistanceMatrix distanceMatrix;

    public final Terminator terminator;
    public final long startTime;
    public final String comment;
    private static final double INFINITY = 1 / 0.0;
    private DoubleInfoElement infoElement;
    private final MolecularSystem molecularSystem;
    private boolean on;


    public TotalEnergy(AtomList atomList,
                       AbstractEnergy[] energyTerms, String comment) throws UpdateableException{
        this(atomList, comment);
        for (int i = 0; i < energyTerms.length; i++)
            this.energyTerms.add(energyTerms[i]);
    }

    public TotalEnergy(Protein protein,
                       EnergyCreator[] energyCreators,
                       CommandList commands, String comment) {
        this(protein, energyCreators, DistanceMatrix.DistanceMatrixType.STANDARD, commands, comment);
    }
    public TotalEnergy(Protein protein,
                       EnergyCreator[] energyCreators, DistanceMatrix.DistanceMatrixType distanceMatrixType,
                       CommandList commands, String comment){
        this(protein.atoms(), distanceMatrixType, comment);
        energyTerms = new ArrayList<AbstractEnergy>();
        for (int i = 0; i < energyCreators.length; i++) {
            if (!energyCreators[i].weightWasSet()) {
                energyCreators[i].getWeight(commands);
            }
            if (energyCreators[i].weight() != 0) {
                AbstractEnergy term = energyCreators[i].createEnergyTerm(protein,
                        distanceMatrix,
                        commands);
                term.setNumberOfUpdates(0);   // chen 3/3/12
                energyTerms.add(term);
            }
        }
    }

    public TotalEnergy(AtomList atomList,
                       String comment) throws UpdateableException{
        this(atomList,DistanceMatrix.DistanceMatrixType.STANDARD,comment);
     }
    public TotalEnergy(AtomList atomList, DistanceMatrix.DistanceMatrixType distanceMatrixType,
                       String comment) {
        molecularSystem = atomList.molecularSystem();
        //This object is turned on and any other TotalEnergy
        //Object associated with the molecular system is turned off
        molecularSystem.register(this);
        on = true;
        this.atomList = atomList;
        coordinates = getCoordinates(atomList);
        molecularSystem.createDistanceMatrix("Some atoms might have gone to nowhere");
        this.distanceMatrix = molecularSystem.getDistanceMatrix();
        if (molecularSystem.terminator().dead())
            throw new RuntimeException("This is weird.\n"+molecularSystem.terminator().message());
        terminator = molecularSystem.terminator();

        numberOfEvaluations = 0;
        numberOfUpdates = 0;
        numberOfReports = 0;
        energyTerms = new ArrayList<AbstractEnergy>();
        resetAtomEnergies();
        startTime = (new Date()).getTime();
        this.comment = comment;
        Utils.println("Created "+this+" and made it current");
    }


    public TotalEnergy(Protein protein,
                       //DistanceMatrix distanceMatrix,
                       EnergyCreator[] energyCreators,
                       CommandList commands, AtomList atoms,String comment) throws UpdateableException{
        this(atoms, comment);
        if (distanceMatrix.terminator.dead()) throw new RuntimeException("This is weird. The message is: "+distanceMatrix.terminator.message());
        energyTerms = new ArrayList<AbstractEnergy>();
        for (int i = 0; i < energyCreators.length; i++) {
            energyCreators[i].getWeight(commands);
            if (energyCreators[i].weight() != 0) {        //chen 3/3/12
                AbstractEnergy term = energyCreators[i].createEnergyTerm(protein,
                        distanceMatrix,
                        commands);
                term.setNumberOfUpdates(0);
                energyTerms.add(term);
            }
        }
    }





    public String toString() {return "TotalEnergy - "+comment;}

    public DistanceMatrix distanceMatrix() {
        return distanceMatrix;
    }


    public void add(TotalEnergy energy) {
        for (Atom atom : energy.atomList)
            atomList.add(atom);
        for (AbstractEnergy ae : energy.energyTerms)
            energyTerms.add(ae);
    }

    public void test(Exception exception) throws Exception {
        // if(distanceMatrix!=null){
        // distanceMatrix.testNonBondedList();
        // }

        System.out.println("**************************************************************************\n" +
                "            Total Energy Test after exception " + exception + "\n");
        StackTraceElement[] trace = exception.getStackTrace();
        for (StackTraceElement element : trace) System.out.println(element);
        System.out.println("**************************************************************************\n" + "Now testing differentiation\n");
        test();
        throw exception;
    }

    public void test () {
        System.out.println("===> TotalEnergy test");
        System.out.println("===>NonBondedListTest");
        Utils.testNonBondedList(distanceMatrix.nonBondedList());
        System.out.println("===>end of NonBondedListTest");
        Atom criminal = findCriminalAtom();
        try {
            if (criminal != null) {
                System.out.println("*** Criminal found:");
                System.out.println(criminal);
                AbstractEnergy criminalEnergy = findCriminalEnergyTerm(criminal);

                if (criminalEnergy != null) {
                    System.out.println("****************** Criminal " + criminalEnergy + " *******************");
                    criminalEnergy.test(this, criminal);
                }
                //      else{

                for (Iterator ti = energyTerms.iterator(); ti.hasNext(); ) {
                    AbstractEnergy energyTerm = (AbstractEnergy) ti.next();
                    System.out.println("****************** Testing " + energyTerm + " *******************");
                    energyTerm.setDebug(criminal);
                    energyTerm.test(this, criminal);
                    energyTerm.resetDebug();
                }
                //}
            } else {
                System.out.println("*** Criminal atom not found");
                System.out.println("*** Test terminated");
            }
        } catch (EvaluationException ex) {throw new RuntimeException(ex);}
    }

    final static double DX = 1e-9;
    final static String[] XYZ = new String[]{"x", "y", "z"};
    final static double verySmall = Math.exp(-15);

    /**
     * Searching for a `criminal' atom
     *
     * @return a criminal <code>Atom</code>
     */
    protected Atom findCriminalAtom()  {
        System.out.println("*** Criminal atom searching");

        double[][] coordinates = new double[3][];
        double x, e1 = 999999.9999999, e2 = 999999.9999999;
        double analiticalForce = 99999.99, numericalForce = -999999.99;
        double diff = 99999.9999, maxDiff = 0;
        Atom criminal = null;
        criminal_found:
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (!atom.frozen()) {
                coordinates[0] = atom.X();
                coordinates[1] = atom.Y();
                coordinates[2] = atom.Z();
                for (int i = 0; i < 3; i++) {
                    // Whatever should be updated ( such as distance matrix torsion list etc. )
                    updateDebug();
                    x = coordinates[i][0];
                    coordinates[i][1] = 0;
                    e1 = evaluate();
                    analiticalForce = coordinates[i][1];
                    coordinates[i][0] += DX;

                    updateDebug();
                    e2 = evaluate();
                    numericalForce = -1 * (e2 - e1) / DX;
//                    if ((numericalForce == 0) && (Math.abs(analiticalForce) > 1))
//                        throw new RuntimeException("Numerical force On X axis of "+atom+" is zero\n"+
//                                "e1 = "+e1+" e2 = "+e2 +" DX = "+DX+"\n"+
//                                "Analytic force is " + analiticalForce);
                    coordinates[i][0] -= DX;
                    updateDebug();

                    diff = Math.abs(analiticalForce - numericalForce);
                    if (maxDiff < diff) {// diff is maximal
                        maxDiff = diff;
                        System.out.println();
                        System.out.println("Atom[" + atom.number() + " " + atom.residueNumber() + "]." + XYZ[i] + " = " + x);
                        System.out.println("Analytical force = " + analiticalForce);
                        System.out.println("Numerical force  = " + numericalForce);
                        System.out.println("maxDiff = " + maxDiff);
                        System.out.println("tolerance = " + 2 * maxDiff / (Math.abs(analiticalForce) +
                                Math.abs(numericalForce) + verySmall));
                        criminal = atom;
                        // break criminal_found;
                    }

                    if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                        System.out.println("e1 = " + e1);
                    if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                        System.out.println("e2 = " + e2);
                    if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                        System.out.println("analiticalForce = " + analiticalForce);
                }
                if (atom.number() != 0 && atom.number() % 100 == 0)
                    System.out.print(atom.number());
                else
                    System.out.print(".");
            }
        }
        System.out.println();
        return criminal;
    }

    /**
     * Searching for a `criminal' energy term , with specific atom
     *
     * @param atom a `criminal' <code>Atom</code>
     * @return a criminal energy term
     */
    protected AbstractEnergy findCriminalEnergyTerm(Atom atom){
        //CooperativeSumma5PolarTypesEnergy[] cooperativeZSumma = (CooperativeSumma5PolarTypesEnergy[])    getEnergyTerms(new CooperativeSumma5PolarTypesEnergy());
        AtomicPairwisePMFSumma atomicPairwisePMFSumma = (AtomicPairwisePMFSumma) getEnergyTerm(new AtomicPairwisePMFSumma());

        CooperativeZPropensityEnergy cooperativeZPropensity = (CooperativeZPropensityEnergy) getEnergyTerm(new CooperativeZPropensityEnergy());
        CooperativeZStdPropensityEnergy cooperativeZStdPropensity = (CooperativeZStdPropensityEnergy) getEnergyTerm(new CooperativeZStdPropensityEnergy());
        CompositePropensityEnergy propensityEnergy = (CompositePropensityEnergy) getEnergyTerm(new CompositePropensityEnergy());

        CooperativeZRamachandran cooperativeZRamachandran = (CooperativeZRamachandran) getEnergyTerm(new CooperativeZRamachandran());
        CooperativeZStdRamachandran cooperativeZStdRamachandran = (CooperativeZStdRamachandran) getEnergyTerm(new CooperativeZStdRamachandran());
        RamachandranSidechainEnergy ramachandranEnergy = (RamachandranSidechainEnergy) getEnergyTerm(new RamachandranSidechainEnergy());

        System.out.println("*** Criminal energy term searching");

        double[][] coordinates = new double[3][];
        double x, e1 = 999999.9999999, e2 = 999999.9999999;
        double analiticalForce = 99999.99, numericalForce = -999999.99;
        double diff = 99999.9999, maxDiff = 0;
        AbstractEnergy criminal = null;
        for (Iterator ti = energyTerms.iterator(); ti.hasNext();) {
            AbstractEnergy energyTerm = (AbstractEnergy) ti.next();
            if (!energyTerm.on) {
                System.out.println(energyTerm + " is OFF");
                continue;
            }
            System.out.println("Now testing " + energyTerm);
            coordinates[0] = atom.X();
            coordinates[1] = atom.Y();
            coordinates[2] = atom.Z();
            for (int i = 0; i < 3; i++) {
                // Whatever should be updated ( such as distance matrix torsion list etc. )
                updateDebug();
                x = coordinates[i][0];
                coordinates[i][1] = 0;
                try {

                //if ((cooperativeZSumma != null) && energyTerm.toString().equals(cooperativeZSumma.toString()))
                if (energyTerm instanceof CooperativeSumma5PolarTypesEnergy) {
                    atomicPairwisePMFSumma.evaluate();
                    coordinates[i][1] = 0;
                } else if ((energyTerm instanceof CooperativeZPropensityEnergy) || (energyTerm instanceof CooperativeZStdPropensityEnergy))
                    propensityEnergy.evaluate();
                else if ((energyTerm instanceof CooperativeZRamachandran) || (energyTerm instanceof CooperativeZStdRamachandran))
                        ramachandranEnergy.evaluate();


                    e1 = energyTerm.evaluate().energy();
                } catch (EvaluationException e) {
                    Utils.throwException(this, e, "Failed to evaluate " + energyTerm);
                }
                analiticalForce = coordinates[i][1];
                coordinates[i][0] += DX;

                updateDebug();

                try {
                    if (energyTerm instanceof CooperativeSumma5PolarTypesEnergy)
                    //                if ((cooperativeZSumma != null) && energyTerm.toString().equals(cooperativeZSumma.toString()))
                    atomicPairwisePMFSumma.evaluate();
                else if ((energyTerm instanceof CooperativeZPropensityEnergy) || (energyTerm instanceof CooperativeZStdPropensityEnergy))
                    propensityEnergy.evaluate();
                else if ((energyTerm instanceof CooperativeZRamachandran) || (energyTerm instanceof CooperativeZStdRamachandran))
                        ramachandranEnergy.evaluate();

                    e2 = energyTerm.evaluate().energy();
                } catch (EvaluationException e) {
                    Utils.throwException(this, e, "Failed to evaluate " + energyTerm + " in e2.");
                }
                numericalForce = -1 * (e2 - e1) / DX;
                coordinates[i][0] -= DX;

                updateDebug();

                diff = Math.abs(analiticalForce - numericalForce);
                if (maxDiff < diff) {// diff is maximal
                    maxDiff = diff;
                    System.out.println();
                    System.out.println("****** test in " + energyTerm);
                    //System.out.println("Atom["+atom.number()+" "+atom.residueNumber()+"]."+XYZ[i]+" = "+x);
                    System.out.println("Atom[" + atom + "]." + XYZ[i] + " = " + x);
                    System.out.println("Analytical force = " + analiticalForce);
                    System.out.println("Numerical force  = " + numericalForce);
                    System.out.println("maxDiff = " + maxDiff);
                    System.out.println("tolerance = " + 2 * maxDiff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + verySmall));
                    criminal = energyTerm;
                }

                if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                    System.out.println("e1 = " + e1);
                if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                    System.out.println("e2 = " + e2);
                if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                    System.out.println("analiticalForce = " + analiticalForce);
            }
        }
        return criminal;
    }

    /**
     * Finds the maximal component (in magnitude) of the gradient vecor in coordinates ( coordinates[][1] ).
     */
    public double getGradMagnitude() {
        double sumGrad = 0;
        int n = 0;
        double max = -1000;
        for (int i = 0; i < coordinates.length; i++) {
            if (coordinates[i][1] != INFINITY) {
                sumGrad += coordinates[i][1] * coordinates[i][1];
                n++;
            }
        }
        if (n != 0) {
            return Math.sqrt(sumGrad / n);
        }
        return 0;
    }

    public double[][] coordinates() {
        return coordinates;
    }

    public int numberOfEvaluations() {
        return numberOfEvaluations;
    }

    /**
     * Sets all forces to zero.
     */
    public static void resetForces(double[][] coordinates) {
        int length = coordinates.length;
        for (int i = 0; i < length; i++) {
            coordinates[i][1] = 0.0;
        }
    }

    public void setCoordinates(AtomList atomList) {
        coordinates = getCoordinates(atomList);
    }

    public void getCoordinates() {
        coordinates = getCoordinates(atomList);
    }

    /**
     * Reduce the Atom coordinates of an atom list to an array.
     * <X coordiante of atom 1><force reset X coordiante of atom 1>
     * <Y coordiante of atom 1><force reset Y coordiante of atom 1>
     * <Z coordiante of atom 1><force reset Z coordiante of atom 1>
     * :
     * :
     * :
     * :
     * <X coordiante of atom N><force reset X coordiante of atom N>
     * <Y coordiante of atom N><force reset Y coordiante of atom N>
     * <Z coordiante of atom N><force reset Z coordiante of atom N>
     */
    public static double[][] getCoordinates(AtomList atomList) {
        AtomList tempList = atomList.filter(new FreeAtom());
        int numberOfCoordinates = tempList.size() * 3;
        double[][] coordinates = new double[numberOfCoordinates][];
        int i = 0;
        for (Atom atom : tempList) {
            coordinates[i++] = atom.X();
            coordinates[i++] = atom.Y();
            coordinates[i++] = atom.Z();
        }
        return coordinates;
    }

    private static class FreeAtom implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            //return ((atom.core.status() == AtomStatus.NORMAL) || (atom.core.status() == AtomStatus.NOWHERE)|| (atom.core.status() == AtomStatus.HIDDEN));
            return ((Atom) obj).normal();
        }
    }

    public void evaluateAtoms() throws UpdateableException, EvaluationException {
        if (!on) throw new RuntimeException("This Total energy "+this+" is off");
        Utils.println("evaluateAtoms");
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            ((Atom) atoms.next()).resetEnergy();
        }

        evaluate();

        for (AbstractEnergy energyTerm : energyTerms) {
            energyTerm.evaluateAtoms();
            double sum = 0, prevSum = 0;
            double atomEnergy;
            for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
                atomEnergy = ((Atom) atoms.next()).energy();
                sum += atomEnergy;
            }
        }
    }

    public void evaluateAtoms(AbstractEnergy TestedEnergy) throws UpdateableException, EvaluationException {
        Utils.println("evaluateAtoms");
        for (Object anAtomList1 : atomList) {
            ((Atom) anAtomList1).resetEnergy();
        }

        evaluate();

        TestedEnergy.evaluateAtoms();
        for (Object anAtomList : atomList) {
            ((Atom) anAtomList).energy();

        }

    }


    public double evaluateAll(ChainsInfo chainsInfo)  {
        return evaluate(true, chainsInfo); // Evaluate also non-differential energy terms.
    }
    public double evaluateAll()  {
        return evaluate(true, null); // Evaluate also non-differential energy terms.
    }



    public double evaluate()   {
        return evaluate(false, null); // Do not evaluate non-differential energy terms.
    }
    private  double evaluate(boolean nonDifferentialFlag, ChainsInfo chainsInfo)  {

        if (!on) throw new RuntimeException("This Total energy "+this.toString()+" is off. Terminator message is: "+terminator.message());
        boolean differential;

        if (molecularSystem.currentTotalEnergy() != this)
            throw new RuntimeException("Only one instance of " +
                    "(a class that extends) TotalEnergy " +
                    "is allowed per molecular system. \n"+
                    "This energy: "+this+"\nCurrent energy:"+molecularSystem.currentTotalEnergy());

        if (terminator.dead())
            throw new RuntimeException("The terminator killed this job with the message\n" + terminator.message());

        energyValues = new ArrayList<Double>();
        double e;
        totalEnergy = 0;

        if (distanceMatrix.terminator.dead()) {
            throw new RuntimeException ("This is very weird. "+distanceMatrix+" "+distanceMatrix.terminator);
        }

        update();

        for (AbstractEnergy energyTerm : energyTerms) {

            if (energyTerm.type == EnergyType.UNKNOWN)
                throw new RuntimeException("Weird energyTerm "+energyTerm+" with type "+energyTerm.type);
            differential = (energyTerm.type == EnergyType.DIFFERENTIAL);
            if (nonDifferentialFlag || differential) {
                try {
                    e = energyTerm.genericEvaluate().energy();
                } catch (EvaluationException e1) {
                    System.out.println("*******************************************\nError in Total energy/update while updating "+energyTerm);
                    e1.printStackTrace();
                    throw new RuntimeException("Quiting update");
                }
                if ((e != 0) & (!energyTerm.isOn()))
                        throw new RuntimeException("This is weird: " + energyTerm +
                            " is off but does not return zero (" + e + ")");
                    if ((!(e < 0)) & (!(e == 0)) & (!(e > 0)))
                        throw new RuntimeException("weird energy #1 " + energyTerm + " " + e);
                  /*  if ((!nonDifferentialFlag) && ((e > 1.0E18) ||   (e < -1.0E18)))
                        throw new RuntimeException("weird energy 2 " + energyTerm + " " + e);     */
                    energyValues.add(new Double(e));
                if (differential) {
                    if (energyTerm.type == EnergyType.NON_DIFFERENTIAL)
                        throw new RuntimeException("This is exceptionally weird "+energyTerm+" "+energyTerm.type);
                     totalEnergy += e;
                }
            }
            if ((chainsInfo != null) & energyTerm.evaluatesResidues())
                ((EvaluatesResidues) energyTerm).evaluateResidues(chainsInfo);
        }
        numberOfEvaluations++;
        return totalEnergy;
    }

    /**
     * Returns the last enrgy value that was calculated
     */
    public double getLastEnergy() {
        return (totalEnergy);
    }

    /**
     * Updates all factors related to the energy function (for example the distance matrix).
     */
    public void update() {
        resetForces();
        numberOfUpdates++;
        for (AbstractEnergy energyTerm : energyTerms) {
            try {
                if (energyTerm == null) {
                    System.out.println("This is weird.");
                    throw new RuntimeException("energy term is null.");
                }
                energyTerm.update(numberOfUpdates);
            }
            catch (UpdateableException ex) {
                System.out.println("************************************\nFailed to update " + energyTerm + "in "+this+"\n" + ex);
                ex.printStackTrace();
                throw new RuntimeException("Quiting update");
            }

        }
    }

    public void updateDebug() {
        if (distanceMatrix != null) {
            distanceMatrix.debugON();
        }
            update();
    }


    public static double getAverageForce(double[][] coordinates) {
        double averageForce = 0;
        for (int i = 0; i < coordinates.length; i++)
            averageForce += coordinates[i][1] * coordinates[i][1];
        averageForce = Math.sqrt(averageForce / coordinates.length);
        return averageForce;
    }

    public void resetForces() {
        resetForces(coordinates);
    }

    public String reportHeader() {
        String report = "";
        report += sformat.f("#") + sformat.f(" E total ") + " ";
        for (AbstractEnergy energyTerm : energyTerms) {
            if (energyTerm.differential()){
                report += sformat.f(energyTerm.comment()) + " ";
            }
        }
        report += (sformat.f("meanF") + " " + sformat.f("radius") +
                " " + sformat.f("#NB") + sformat.f("#Eeval") + sformat.f("Time"));
        return report;
    }

    public String report(int step) {

        boolean allZero = true;
        for (Double value : energyValues) {
            if (value.doubleValue() != 0) allZero = false;
        }
        if (allZero)
            throw new RuntimeException("All energy Values are zero.");
        String report = "";
        if (numberOfReports % 10 == 0) {
            report += reportHeader() + "\n";
        }
        numberOfReports++;
        double averageForce = getAverageForce(coordinates);
        double maxForce = getGradMagnitude();
        report += iformat.f(step) + dformat.f(totalEnergy) + " ";
        for (Double value : energyValues) {
            report += dformat.f(value.doubleValue()) + " ";
        }
        report += Fdouble.SHORTER.f(maxForce) + dformat.f(atomList.radius()) + " ";
        if (distanceMatrix != null) {
            report += iformat.f(distanceMatrix.nonBondedListSize()) + " ";
        }
        report += iformat.f(numberOfEvaluations);
        report += "     " + ((new Date()).getTime() - startTime);
        return report;
    }

    public double energy() {
        return totalEnergy;
    }

    public boolean frozenAtomsExist() {
        return atomList.frozenAtomsExist();
    }

    private static class IsEnergy implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof AbstractEnergy);
        }
    }


    public AtomList atomList() {
        return atomList;
    }

    public ArrayList<Double> energyValues() {
        return energyValues;
    }

    /*
     * Returns a specific energy term, according to its type
     */

    public AbstractEnergy getEnergyTerm(AbstractEnergy ae) {
        for (AbstractEnergy energyTerm : energyTerms) {
            if (energyTerm.getClass().equals(ae.getClass())) {
                return energyTerm;
            }
        }
        return null;
    }

    public DistanceMatrix getDistanceMatrix() {
        return distanceMatrix;
    }

    public final int numberOfUpdates() {
        return numberOfUpdates;
    }

    public void addTerm(AbstractEnergy term) {
        energyTerms.add(term);
    }

    public void resetAtomEnergies() {
        for (Atom atom : atomList) {
            atom.resetEnergy();
        }
    }

    /*
    * Returns a all instances of a specific energy term, according to its type
    */

    public AbstractEnergy[] getEnergyTerms(AbstractEnergy ae) {
        int counter = 0;
        for (AbstractEnergy energyTerm : energyTerms) {
            if (energyTerm.getClass().equals(ae.getClass())) {
                counter++;
            }
        }
        AbstractEnergy[] result = new AbstractEnergy[counter];
        counter = 0;
        for (AbstractEnergy energyTerm : energyTerms) {
            if (energyTerm.getClass().equals(ae.getClass())) {
                result[counter] = energyTerm;
                counter++;
            }
        }
        return result;
    }

    public double totalEnergy() {
        return totalEnergy;
    }

    public double avgEnergy() {
        return totalEnergy / atomList.size();
    }

    public ArrayList<AbstractEnergy> energyTerms() {
        return energyTerms;
    }

    public void summary() throws UpdateableException, EvaluationException {
        if (Utils.verbose()) {
            evaluate();
            Iterator terms = energyTerms.iterator();
            System.out.printf("\nENERGY SUMMARY\n");
            System.out.printf("%30s = %-15.8f\n", "TotalEnergy", totalEnergy);
            System.out.printf("%30s = %-15.8f\n", "avgEnergy", avgEnergy());

            for (Double value : energyValues) {
                String name = ((AbstractEnergy) terms.next()).comment();
                System.out.printf("%30s = %-15.8f\n", name, value);
            }
            SortedAtoms sortedAtoms = new SortedAtoms();
            System.out.printf("%30s = %-15.8f\n", "\"cold1\" getAtoms energy (avg+std) = ", sortedAtoms.sum(1));
            System.out.printf("%30s = %-15.8f\n", "\"cold2\" getAtoms energy (avg+2*std) = ", sortedAtoms.sum(2));
            System.out.printf("%30s = %-15.8f\n", "average \"cold1\" getAtoms energy (avg+std)",
                    sortedAtoms.sum(1) / sortedAtoms.lowestEnergyAtoms(1).size());
            System.out.printf("%30s = %-15.8f\n", "average \"cold2\" getAtoms energy (avg+2*std)",
                    sortedAtoms.sum(2) / sortedAtoms.lowestEnergyAtoms(2).size());
        }
    }


    /*  public double filteredEnergy(double stdThreshold) {
      SortedAtoms sortedAtoms = new SortedAtoms();
      return sortedAtoms.sum(stdThreshold);
  }

  public double avgFilteredEnergy(double stdThreshold) {
      SortedAtoms sortedAtoms = new SortedAtoms();
      return sortedAtoms.sum(stdThreshold)/sortedAtoms.lowestEnergyAtoms(stdThreshold).size();
  }

  public AtomList highestEnergyAtoms(double stdThreshold) {
      SortedAtoms sortedAtoms = new SortedAtoms();
      return sortedAtoms.highestEnergyAtoms(stdThreshold);
  }  */


    public static class EnergyComparator implements Comparator {
        public int compare(Object obj1, Object obj2) {
            Atom atom1 = (Atom) obj1;
            Atom atom2 = (Atom) obj2;
            if (atom1.energy() > atom2.energy()) return 1;
            if (atom1.energy() < atom2.energy()) return -1;
            return 0;
        }
    }


    public void off() {
        on = false;
        for (Iterator terms = energyTerms.iterator();
             terms.hasNext();) {
            AbstractEnergy term = (AbstractEnergy) terms.next();
            term.off();
        }
    }
    public void offAllBut(AbstractEnergy term){
        for(AbstractEnergy ae : energyTerms)
           ae.off();
        term.on();
    }

    public void on() {
        molecularSystem.register(this);
        on = true;
        distanceMatrix.setNumberOfUpdates(numberOfUpdates());
        for (AbstractEnergy energyTerm : energyTerms) {
            energyTerm.setNumberOfUpdates(numberOfUpdates);
            energyTerm.on();
        }
    }


    private class SortedAtoms {
        private Atom[] atoms = new Atom[atomList.size()];
        double avg, std;

        public SortedAtoms() throws UpdateableException, EvaluationException {
            for (int i = 0; i < atoms.length; i++)
                atoms[i] = (Atom) atomList().get(i);
            evaluateAtoms();
            Arrays.sort(atoms, new EnergyComparator());
            double sum = 0;
            double sum2 = 0;
            double e;
            for (int i = 0; i < atoms.length; i++) {
                e = atoms[i].energy();
                sum += e;
                sum2 += e * e;
            }
            avg = sum / atoms.length;
            std = Math.sqrt(sum2 / atoms.length - avg * avg);
        }

        public double avg() {
            return avg;
        }

        public double std() {
            return std;
        }


        public double sum(double stdThreshold) {
            double sum = 0;
            AtomList lowestEnergyAtoms = lowestEnergyAtoms(stdThreshold);
            Atom atom;
            for (int i = 0; i < lowestEnergyAtoms.size(); i++) {
                atom = lowestEnergyAtoms.atomAt(i);
                sum += atom.energy();
            }
            return sum;
        }

        private AtomList lowestEnergyAtoms(double stdThreshold) {
            AtomList out = new AtomList(atoms[0].molecularSystem);
            int iAtom = 0;
            while ((iAtom < atoms.length) && (atoms[iAtom].energy() < avg + stdThreshold * std)) {
                out.add(atoms[iAtom]);
                iAtom++;
            }
            return out;
        }

        /*private AtomList highestEnergyAtoms(double stdThreshold) {
            AtomList out = new AtomList();
            int iAtom = getAtoms.length - 1;
            while ((iAtom >= 0) & (getAtoms[iAtom].energy() > avg + stdThreshold * std)) {
                out.add(getAtoms[iAtom]);
                iAtom--;
            }
            return out;
        }*/
    }

    public void setCurrent()  {
//        molecularSystem.currentTotalEnergy().evaluate();
        Utils.println("replacing ("+molecularSystem.currentTotalEnergy());
//        Utils.println(molecularSystem.currentTotalEnergy().report(8888));
        molecularSystem.register(this);
        Utils.println("Setting current energy to - " + this);
        on();
        evaluate();
        Utils.println(molecularSystem.currentTotalEnergy().report(9999));
    }

    public MeshiInfo energyInfo()  {
        double energy = evaluateAll();
//        MeshiInfo[] meshiInfos =new MeshiInfo[energyTerms().size()];
        MeshiInfo energyInfo = new MeshiInfo(InfoType.ENERGY,new Double(energy),"Total energy");
 //       int iTerm = 0;
        for (AbstractEnergy term : energyTerms()) {
            energyInfo.getChildren().add(term.energy());
        }
        return energyInfo;
    }

//    public double energyByEnergyInfoType(InfoType infoType) throws UpdateableException, EvaluationException{
//        EnergyInfoElement energyInfoElement;
//        evaluateAll();
//        for (AbstractEnergy term : energyTerms()) {
//            energyInfoElement = new EnergyInfoElement(term.energy());
//            if(energyInfoElement.type == infoType) return energyInfoElement.energy();
//            for (DoubleInfoElement element : energyInfoElement.infoList())
//                if (element.type == infoType) return element.value();
//        }
//        throw new RuntimeException("This is weird. Failed to find "+infoType);
//    }
}



