package meshi.molecularElements.loops;

import java.io.IOException;
import java.util.*;

import meshi.applications.prediction.OriginalAtom;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.Coordinates;
import meshi.geometry.DistanceMatrix.DistanceMatrixType;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.CaFilter;


public class Loop extends ResidueList {
    public static void setNowhere(Residue tempQResidue, MeshiWriter outWriter) {
        for (Atom atom : tempQResidue.getAtoms()) {
            atom.setNowhere();
        }
        if (outWriter != null) outWriter.print(tempQResidue.name + " " + tempQResidue.number()
                                                + "  is nowhere\n");
        System.out.println(tempQResidue.name + " " + tempQResidue.number()
                + "  is nowhere");
    }

    public enum LoopMode {NTERM,LOOP,CTERM}
    public static final double CA_CA_DISTANCE = 3.8;

    public final static String BOND_PARAMETERS = "meshiPotential/bondEnergyParameters.dat";
    public final static String ANGLE_PARAMETERS = "meshiPotential/angleEnergyParameters.dat";
    public final static String PLANE_PARAMETERS = "meshiPotential/planeEnergyParameters.dat";
    private int rank = -1;
    private static boolean debug = true;
    private static int nAddAtoms = 0;
	/*
	 * [ ][ ][breakStart]--------[breakEnd][ ][ ] 'breakStart' is the residue
	 * number of the last residue before the break. 'breakEnd' is the residue
	 * number of the first residue after the break.
	 */
	public static final double[][] data = {
			{ 0.0000, 5.9000, 6.0000, 6.1000, 6.1500 },
			{ 1.0000, 5.9000, 6.0000, 6.1000, 6.1500 },
			{ 2.0000, 5.9000, 6.0000, 6.1000, 6.1500 },
			{ 3.0000, 9.0000, 9.4000, 9.6000, 9.7000 },
			{ 4.0000, 12.0000, 12.7000, 13.0000, 13.1000 },
			{ 5.0000, 14.9000, 15.9000, 16.3000, 16.5000 },
			{ 6.0000, 17.5000, 19.0000, 19.6000, 19.9000 },
			{ 7.0000, 19.7000, 22.1000, 22.8000, 23.1000 },
			{ 8.0000, 21.7000, 24.9000, 25.8500, 26.4000 },
			{ 9.0000, 23.4000, 27.7000, 28.9000, 29.4000 },
			{ 10.0000, 24.8000, 30.1000, 31.7000, 32.4000 },
			{ 11.0000, 25.9000, 32.2000, 34.5000, 35.4500 },
			{ 12.0000, 26.8000, 34.1000, 37.1500, 38.2000 },
			{ 13.0000, 27.4000, 35.7500, 39.5000, 40.9000 },
			{ 14.0000, 28.1000, 37.1000, 41.5000, 43.3000 },
			{ 15.0000, 28.6000, 38.0000, 43.0500, 45.6000 },
			{ 16.0000, 29.0000, 39.0000, 44.5000, 48.4000 },
			{ 17.0000, 29.5000, 39.7500, 45.9000, 50.4500 },
			{ 18.0000, 29.9000, 40.4000, 47.0000, 52.8000 },
			{ 19.0000, 30.4000, 41.1000, 47.9000, 54.2000 } };

    public final static EnergyCreator[] energyCreatorsBondedTermsEVandRamachandran = {
            new BondCreator(),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new AtomicPairwisePMFSummaCreator(),
            new RamachandranSidechainEnergyCreator()
    };


    public final static	EnergyCreator[] energyCreatorsBondedTermsPlusEV = {
            new BondCreator(),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new AtomicPairwisePMFSummaCreator(),
            new RamachandranSidechainEnergyCreator()
    };

    public int rank() {
		if (rank == -1) {
			Residue first = get(0);
			Residue last = get(size() - 1);
			rank = last.number() - first.number()
					+ MeshiProgram.randomNumberGenerator().nextInt(3);
		}
		return rank;
	}

	public static void deletion(Iterator residues, List<Loop> tempList,
			Loop loop, Three three) {
		if (!residues.hasNext()) {
			tempList.add(loop);
			return;
		}
	step(residues, three);
	three.getFirst().addAttribute(LoopAttribute.attribute);
	loop.add(three.getFirst());
	if (three.getFirst().ca().nowhere())
			insertion(residues, tempList, loop, three);
		else {
			tempList.add(loop);
			getLoopsBefore(residues, tempList, three);
		}
	}

    private static void step(Iterator residues, Three three) {
        while ((three.getFirst() == null) || (three.getFirst().dummy())) three.setFirst((Residue) residues.next());
        three.setThird(three.getSecond());
        three.setSecond(three.getFirst());
        three.setFirst((Residue) residues.next());
    }

    private static class LoopAttribute implements MeshiAttribute {
        public static LoopAttribute attribute = new LoopAttribute();
        public int key() {return LOOP_RESIDUE;}
    }

	public static void insertion(Iterator residues, List<Loop> tempList, Loop loop, Three three) {
		if (!residues.hasNext()) {
			tempList.add(loop);
			return;
		}
		step(residues, three);
	three.getFirst().addAttribute(LoopAttribute.attribute);
	loop.add(three.getFirst());
	if (three.getFirst().ca().nowhere())
			insertion(residues, tempList, loop, three);
		else
			Loop.deletion(residues, tempList, loop, three);
	}

    private static void getLoopsBefore(Iterator residues,  List<Loop> tempList, Three three) {
        if (!residues.hasNext()) return;
        step(residues, three);
        if (three.getFirst().ca().nowhere()) {
            Loop loop = new Loop();
            if (three.getThird()  != null) {
                three.getThird().addAttribute(LoopAttribute.attribute);
                loop.add(three.getThird());
            }
            if (three.getSecond() == null) throw new RuntimeException("This is weird");
            three.getSecond().addAttribute(LoopAttribute.attribute);
            loop.add(three.getSecond());
            three.getFirst().addAttribute(LoopAttribute.attribute);
            loop.add(three.getFirst());
            insertion(residues, tempList, loop, three);
        }
        if ((three.getSecond() != null) && (!three.getSecond().ca().nowhere()))
            if (three.getFirst().ca().distanceFrom(three.getSecond().ca()) > 4.5) {
                Loop loop = new Loop();
                if (three.getThird()  != null) {
                    three.getThird().addAttribute(LoopAttribute.attribute);
                    loop.add(three.getThird());
                }
                if (three.getSecond() == null) throw new RuntimeException("This is weird");
                three.getSecond().addAttribute(LoopAttribute.attribute);
                loop.add(three.getSecond());
                three.getFirst().addAttribute(LoopAttribute.attribute);
                loop.add(three.getFirst());
                deletion(residues, tempList, loop,three);
            }
        getLoopsBefore(residues, tempList, three);
    }
    public static void getLoops(Protein model, List<Loop> loops){
        Three three = new Three();
        List<Loop> tempList = new ArrayList();
        Iterator residues = model.chain().iterator();
        getLoopsBefore(residues, tempList, three);
        Object[] tempArray = tempList.toArray();
        Arrays.sort(tempArray,new LoopComperator());
        for (Object o:tempArray) {
            Loop loop = (Loop) o;
            System.out.println(loop);
            for (residues = loop.iterator(); residues.hasNext();) {
                Residue residue = (Residue) residues.next();
            }
            loops.add((Loop) o);
        }
    }

    private static class LoopComperator implements Comparator {
        public int compare(Object a, Object b) {
            Loop l1 = (Loop) a;
            Loop l2 = (Loop) b;
            if (l1.rank()>l2.rank()) return 1;
            if (l1.rank()<l2.rank()) return -1;
            return 0;
        }
        public boolean equals(Object a, Object b) {
            return (compare(a,b) == 0);
        }
    }


//
//    public static void getLoopsBefore(Iterator residues, List<Loop> tempList,
//			HomologyModeling.Three three) {
//		if (!residues.hasNext())
//			return;
//
//		HomologyModeling.step(residues, three);
//		if (three.first.ca().nowhere()) {
//			Loop loop = new Loop();
//			if (three.third != null) {
//				three.third
//						.addAttribute(HomologyModeling.LoopAttribute.attribute);
//				loop.add(three.third);
//			}
//			if (three.second == null)
//				throw new RuntimeException("This is weird");
//			three.second.addAttribute(HomologyModeling.LoopAttribute.attribute);
//			loop.add(three.second);
//			three.first.addAttribute(HomologyModeling.LoopAttribute.attribute);
//			loop.add(three.first);
//			Loop.insertion(residues, tempList, loop, three);
//		}
//		if ((three.second != null) && (!three.second.ca().nowhere()))
//			if (three.first.ca().distanceFrom(three.second.ca()) > 4.5) {
//				Loop loop = new Loop();
//				if (three.third != null) {
//					three.third
//							.addAttribute(HomologyModeling.LoopAttribute.attribute);
//					loop.add(three.third);
//				}
//				if (three.second == null)
//					throw new RuntimeException("This is weird");
//				three.second
//						.addAttribute(HomologyModeling.LoopAttribute.attribute);
//				loop.add(three.second);
//				three.first
//						.addAttribute(HomologyModeling.LoopAttribute.attribute);
//				loop.add(three.first);
//				Loop.deletion(residues, tempList, loop, three);
//			}
//		getLoopsBefore(residues, tempList, three);
//	}

//	public static void getLoops(Protein model, List<Loop> loops,
//			MeshiWriter outWriter) {
//		HomologyModeling.Three three = new HomologyModeling.Three();
//		List<Loop> tempList = new ArrayList<Loop>();
//		Iterator<Residue> residues = model.chain().getNonDummyResidues()
//				.iterator();
//		// Iterator<Residue> residues = model.chain().iterator();
//		model.chain().getNonDummyResidues();
//		Loop.getLoopsBefore(residues, tempList, three);
//		Object[] tempArray = tempList.toAlignmentArray();
//		Arrays.sort(tempArray, new HomologyModeling.LoopComperator());
//		for (Object o : tempArray) {
//			Loop loop = (Loop) o;
//			outWriter.print(loop.toString() + "\n");
//			System.out.println(loop.toString());
//			for (residues = loop.iterator(); residues.hasNext();) {
//				Residue residue = (Residue) residues.next();
//			}
//			loops.add((Loop) o);
//		}
//	}
//
	@Override
	public String toString() {
		String cont = "LOOP: ";
		for (Residue res : this) {
			cont = cont + res.toString() + "	";
		}
		return cont;
	}

    public static boolean setCAsQD(Protein model, Loop loop, CommandList commands, double clashDistance, int nTrys) {
        Residue first =  loop.get(0);
        Residue last  =  loop.get(loop.size()-1);

        int to = -1,from = -1, direction = -1;
        Atom toAtom, fromAtom, currentAtom;
        double rnd = MeshiProgram.randomNumberGenerator().nextDouble();
        if (rnd < 0.5) {
                from = 0;
                to = loop.size()-1;
                direction = 1;
        }
        else {
                from = loop.size()-1;
                to = 0;
                direction = -1;
        }
        fromAtom = loop.get(from).ca();
        toAtom = loop.get(to).ca();

        double distanceGrade;
        Coordinates save = null;
        System.out.println("Doing loop "+loop);
        System.out.println("first = "+first+" ; last = "+last+" ; direction = "+direction);
        System.out.println("fromAtom = "+fromAtom+" ; toAtom "+toAtom);
        for (int current = from+direction; current != to; current += direction) {
            Residue residue = loop.get(current);
            currentAtom = residue.ca();
            if (currentAtom.nowhere()) {
                System.out.println("Now working isOn "+current);
                double bestGrade = 100000;
                List<Atom> neighbors = Utils.getNeighbors(fromAtom, model, clashDistance);
                for (int i = 0; i < nTrys; i++) {
                    Coordinates temp = new Coordinates(new Coordinates(fromAtom.x(), fromAtom.y(), fromAtom.z()), CA_CA_DISTANCE);
                    currentAtom.setXYZ(temp);
                    int nClashes    = Utils.getClashes(currentAtom,neighbors, clashDistance);
                    distanceGrade = getDistanceGrade(currentAtom,toAtom);
                    double grade = nClashes+distanceGrade;
                    if (grade < bestGrade) {
                        System.out.println("Try No. "+i+" "+currentAtom.residue()+" "+toAtom.residue()+
                                    " distance "+currentAtom.distanceFrom(toAtom)+" clashes "+nClashes+" grade "+grade);
                        save = temp;
                        bestGrade = grade;
                    }
                }
                currentAtom.setXYZ(save.x(),save.y(),save.z());
            }
            fromAtom = currentAtom;
        }
        Atom prevAtom = null;
        for (Residue residue : loop) {
            Atom atom = residue.ca();
            if (prevAtom != null)
                if (prevAtom.distanceFrom(atom) > CA_CA_DISTANCE * 1.1) {
                    for (int i = 1; i < loop.size()-1; i++)
                        loop.get(i).setNowhere();
                    return false;
                }
        }
        return true;
    }

    /*
   * [   ][   ][breakStart]--------[breakEnd][   ][   ]
   * 'breakStart' is the residue number of the last residue before the break.
   * 'breakEnd' is the residue number of the first residue after the break.
   */
    private static double getDistanceGrade(Atom start, Atom end) {
        int len = end.residue().number()-start.residue().number();
        if (len < 0) len = -1*len;
        double distance = start.distanceFrom(end) - 2.4;  // The 2.4 is the length of two peptide bonds
        if (len >= data.length) return 0;
        if (distance>data[len][4])
            return 10.0;
        if (distance>data[len][3])
            return 4.3 + (distance-data[len][3])/(data[len][4]-data[len][3])*2.3;
        if (distance>data[len][2])
            return 2.3 + (distance-data[len][2])/(data[len][3]-data[len][2])*2.3;
        if (distance>data[len][1])
            return (distance-data[len][1])/(data[len][2]-data[len][1])*2.3;
        return 0.0;
    }

//
//    public static void setCAsQD(Protein model, Loop loop, CommandList commands,
//			double clashDistance, int nTrys) {
//		Residue first = (Residue) loop.get(0);
//		Residue last = (Residue) loop.get(loop.size() - 1);
//		HomologyModeling.LoopMode mode = null;
//		if (first.ca().nowhere())
//			mode = HomologyModeling.LoopMode.NTERM;
//		else if (last.ca().nowhere())
//			mode = HomologyModeling.LoopMode.CTERM;
//		else if ((!first.ca().nowhere()) && (!last.ca().nowhere()))
//			mode = HomologyModeling.LoopMode.LOOP;
//		else
//			throw new RuntimeException("Weird loop " + loop);
//
//		int to = -1, from = -1, direction = -1;
//		Atom toAtom, fromAtom, currentAtom;
//		double rnd = MeshiProgram.randomNumberGenerator().nextDouble();
//		if (mode == HomologyModeling.LoopMode.NTERM) {
//			from = loop.size() - 1;
//			to = -1;
//			direction = -1;
//		}
//		if (mode == HomologyModeling.LoopMode.CTERM) {
//			from = 0;
//			to = loop.size();
//			direction = 1;
//		}
//		if (mode == HomologyModeling.LoopMode.LOOP) {
//			if (rnd < 0.5) {
//				from = 0;
//				to = loop.size() - 1;
//				direction = 1;
//			} else {
//				from = loop.size() - 1;
//				to = 0;
//				direction = -1;
//			}
//		}
//		fromAtom = loop.get(from).ca();
//		if (mode == HomologyModeling.LoopMode.LOOP)
//			toAtom = loop.get(to).ca();
//		else
//			toAtom = null;
//
//		double distanceGrade;
//		Coordinates save = null;
//		System.out.println("Doing loop " + loop);
//		System.out.println("first = " + first + " ; last = " + last
//				+ " ; direction = " + direction + " ; mode = " + mode);
//		System.out.println("fromAtom = " + fromAtom + " ; toAtom " + toAtom);
//		for (int current = from + direction; current != to; current += direction) {
//			Residue residue = loop.get(current);
//			currentAtom = residue.ca();
//			if (currentAtom.nowhere()) {
//				System.out.println("Now working isOn " + current);
//				double bestGrade = 100000;
//				List<Atom> neighbors = Utils.getNeighbors(fromAtom, model,
//						clashDistance);
//				for (int i = 0; i < nTrys; i++) {
//					Coordinates temp = new Coordinates(new Coordinates(
//							fromAtom.x(), fromAtom.y(), fromAtom.z()),
//							HomologyModeling.CA_CA_DISTANCE);
//					currentAtom.setXYZ(temp);
//					int nClashes = Utils.getClashes(currentAtom, neighbors,
//							clashDistance);
//					if (toAtom != null)
//						distanceGrade = HomologyModeling.getDistanceGrade(
//								currentAtom, toAtom);
//					else
//						distanceGrade = 0;
//					double grade = 10 * nClashes + distanceGrade;
//					if (grade < bestGrade) {
//						if (toAtom != null)
//							System.out.println("Try No. " + i + " "
//									+ currentAtom.residue() + " "
//									+ toAtom.residue() + " distance "
//									+ currentAtom.distanceFrom(toAtom)
//									+ " clashes " + nClashes + " grade "
//									+ grade);
//						else
//							System.out.println("Try No. " + i + " "
//									+ currentAtom.residue() + " clashes "
//									+ nClashes + " grade " + grade);
//						save = temp;
//						bestGrade = grade;
//					}
//				}
//				currentAtom.setXYZ(save.x(), save.y(), save.z());
//			}
//			fromAtom = currentAtom;
//		}
//	}

	public static void setHOsQD(Protein model, Loop loop, CommandList commands)
			throws UpdateableException, EvaluationException, AlignmentException {
		boolean success = false;
		// for (int i = 0; (i<5) & (!success); i++) {
		model.atoms().freeze("In HomologyModeling");
		for (Residue residue : loop) {
			Atom amideN = residue.amideN();
			Atom carbonylC = residue.carbonylC();
			Atom amideH = residue.amideH();
			Atom carbonylO = residue.carbonylO();
			Coordinates amidNcoor = new Coordinates(amideN.x(), amideN.y(),
					amideN.z());
			Coordinates carbonylCcoor = new Coordinates(carbonylC.x(),
					carbonylC.y(), carbonylC.z());
			if ((amideH != null) && amideH.nowhere())
				amideH.setXYZ(new Coordinates(amidNcoor, 1.5));
			if (carbonylO.nowhere())
				carbonylO.setXYZ(new Coordinates(carbonylCcoor, 1.5));
		}
		AtomList neighbors = Utils.getNeighbors(loop, model, 6);
		System.out.println("****************** Third minimization *****************************");
		//if (1 == 1)     		throw new RuntimeException("???????????????????????setNCCBsQD???????????????????????????????");
        model.molecularSystem.createDistanceMatrix("This ie a risky stage in loop formation.",DistanceMatrixType.STANDARD);
		TotalEnergy te = Utils.relax(model,
                energyCreatorsBondedTermsPlusEV,
				                     commands);
		if (te.evaluate() < 10000) {
			success = true;
		} else
			System.out.println("Failed to relax setHOsQD");

		// System.out.println("Problem # "+i+" in  setHOsQD");
		/*
		 * System.out.println(ex); ex.printStackTrace(); //if (i >= 5) throw ex;
		 * }
		 */
		// }
	}

	public static void setNCCBsQD(Protein model, Loop loop, CommandList commands)
			throws UpdateableException, EvaluationException, AlignmentException {
		model.atoms().freeze("In HomologyModeling");
		for (Object aLoop : loop) {
			Residue residue = (Residue) aLoop;
			Atom ca = residue.ca();
			Atom amideN = residue.amideN();
			Atom carbonylC = residue.carbonylC();
			Atom cb = residue.cb();
			Coordinates caCoor = new Coordinates(ca.x(), ca.y(), ca.z());
			if (amideN.nowhere())
				amideN.setXYZ(new Coordinates(caCoor, 1.5));
			if (carbonylC.nowhere())
				carbonylC.setXYZ(new Coordinates(caCoor, 1.5));
			if ((cb != null) && cb.nowhere())
				cb.setXYZ(new Coordinates(caCoor, 1.5));
		}
		AtomList neighbors = Utils.getNeighbors(loop, model, 6);
		neighbors.defrost();
		System.out
				.println("****************** First minimization *****************************");
		try {
			model.atoms().molecularSystem().getDistanceMatrix();
			// model.molecularSystem.createDistanceMatrix();
			Utils.relax(model,
					    energyCreatorsBondedTermsPlusEV, commands);
		} catch (RuntimeException ex) {
			System.out.println(ex);
			System.out.println("Problem in setNCCBsQD working isOn loop "
					+ loop + "\n" + "neighbor getAtoms are:");
			neighbors.print();
			throw ex;
		}
		model.atoms().defrost();
		model.atoms().freeze(OriginalAtom.filter);
		neighbors = Utils.getNeighbors(loop, model, 6);
		System.out
				.println("****************** Second minimization *****************************");
		Utils.relax(model, energyCreatorsBondedTermsPlusEV,commands);
	}

	public static void addLoopQD(Protein model, Loop loop,
			CommandList commands, double clashDistance, int nTrys)
			throws UpdateableException, EvaluationException, AlignmentException {
		boolean success = false;
		RuntimeException ex = null;
		// for (int i = 0; (i <5) & (!success); i++) {

		Loop.setCAsQD(model, loop, commands, clashDistance, nTrys);
		Loop.setNCCBsQD(model, loop, commands);
		Loop.setHOsQD(model, loop, commands);
		success = true;
		if (!success)
			throw new RuntimeException(ex);
	}

    public static void addLoopsQuickAndDirty(Protein model, List<Loop> loops,
                                             CommandList commands, double clashDistance, int nTrys)
                                             throws UpdateableException, EvaluationException,
                                                    AlignmentException, IOException {
        addLoopsQuickAndDirty(model,loops,commands,clashDistance,nTrys,10000);
    }

    public static void addLoopsQuickAndDirty(Protein model, List<Loop> loops,
			CommandList commands, double clashDistance, int nTrys, int numberOfLoops)
			throws UpdateableException, EvaluationException,
			AlignmentException, IOException {
		boolean success = false;
		RuntimeException ex = null;

		int count=0; //Tommer 5.2.15
		long loopsTime= System.currentTimeMillis();
		long oneLoopTime;
		for (Loop loop : loops) {
			count++;
            if (count > numberOfLoops) return;
			oneLoopTime= System.currentTimeMillis();
			try{
                Loop.addLoopQD(model, loop, commands, clashDistance, nTrys);
                if (debug) {
                    model.printAtomsToFile("debug\\afterAddLoop" + count + "_1.pdb");
                    model.printAtomsToFile("debug\\afterAddLoop" + count + "_2.pdb", true);
                }
            }
			catch(Exception exe){
				exe.printStackTrace();
				throw new RuntimeException("problem in loop number "+count+
						"\ntime on this loop: "+(System.currentTimeMillis()-oneLoopTime)+
						"\ntime on loops so far: "+(System.currentTimeMillis()-loopsTime));
			}
		}
		System.out
				.println("***********************  backbones of loops added  ************************************");
		model.atoms().defrost();
		// model.getAtoms().freeze(OriginalAtom.filter);
		//addAtoms(model, commands);
		MeshiWriter outWriter = new MeshiWriter("modelAfteraddAtoms.pdb");
		outWriter.print(model);
		outWriter.close();
		// Add the rest of the sidechains
		model.atoms().freeze(new CaFilter());
		Utils.relax(model,energyCreatorsBondedTermsEVandRamachandran,commands);


	}


	public static void addAtoms(Protein model, CommandList commands) throws EvaluationException, UpdateableException, AlignmentException, IOException {

        EnergyCreator[] energyCreatorsBondedTermsOnly = {
                new BondCreator(), new AngleCreator(), new PlaneCreator(),
                new OutOfPlaneCreator(),
                new RamachandranSidechainEnergyCreator()};
                EnergyCreator[] energyCreatorsBondAngleOnly = {
                        new BondCreator(), new AngleCreator()};

        Command command = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        BondParametersList  bondParametersList  = new BondParametersList(parametersDirectory+"/" + BOND_PARAMETERS);
        AngleParametersList angleParametersList = new AngleParametersList(parametersDirectory+"/" + ANGLE_PARAMETERS);
        PlaneParametersList planeParametersList = new PlaneParametersList(parametersDirectory+"/" + PLANE_PARAMETERS);
        for (Atom atom : model.atoms()) {
            if ((!atom.nowhere()) && (atom.x() < -999))
             atom.setNowhere();
        }

        for (Atom atom : model.atoms()) {
            if (atom.isHydrogen() && atom.nowhere()) {
                PutHposLog puthLog = PutHpos.pos(atom, bondParametersList, angleParametersList);
            }
        }
        if (debug) {
            model.printAtomsToFile("debug\\afterAddingHydrogens1" + nAddAtoms + ".pdb");
            model.printAtomsToFile("debug\\afterAddingHydrogens2" + nAddAtoms + ".pdb",true);
        }

        AtomFinding.finalFindAtoms(model, bondParametersList, angleParametersList, planeParametersList);
        Scmod.scmod(commands,model,2,true,MeshiProgram.randomNumberGenerator());

        if (debug) {
            model.printAtomsToFile("debug\\afterAddingAtoms"+nAddAtoms+"_1.pdb");
            model.printAtomsToFile("debug\\afterAddingAtoms"+nAddAtoms+"_2.pdb",true);
        }

        double max = 1000;
        int counter = 0;
        while (max > 20) {
            counter++;
            TotalEnergy energy = Utils.relax(model, energyCreatorsBondAngleOnly, commands);
            energy.evaluateAtoms();
            max = -1000;
            Atom maxAtom = null;
            for (Atom atom : model.atoms()) {
                if ((!atom.nowhere()) && (atom.energy()) > max) {
                    max = atom.energy();
                    maxAtom = atom;
                }
            }
            if (max > 20) {
                if (debug)
                    Utils.printDebug("Loop.addAtoms", " sent " + maxAtom.residue() + " to nowhere, due to " + maxAtom + " of energy " + max);
                setNowhere(maxAtom.residue(), null);
            }
        }
//        for (Atom iAtom : model.getAtoms()) {
//            if (iAtom.nowhere()) {
//                boolean found = false;
//                for (Atom jAtom : iAtom.residue().getAtoms())
//                    if (!jAtom.nowhere()) found = true;
//                if (found)
//                    for (Atom neighbor :iAtom.bonded())
//                        if (!neighbor.nowhere())   iAtom.setXYZ(new Coordinates(new Coordinates(neighbor),1.5));
//                }
//        }
        Utils.relax(model, energyCreatorsBondedTermsEVandRamachandran, commands);
//        Scmod.scmod(commands, model, 2, 50);

        if (debug) {
            model.printAtomsToFile("debug\\modelAfterRelax"+counter+"_1.pdb");
            model.printAtomsToFile("debug\\modelAfterRelax"+counter+"_2.pdb",true);
        }
    }
}
