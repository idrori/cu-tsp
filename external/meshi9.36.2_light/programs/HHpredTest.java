package programs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import meshi.applications.HHpred.Condition;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.GoodBonds;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.geometry.AngleList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.QuickAndDirtyTorsionList;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPairList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.Loop;
import meshi.parameters.ResidueMode;
import meshi.parameters.ResidueType;
import meshi.sequences.AlignmentCell;
import meshi.sequences.AlignmentException;
import meshi.applications.HHpred.HhpredAlignment;
import meshi.applications.HHpred.HhpredAlignmentColumn;
import meshi.util.*;
import meshi.util.file.MeshiWriter;

public class HHpredTest {
    private static boolean debug = true;

	/**
	 * @param args
	 * @throws FileNotFoundException
	 * @throws AlignmentException
	 * @throws UpdateableException
	 * @throws EvaluationException
	 */
	public static void main(String[] args) throws FileNotFoundException,
			UpdateableException, AlignmentException, IOException,
			EvaluationException {

		MeshiProgram.initRandom(7);
		String textfile = args[0]; // hhpred file
		String pdb = args[1]; // pdb file
		String param = args[3];
		CommandList commands = new CommandList(args[2]);
		String output = args[4];
		File file = new File(textfile);
		File pdbfile = new File(pdb);
        BondParametersList bondParametersList = new BondParametersList(EnergyCreator.parametersDirectory(commands) + "/"+ Loop.BOND_PARAMETERS);
		HhpredAlignment hhpredalignment = new HhpredAlignment(file);
		AtomList atoms = new AtomList(pdbfile.getAbsolutePath());
		new MolecularSystem();
		Protein templateProtein = new Protein(atoms,
				ResidueExtendedAtomsCreator.creator);
		Protein queryProtein = new Protein("test");
		MolecularSystem molecularSystem = queryProtein.molecularSystem;
		ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();
		ResidueList residues = new ResidueList();

		buildModel(hhpredalignment, templateProtein, queryProtein, creator,
                molecularSystem, residues, commands, param, output);
		MeshiWriter outWriter = new MeshiWriter("out.pdb");
		queryProtein.atoms().print(outWriter);
		outWriter.close();/**/
	}

    private static Protein buildTrevialModel(HhpredAlignment hhpredalignment,
                                             Protein templateProtein, Protein queryProtein,
                                             ResidueExtendedAtomsCreator creator,
                                             MolecularSystem molecularSystem, ResidueList residues,
                                             CommandList commands) throws IOException{
        MeshiWriter alignmentWriter = null;
        if (debug)
            alignmentWriter = new MeshiWriter("debug\\alignment.pdb");

        StringBuilder comment = new StringBuilder("");
        Residue prevQResidue = new Residue("DMY");
        Condition condition = Condition.Q_RES_AND_TRES_NOT_DUMMY;
        for (HhpredAlignmentColumn hhalcol : hhpredalignment) {
            AlignmentCell tCell = hhalcol.getTCell();
            AlignmentCell qCell = hhalcol.getQCell();
            ResidueType qType = (ResidueType) qCell.obj;
            ResidueType tType = (ResidueType) qCell.obj;

            if (debug) alignmentWriter.print("template-> " + tCell + "  query-> " + qCell + "\n");

            comment = comment.append("template-> ").append(tCell)
                    .append("  query-> ").append(qCell).append("\n");
            System.out.println("template-> " + tCell + "  query-> " + qCell);
            Boolean checkingAlignment = new Boolean(false);
            // ***********************************************************************
            // template protein (from pdb) hhpred template alignment
            // THR_2 5
            // GLN_3 R
            // ARG_5 F
            // PHE_6 T
            // THR_7 E
            // GLU_8 E
            // GLU_9 D
            // ASP_10 F
            // PHE_11
            //
            // ************************************************************************
            prevQResidue = alignResidueAtoms(qType, tType, comment, creator, residues,
                    qCell, templateProtein, commands, molecularSystem,
                    condition, hhalcol, prevQResidue, checkingAlignment, tCell,
                    alignmentWriter);

        }
        alignmentWriter.close();
        Chain chain = new Chain(residues, queryProtein);
        queryProtein.addChain(chain);

        if (debug) {// The trivial model based on the template and alignment
            queryProtein.printAtomsToFile("debug\\ProteinAfteraddChain1.pdb");
            queryProtein.printAtomsToFile("debug\\ProteinAfteraddChain2.pdb",true);
        }
        return queryProtein;
    }

    public static void createLoops(Protein model) {
        for (Chain chain : model.chains()) {
            ResidueExtendedAtoms prevResidue = null;
            for (Residue residue : chain) {
                if ((!residue.dummy()) && (!residue.nowhere())) {
                    if ((prevResidue != null) &&
                            (!prevResidue.dummy()) &&
                            (!prevResidue.nowhere())) {
                        if (residue.number() != prevResidue.number() + 1)
                            throw new RuntimeException("This is weird " + prevResidue + "  " + residue);
                        Atom tail = ((ResidueExtendedAtoms) residue).N;
                        Atom prevHead = prevResidue.C;
                        if (prevHead.distanceFrom(tail) > 3) {
                            Utils.println("A nonphysical peptide bond: " + prevHead.distanceFrom(tail) + " " + prevHead + "  " + tail);
                            residue.setNowhere();
                            prevResidue.setNowhere();
                        }
                    }
                }
                if (!residue.dummy()) prevResidue = (ResidueExtendedAtoms) residue;
            }
        }
    }
	private static void buildModel(HhpredAlignment hhpredalignment,
                                   Protein templateProtein, Protein queryProtein,
                                   ResidueExtendedAtomsCreator creator,
                                   MolecularSystem molecularSystem, ResidueList residues,
                                   CommandList commands, String param, String output)
			throws UpdateableException, AlignmentException, IOException, EvaluationException {

       Protein model = buildTrevialModel(hhpredalignment,templateProtein, queryProtein,creator,
                                              molecularSystem, residues, commands);

        createLoops(model);

        if (debug) {// The trivial model based on the template and alignment
            queryProtein.printAtomsToFile("debug\\ModelAfterRemovingNonPhysical.pdb");
            queryProtein.printAtomsToFile("debug\\ModelAfterRemovingNonPhysical2.pdb",true);
        }

        Loop.addAtoms(model, commands);

        List<Loop> loops = new ArrayList();
		Loop.getLoops(model, loops);
		Loop.addLoopsQuickAndDirty(model, loops, commands, 3, 1000);

        if (debug)
            queryProtein.printAtomsToFile("debug\\afterAddLoops.pdb");

        Loop.addAtoms(queryProtein, commands);

        if (debug)
            queryProtein.printAtomsToFile("debug\\afterAddLoopsAndAddAtoms.pdb");

        System.out.println("*********************** last SCMOD  ************************************");

        Scmod.scmod(commands, model, 2, 50);
        Utils.relax(model, Loop.energyCreatorsBondedTermsEVandRamachandran, commands);

	}


	private static Residue alignResidueAtoms(ResidueType qType, ResidueType tType,
                                             StringBuilder comment, ResidueExtendedAtomsCreator creator,
                                             ResidueList residues, AlignmentCell qCell, Protein templateProtein,
                                             CommandList commands, MolecularSystem molecularSystem,
                                             Condition condition, HhpredAlignmentColumn hhalcol,
                                             Residue prevQResidue, Boolean checkingAlignment,
                                             AlignmentCell tCell, MeshiWriter outWriter) {
		for (Residue tResidue : templateProtein.residues()) {
			if (tResidue.number() < tCell.number) {
				checkingAlignment = true;
				continue;
			} else if (tResidue.number() > tCell.number) {
				if (checkingAlignment == true
						&& residues.get(residues.size() - 1).number() != qCell.number) {
					// problem in template, create and add nowhere residue to
					// protein
					// template protein hhpred template alignment
					// ASP_51_chain:A ASP 51
					// GLU_52_chain:A GLU 52
					// PHE_54_chain:A MET 53
					// PRO_55_chain:A PHE 54
					// ***************************************************************
					createNowhereResidue(residues, qType, qCell, creator,
							comment, molecularSystem, outWriter);
				}
				checkingAlignment = false;
				break;
			}
			if (qType == ResidueType.DMY && tType == ResidueType.DMY)
				break;
			if (condition == Condition.Q_RES_AND_TRES_NOT_DUMMY) {
				condition = Condition.getCondition(hhalcol);
			}
			condition = alignmentOperation(qType, tType, comment, creator,
					residues, qCell, tResidue, commands, molecularSystem,
					condition, hhalcol, prevQResidue, outWriter);
			prevQResidue = residues.get(residues.size() - 1);
		}
		return prevQResidue;

	}

	private static Condition alignmentOperation(ResidueType qType,
			ResidueType tType, StringBuilder comment,
			ResidueExtendedAtomsCreator creator, ResidueList residues,
			AlignmentCell qCell, Residue tResidue, CommandList commands,
			MolecularSystem molecularSystem, Condition condition,
			HhpredAlignmentColumn hhalcol, Residue prevQResidue,
			MeshiWriter outWriter) {
		switch (condition) {
		case Q_RES_AND_TRES_NOT_DUMMY:
			match(qType, tType, comment, creator, residues, qCell, tResidue,
					commands, molecularSystem, outWriter);
			break;
		case Q_RES_DUMMY_AND_TRES_NOT_DUMMY:
			if (qType == ResidueType.DMY) {
				if (!isNowhere(prevQResidue))
					Loop.setNowhere(prevQResidue, outWriter);
			} else if (qType != ResidueType.DMY) {

				setNowhereAndCreateQres(prevQResidue, residues,
						molecularSystem, creator, qType, qCell, outWriter);
				condition = Condition.Q_RES_AND_TRES_NOT_DUMMY;
			}
			break;
		case Q_RES_NOT_DUMMY_AND_TRES_DUMMY:

			if (tType == ResidueType.DMY) {
				setNowhereAndCreateQres(prevQResidue, residues,
						molecularSystem, creator, qType, qCell, outWriter);
			} else if (tType != ResidueType.DMY) {
				setNowhereAndCreateQres(prevQResidue, residues,
						molecularSystem, creator, qType, qCell, outWriter);
				condition = Condition.Q_RES_AND_TRES_NOT_DUMMY;
			}
			break;
		default:
			throw new RuntimeException("this is weird");
		}
		return condition;

	}

	private static void createNowhereResidue(ResidueList residues,
			ResidueType qType, AlignmentCell qCell,
			ResidueExtendedAtomsCreator creator, StringBuilder comment,
			MolecularSystem molecularSystem, MeshiWriter outWriter) {
		if (!exist(qType, qCell.number, residues)) {

			Residue qResidue = creator.create(qType, new ResidueIdentifier(
					qCell.number), ResidueMode.NORMAL, molecularSystem);
			Loop.setNowhere(qResidue, outWriter);
			residues.add(qResidue);

		}

	}

	public static void setNowhereAndCreateQres(Residue prevQResidue,
			ResidueList residues, MolecularSystem molecularSystem,
			ResidueExtendedAtomsCreator creator, ResidueType qType,
			AlignmentCell qCell, MeshiWriter outWriter) {
		if (!isNowhere(prevQResidue))
			Loop.setNowhere(prevQResidue, outWriter);
		Residue qResidue = creator.create(qType, new ResidueIdentifier(
				qCell.number), ResidueMode.NORMAL, molecularSystem);
		Loop.setNowhere(qResidue, outWriter);
		residues.add(qResidue);

	}

    public static boolean isNowhere(Residue residue) {
        for (Atom atom :residue.getAtoms())
            if (!atom.nowhere()) return false;
        return true;
    }

	private static void match(ResidueType qType, ResidueType tType,
			StringBuilder comment, ResidueExtendedAtomsCreator creator,
			ResidueList residues, AlignmentCell qCell, Residue tResidue,
			CommandList commands, MolecularSystem molecularSystem,
			MeshiWriter outWriter) {


		if (!exist(qType, qCell.number, residues)) {

			Residue qResidue = creator.create(qType, new ResidueIdentifier(
					qCell.number), ResidueMode.NORMAL, molecularSystem);
			residues.add(qResidue);
			String atomsCopied = copyCoordinates(tResidue, qResidue, commands,
					outWriter);
			comment = comment.append(atomsCopied);
		}

	}

	private static boolean exist(ResidueType qType, int number,
			ResidueList residues) {// why is qType needed? Tommer 16.11.14
		for (Residue res : residues) {
			if (res.number() == number)
				return true;
		}
		return false;
	}

//	private static void Reconstruct(Protein queryProtein, CommandList commands)
//			throws UpdateableException, AlignmentException, IOException {
//		Utils.addAtoms(queryProtein, true, commands, new PlaneParametersList(
//				EnergyCreator.parametersDirectory(commands) + "/"
//						+ MeshiPotential.PLANE_PARAMETERS));
//
//	}



	public static TorsionList createQuickAndDirtyTorsionList(Residue residue,
			DistanceMatrix distanceMatrix) {
		AtomPairList bondList = (AtomPairList) residue.bonds().filter(
				new GoodBonds());
		// AtomPairList bondList = (AtomPairList) protein.bonds().filter(new
		// GoodBonds());
		AngleList angleList = new AngleList(bondList, distanceMatrix);
		return new QuickAndDirtyTorsionList(angleList, distanceMatrix);
	}

	private static String copyCoordinates(Residue tres, Residue qres,
										  CommandList commands, MeshiWriter outWriter) {
		String content = "";
		int i = 0;
		// System.out.println(tempResidue+ " temp residue in copy Coord");
		for (i = 0; i < tres.getAtoms().size() && i < qres.getAtoms().size(); i++) {
			Atom atom1 = tres.getAtoms().atomAt(i);     // Does it save time?;
			Atom atom2 = qres.getAtoms().atomAt(i);
			if (!atom1.nowhere()) {
				if (tres.name.equals(qres.name))
					setSameCoordinates(atom1, atom2);
				else
					setCoordinates(tres, qres, commands); // Todo if we should no add continue;
				outWriter.print(tres.getAtoms().atomAt(i).toString() + "   "
						+ qres.getAtoms().atomAt(i).toString() + "\n");
				content = content + tres.getAtoms().atomAt(i).toString() + "   "
						+ qres.getAtoms().atomAt(i).toString() + "\n";
				System.out.println(atom1.toString() + "   " + atom2.toString());
			}
		}
		if (tres.getAtoms().size() < qres.getAtoms().size()) {
			for (int j = i; j < qres.getAtoms().size(); j++) {
				outWriter.print(qres.getAtoms().atomAt(j).toString() + "\n");
				content = content
						+ "                                                                                   "
						+ qres.getAtoms().atomAt(j).toString() + "\n";
				System.out
						.println("                                                                                   "
								+ qres.getAtoms().atomAt(j).toString());
			}
		} else {
			for (int j = i; j < tres.getAtoms().size(); j++) {
				outWriter.print(tres.getAtoms().atomAt(j).toString() + "\n");
				content = content + tres.getAtoms().atomAt(j).toString() + "\n";
				System.out.println(tres.getAtoms().atomAt(j).toString());
			}
		}
		return content;

	}

	               //Todo may be replace with a black list
	private static boolean allowedCoordCopy(Residue tres, Atom tatom,
			Residue qres, Atom qatom) {
		if (qres.name.equals("PRO") && (qatom.name().equals("CD")))
			return false;// was already present (22.11.14)
		if (qres.name.equals("PRO") && (qatom.name().equals("CG")))
			return false;// was already present (22.11.14)
		if (tres.name.equals("PRO") && (tatom.name().equals("CD")))
			return false;// to take care of
		if (tres.name.equals("PRO") && (tatom.name().equals("CG")))
			return false;// the other direction.
		if (qres.name.equals("ARG") && (qatom.name().equals("CZ")))
			return false;// these two line should replace
		if (tres.name.equals("ARG") && (tatom.name().equals("CZ")))
			return false;// the handling of CZ.
        if (qres.name.equals("TRP") && (qatom.name().startsWith("CE")))
            return false;// these two line should replace
        if (qres.name.equals("TRP") && (qatom.name().equals("NE1")))
            return false;// these two line should replace
        if (qres.name.equals("TRP") && (qatom.name().startsWith("CZ")))
            return false;// these two line should replace
        if (qres.name.equals("TRP") && (qatom.name().startsWith("CH")))
            return false;// these two line should replace
        if (tres.name.equals("TRP") && (qatom.name().startsWith("CE")))
            return false;// these two line should replace
        if (tres.name.equals("TRP") && (qatom.name().equals("NE1")))
            return false;// these two line should replace
        if (tres.name.equals("TRP") && (qatom.name().startsWith("CZ")))
            return false;// these two line should replace
        if (tres.name.equals("TRP") && (qatom.name().startsWith("CH")))
            return false;// these two line should replace
		// if(qatom.name.equals("CZ")) return false; Tommer 22/11/14
		return true;
	}

	private static void setCoordinates(Residue tResidue, Residue qResidue,
			CommandList commands) {
		for (Atom tAtom : tResidue.getAtoms()) {
			// if(!(tatom.isCarbon() || tatom.i)) continue;
			for (Atom qAtom : qResidue.getAtoms()) {
				if (qAtom.name().equals(tAtom.name()) && !tAtom.nowhere()) {
					if (allowedCoordCopy(tResidue, tAtom, qResidue, qAtom))
						setSameCoordinates(tAtom, qAtom);
				}
			}
		}

		for (Atom qatom : qResidue.getAtoms()) {
			boolean found = false;
			if (!qatom.nowhere()) {
				for (Atom neighbor : qatom.bonded())
					if (!neighbor.nowhere())
						found = true;
				if (!found) {
					qatom.setNowhere();
				}
			}
		}
	}

	private static void setSameCoordinates(Atom atom1, Atom atom2) {
		atom2.setXYZ(atom1.x(), atom1.y(), atom1.z());
	}

}