package meshi.applications.conserv.statistic;

import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.molecularElements.ConSeq;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.AtomType;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.info.InfoType;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**

 */
public class ConservationStatistics {

	private static final String 	RMS 				= "RMS";
	private static final String 	GDT 				= "GDT";

	private static Protein 			mNativeStructure 	= null;
	private static int 				mScoreType 			= 0; //ConservationContactsRatio --> 0, ConservationContactsHrRatio --> 1
	private static Vector<Double> 	mRadius; //Default for ConservationContactsRatio: 8, Default for ConservationContactsHrRatio: 6
	private static String 			mDirName;
	private static OutputHandler	sRMSConOutputHandler; // = new OutputHandler("HrResults.csv", headers);
	private static OutputHandler	sRMSRgOutputHandler;
	private static OutputHandler	sGDTConOutputHandler; // = new OutputHandler("HrResults.csv", headers);
	private static OutputHandler	sGDTRgOutputHandler;
	private static String 			mXValue;

	/**
	 * 
	 * @param argv
	 * @throws IOException
	 */
	public static void main(String[] argv) throws IOException,AlignmentException {

		/*Add to Run -> Configuration -> Arguments :
			CASP9_FM_decoys execOrders.txt

		execOrders.txt contains:
		<output file name> <ScoreType (0 for regular and 1 for HR)> <X_label value (o for RMS and 1 for GDT)> <list of radius values>*/

		mDirName = argv[0];
//        File dir = new File(mDirName);
//        System.out.println(dir.getAbsolutePath());
//        String[] files = dir.list();
//        for (String fn : files){
//            System.out.println(fn);
//        }
		MeshiProgram.initRandom(0);
		BufferedReader in = new BufferedReader(new FileReader(argv[1]));
		String str;
		while ((str = in.readLine()) != null) {
			String[] params = str.split(" ");
			mScoreType = Integer.valueOf(params[1]);
			if (Integer.valueOf(params[2]) == 0) {
				mXValue = RMS;
			} else {
				mXValue = GDT;
			}
			mRadius = new Vector<Double>();
			List<String> headers = new ArrayList<String>();
			for (int i=3; i<params.length; i++) {
				mRadius.add(Double.valueOf(params[i]));
				headers.add("radius " + mRadius.elementAt(i-3));
			}

			sRMSConOutputHandler = new OutputHandler("con"+RMS+params[0], headers, RMS);
			sRMSRgOutputHandler = new OutputHandler("rg"+RMS+params[0], headers, RMS);
			sGDTConOutputHandler = new OutputHandler("con"+GDT+params[0], headers, GDT);
			sGDTRgOutputHandler = new OutputHandler("rg"+GDT+params[0], headers, GDT);
			runCalculating();
		}
		in.close();
		sGDTConOutputHandler.close();
		sGDTRgOutputHandler.close();
		sRMSConOutputHandler.close();
		sRMSRgOutputHandler.close();


	}

	private static void runCalculating() throws IOException,AlignmentException {
		File  dir = new File(mDirName);
		File[] files = dir.listFiles();
		String[] args = {"commands.txt"};
		CommandList commands = new CommandList(args,
				new CommandsException("This is weird"));
		int i = 0;
		for (File file: files) {
			i++;
			String fileName = file.getAbsolutePath();
			if (fileName.endsWith("pdb")) {
				System.out.println(fileName+ " #" + i + " out of estimated " + files.length);
				BondParametersList bondParametersList = Utils.getBondParameters(commands);
				AngleParametersList angleParametersList = Utils.getAngleParameters(commands);
				Protein protein = Protein.getExtendedAtomsProteinFromPdbFile(file, bondParametersList, angleParametersList);
				    getNative(file);

				if (mScoreType==0) {
					getConservationContactsRatioScores(protein, fileName);
				} else if (mScoreType == 1) {
                    if (!hasNowhwerAtoms(protein))  {
					      getConservationContactsHrRatioScores(protein, fileName);
                    }
                    else {
                        System.out.println("Ignoring "+fileName+" incomplete model.");
                    }
				} else {
					System.out.printf("*** Unknows Score type: " + mScoreType);
				}
			}
		}
	}

    public static boolean hasNowhwerAtoms(Protein protein) {
        AtomList atoms = protein.atoms();
        for (Atom atom: atoms)
            if (atom.nowhere()&& (atom.type() != AtomType.TRO))
                return true;
        return false;
    }
	private static void getConservationContactsRatioScores(Protein protein, String fileName) throws IOException,AlignmentException {
		ConservationContactsRatio score = getConservationContactsRatio(protein,fileName);
		double rms = getRms(mNativeStructure,protein);
		double[] gdt = getGdt(mNativeStructure,protein);
		Map<String, Double> conResults = new HashMap<String, Double>();
		Map<String, Double> rgResult = new HashMap<String, Double>();
		for (Double radius : mRadius) {
			score.setRadius(radius.doubleValue());
			EnergyInfoElement element = score.evaluate();
			conResults.put("radius " + radius.doubleValue(), new Double(element.energy()));
			rgResult.put("radius " + radius.doubleValue(), new Double((Double) element.getChildren().get(1).getValue()).doubleValue());
		}

		//System.out.printf("%-30f %-30f %-30f %-30f\n", contactRation,rgRatio, rms,gdt[0]);
		//if (mXValue == RMS) {
		sRMSConOutputHandler.writeLine(rms, conResults, fileName);
		sRMSRgOutputHandler.writeLine(rms, rgResult, fileName);
		//} else {
		sGDTConOutputHandler.writeLine(gdt[0], conResults, fileName);
		sGDTRgOutputHandler.writeLine(gdt[0], rgResult, fileName);
		//}
	}

	private static void getConservationContactsHrRatioScores(Protein protein, String fileName) throws IOException,AlignmentException {
		ConservationContactsHrRatio score = getConservationContactsHrRatio(protein,fileName);
		double rms = getRms(mNativeStructure,protein);
		double[] gdt = getGdt(mNativeStructure,protein);
		Map<String, Double> conResults = new HashMap<String, Double>();
		Map<String, Double> rgResult = new HashMap<String, Double>();
		for (Double radius : mRadius) {
			score.setRadius(radius.doubleValue());
			EnergyInfoElement element = score.evaluate();
			conResults.put("radius " + radius.doubleValue(), new Double(element.energy()));
			rgResult.put("radius " + radius.doubleValue(), new Double((Double) element.getChildren().get(1).getValue()).doubleValue());
		}

		//System.out.printf("%-30f %-30f %-30f %-30f\n", contactRation,rgRatio, rms,gdt[0]);
		sRMSConOutputHandler.writeLine(rms, conResults, fileName);
		sRMSRgOutputHandler.writeLine(rms, rgResult, fileName);
		//} else {
		sGDTConOutputHandler.writeLine(gdt[0], conResults, fileName);
		sGDTRgOutputHandler.writeLine(gdt[0], rgResult, fileName);
		//}
	}

	/*private static String getName(String fileName) {
		int loc = fileName.indexOf(mDirName);
		return fileName.substring(loc+mDirName.length()+1);
	}*/

	/**
	 * 
	 * @param file
	 * @return
	 */
	public static Protein getProtein(File file) {
		new MolecularSystem() ;// Resets atom numbers
		AtomList atomList = new AtomList(file.getAbsolutePath());
		new MolecularSystem();
		Protein protein = new Protein(atomList, ResidueExtendedAtomsCreator.creator);
		return protein;
	}

	/**
	 * 
	 * @param protein
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	public static ConservationContactsRatio getConservationContactsRatio(Protein protein, String fileName) throws IOException{
		ConservationContactsCreator creator;
		CommandList commands = new CommandList("commands.txt",new CommandsException("Problem with commands"));
		ConSeq.setConSeq(protein.chain(), fileName); //Assign conservation to the residues
		creator = new ConservationContactsCreator(InfoType.CONTACTS11);
		creator.getWeight(commands);
		creator.setFileFound(true);
		ConservationContactsRatio score = ( ConservationContactsRatio)creator.createEnergyTerm(protein,null,null);
		return score;
	}


	public static ConservationContactsHrRatio getConservationContactsHrRatio(Protein protein, String fileName) throws IOException{
		ConservationContactsHrCreator creator;
		CommandList commands = new CommandList("commands.txt",new CommandsException("Problem with commands"));
		ConSeq.setConSeq(protein.chain(), fileName); //Assign conservation to the residues
		creator = new ConservationContactsHrCreator();
		creator.getWeight(commands);
		creator.setFileFound(true);
		ConservationContactsHrRatio score = (ConservationContactsHrRatio)creator.createEnergyTerm(protein,null,null);
		return score;
	}

	/**
	 * 
	 * @param file
	 */
	public static void getNative(File file) {
		String nativeFileName;
		String fileName = file.getAbsolutePath();
		if (fileName.endsWith("N.pdb")) nativeFileName = fileName;
		else {
			String[] subs = fileName.split("\\.");
			int index = fileName.indexOf(subs[subs.length-2]);
			nativeFileName = fileName.substring(0,index)+"N.pdb";
			System.out.println(nativeFileName);
		}
		if ((mNativeStructure == null)||
				nativeFileName.indexOf(mNativeStructure.name())== -1)
			mNativeStructure = getProtein(new File(nativeFileName));
	}

	/**
	 * 
	 * @param protein1
	 * @param protein2
	 * @return
	 */
	public static double getRms(Protein protein1, Protein protein2) {
        try {
		    return Rms.rms(protein1,protein2,ResidueAlignmentMethod.IDENTITY, Rms.RmsType.CA);// protein1,protein2,alignment);
	    }catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
    }

	/**
	 * 
	 * @param protein1
	 * @param protein2
	 * @return
	 */
	public static double[] getGdt(Protein protein1, Protein protein2)throws AlignmentException {
		return Rms.gdt(protein1,protein2);
	}

}
