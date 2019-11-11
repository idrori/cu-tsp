package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsEnergy;
import meshi.energy.hydrogenBondsPairs.PairOfHydrogenBondsElements;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentException;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.Iterator;

/**
 * Created by chen on 11/12/2014.
 */
public class STR2Analyzer {
    static MeshiWriter output;
    static MeshiWriter result;
    static HydrogenBondsCreator hydrogenBondsCreator;
    static HydrogenBondsPairsCreator hydrogenBondsPairsCreator;

    public static void main(String[] args) throws UpdateableException,EvaluationException, AlignmentException, IOException{
        // 1 - initialize writers
        output = new MeshiWriter("WriterTest.siditom.txt");
        result = new MeshiWriter(args[0]+".str2");
        // 2 - create protein model
        Protein model = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        result.print(model.chain().sequence());                  // print chain sequence to result file
        // 3 - assign dssp to model
        Utils.AssignDSSP(model, args[1]);
        result.print(getChainSecondaryStructure(model.chain())); // print secondary structure to result file
        // 4 - calculate hydrogen bonds: add hydrogen atoms via minimization // 5 - calculate hydrogen bond pairs

        /** Assumes PairOfHydrogenBondsElements.MAX_SEQ_PAIRS_DISTANCE = 5  for pairs in the same beta orientation **/
        CommandList commands = new CommandList(args[2]);
        STR2Analyzer.calculateHydrogenBonds(model, commands);

        // 6 - assign str2 alphabet on top dssp
        STR2Analyzer.assignSTR2(model,hydrogenBondsPairsCreator);

        result.close();
        output.close();
      }
     private static boolean accept(Distance distance) {
         int[] range = {3,504};
         Atom atom1 = distance.atom1();
         Atom atom2 = distance.atom2();
         int n1 = atom1.residueNumber();
         int n2 = atom2.residueNumber();
         if ((n1 >= range[0]) && (n2 >= range[0]) && (n1 <= range[1]) && (n2 <= range[1]) && (distance.distance() < 2.5) ) return true;
         return false;
     }
    private static void display(Distance distance) {
        Atom atom1 = distance.atom1();
        Atom atom2 = distance.atom2();
        output.println(distance+" "+atom1.residue()+" "+atom1.type()+" ; "+atom2.residue()+" "+atom2.type());
    }

    public static String getChainSecondaryStructure(Chain c){
        String sS = "";
        Iterator residues = c.iterator();
        Residue r = (Residue) residues.next();
        String s = "";
        while (residues.hasNext()){
            r = (Residue) residues.next();
            sS+=r.getSecondaryStructure().getNameOneLetter();
        }
        sS+="\n";
        return sS;
    }

    public static boolean calculateHydrogenBonds(Protein model, CommandList commands) throws UpdateableException,EvaluationException, AlignmentException, IOException{
        hydrogenBondsCreator = new HydrogenBondsCreator();
        hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
        AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();

        EnergyCreator[] energyCreators = {
                summaCreator,
                hydrogenBondsCreator,
                hydrogenBondsPairsCreator,
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
        };
        model.atoms().molecularSystem.createDistanceMatrix(DistanceMatrix.DistanceMatrixType.STANDARD.toString(),DistanceMatrix.DistanceMatrixType.STANDARD);
        model.atoms().molecularSystem.createDistanceMatrix(DistanceMatrix.DistanceMatrixType.STANDARD.toString(),DistanceMatrix.DistanceMatrixType.STANDARD);
        MeshiProgram.initRandom(0);
        System.out.println(" --------------------- Add getAtoms starts ----------------------------------------") ;
        Utils.addAtoms(model, false, commands,
                new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
                        "/" + MeshiPotential.PLANE_PARAMETERS), false);
        System.out.println(" --------------------- Add getAtoms done ----------------------------------------") ;

        //TotalEnergy minimizationEnergy = new TotalEnergy(model,model.atoms().molecularSystem.getDistanceMatrix(),energyCreators,commands,model.atoms(),"My energy function.");
        TotalEnergy minimizationEnergy = new TotalEnergy(model,energyCreators,commands,model.atoms(),"My energy function.");
        minimizationEnergy.evaluate();

        HydrogenBondsEnergy hydrogenBondsEnergy = (HydrogenBondsEnergy) hydrogenBondsCreator.term();
        output.println(hydrogenBondsEnergy);
        output.println("\n\n H-Bonds list\n");
        for (Distance distance : hydrogenBondsEnergy.hBondList())
            if (accept(distance)) display(distance);

        output.println("\n\n H-Bond pairs list\n");
        for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()){
            Distance hb1 = element.HOelement1();
            Distance hb2 = element.HOelement2();
            if (accept(hb1) & accept(hb2)) {
                output.println(element);
                display(hb1);
                display(hb2);
                output.println();
            }
        }
        return true;
    }
    public static boolean isResidueConnected(Residue r, HydrogenBondsPairsCreator hydrogenBondsPairsCreator){
        for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
            Distance hb1 = element.HOelement1();
            Distance hb2 = element.HOelement2();
            if ((element.isBetaPair() & accept(hb1) & accept(hb2)) && (hb1.atom1().residue().number() == r.number() || hb1.atom2().residue().number() == r.number() ||
                    hb2.atom1().residue().number() == r.number() || hb2.atom2().residue().number() == r.number())){ //its a beta pair, with proper distances for a hydrogen bond, and one of them is with the residue r.
                return true;
            }
        }
        return false;
    }

    //Assumption: the Ss orientation of r is beta sheet.
    public static PairOfHydrogenBondsElements.SsElement getResidueBetaSs(Residue r, HydrogenBondsPairsCreator hydrogenBondsPairsCreator){
        PairOfHydrogenBondsElements.SsElement  residueSs = PairOfHydrogenBondsElements.SsElement.BETA;
        boolean changed = false;
        for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
            Distance hb1 = element.HOelement1();
            Distance hb2 = element.HOelement2();

            if ((element.isBetaPair() & accept(hb1) & accept(hb2)) && (hb1.atom1().residue().number() == r.number() || hb1.atom2().residue().number() == r.number() ||
                    hb2.atom1().residue().number() == r.number() || hb2.atom2().residue().number() == r.number()) ) { //its a beta pair, with proper distances for a hydrogen bond, and one of them is with the residue r.
                if (!changed && element.getSsElement() == PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA) {
                    changed = true;
                    residueSs = PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA;
                } else if (!changed && element.getSsElement() == PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA) {
                    changed = true;
                    residueSs = PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA;
                } else if (changed && (element.getSsElement() != PairOfHydrogenBondsElements.SsElement.BETA) && (residueSs != element.getSsElement())) {
                    throw new RuntimeException("The residue secondary structure orientation exists in both: ANTI_PARALLEL_BETA & PARALLEL_BETA");
                }//if
            }//if
        }//for
        return residueSs;
    }


    public static char assignResidueSTR2(Protein model, Residue r, HydrogenBondsPairsCreator hydrogenBondsPairsCreator){
        boolean isOneSided = false;
        boolean isTwoSided = false;
        int residueNumber = r.ID().number(); //getResidueNumber
        Residue prevResidue = model.residue(residueNumber-1); //get previous residue in chain
        Residue nextResidue = model.residue(residueNumber+1); //get next residue in chain
        //if residue (i-1) or (i+1) is connected -> two-sided, else one sided.
        boolean rConnected =STR2Analyzer.isResidueConnected(r,hydrogenBondsPairsCreator);
        boolean prevConnected = STR2Analyzer.isResidueConnected(prevResidue,hydrogenBondsPairsCreator);
        boolean nextConnected = STR2Analyzer.isResidueConnected(nextResidue,hydrogenBondsPairsCreator);
        if ((rConnected) && ( prevConnected || nextConnected)) {
            isTwoSided = true;
        } else if ((rConnected && !prevConnected && !nextConnected) ||
                (!rConnected && (prevConnected || nextConnected))){
            isOneSided = true;
        }
        try {
            PairOfHydrogenBondsElements.SsElement residueSs = STR2Analyzer.getResidueBetaSs(r, hydrogenBondsPairsCreator);
            PairOfHydrogenBondsElements.SsElement prevResidueSs = STR2Analyzer.getResidueBetaSs(prevResidue, hydrogenBondsPairsCreator);
            PairOfHydrogenBondsElements.SsElement nextResidueSs = STR2Analyzer.getResidueBetaSs(nextResidue, hydrogenBondsPairsCreator);
            if (isOneSided){//One sided Sheet
                    if (residueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA) || prevResidueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA) || nextResidueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA))
                        return 'Z';//One sided Anti-Parallel (Z)
                    else if (residueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA) || prevResidueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA) || nextResidueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA))
                        return 'Q';//One sided Parallel (Q)
                    else return '_';
            }else if (isTwoSided){ //Two sided Sheet
                    if (residueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA) && (prevResidueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA) || nextResidueSs.equals(PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA))){
                        return 'A';//two sided Anti-Parallel (A)
                    } else if (residueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA) && (prevResidueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA) || nextResidueSs.equals(PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA))){
                        return 'P';//two sided Parallel (P)
                    } else return 'M';//two sided Mixed (M)
            } else return 'E';
        } catch(Exception e){
            System.out.println(e.getMessage());
            return '-';
        }//try
    }

    public static boolean assignSTR2(Protein model, HydrogenBondsPairsCreator hydrogenBondsPairsCreator) {
        int sizeOfChain = model.chain().size();
        for (int iResidue=1; iResidue < sizeOfChain; iResidue++) {
            Residue r = model.chain().residueAt(iResidue);
            SecondaryStructure sSResidue = r.getSecondaryStructure();
            if (!sSResidue.equals(SecondaryStructure.SHEET)){
                result.print(sSResidue.getNameOneLetter());
            }else result.print(STR2Analyzer.assignResidueSTR2(model, r,hydrogenBondsPairsCreator));

        }
        return true;
    }

}
