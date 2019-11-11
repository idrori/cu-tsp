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
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsEnergyElement;
import meshi.energy.hydrogenBondsPairs.PairOfHydrogenBondsElements;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.MeshiPotential;
import meshi.sequences.AlignmentException;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;
import meshi.util.Utils;

/**
 * Created by chen on 11/12/2014.
 */
public class HBondAnalyzer {
    public static void main(String[] args) throws UpdateableException,EvaluationException, AlignmentException{
        Protein model = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        System.out.println("new model "+ model);
        Utils.AssignDSSP(model, args[1]);
        model.chain().print();
        CommandList commands = new CommandList(args[2]);
        HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
        HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
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
        System.out.println(hydrogenBondsEnergy);
        System.out.println("\n\n H-Bonds list\n");
        for (Distance distance : hydrogenBondsEnergy.hBondList())
                if (accept(distance)) display(distance);


        System.out.println("\n\n H-Bond pairs list\n");
        for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()){
            Distance hb1 = element.HOelement1();
            Distance hb2 = element.HOelement2();
            if (accept(hb1) & accept(hb2)) {
                System.out.println(element);
                display(hb1);
                display(hb2);
                System.out.println();
            }
        }
      }
     private static boolean accept(Distance distance) {
         int[] range = {40,69};
         Atom atom1 = distance.atom1();
         Atom atom2 = distance.atom2();
         int n1 = atom1.residueNumber();
         int n2 = atom2.residueNumber();
         if ((n1 >= range[0]) && (n2 >= range[0]) && (n1 <= range[1]) && (n2 <= range[1])) return true;
         return false;
     }
    private static void display(Distance distance) {
        Atom atom1 = distance.atom1();
        Atom atom2 = distance.atom2();
        System.out.println(distance+" "+atom1.residue()+" "+atom1.type()+" ; "+atom2.residue()+" "+atom2.type());
    }
}
