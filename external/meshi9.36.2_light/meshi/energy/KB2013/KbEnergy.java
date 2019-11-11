package meshi.energy.KB2013;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.scoringFunctions.Score;
import meshi.util.UpdateableException;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.*;

import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 20/07/13
 * Time: 21:36
 * To change this template use File | Settings | File Templates.
 */
public class KbEnergy implements Score{
    private KbPotential potential;
    private Protein protein;
    private DistanceMatrix distanceMatrix;
    public KbEnergy(Protein protein, KbPotential kbPotential){
        MolecularSystem molecularSystem = protein.atoms().molecularSystem();
        molecularSystem.terminator().reset();
        distanceMatrix = new DistanceMatrix(molecularSystem,
                                            DistanceMatrix.DistanceMatrixType.LONG_DISTANCES);
        this.protein = protein;
        this.potential = kbPotential;
    }

    public MeshiInfo score(MeshiInfo energyInfo){
        AtomPairIterator atomPairIterator = new AtomPairIterator(protein.atoms(),4);

        double score = 0;
        while(atomPairIterator.hasNext()) {
            AtomPair pair = atomPairIterator.next();
            Atom atom1 = pair.atom1();
            Atom atom2 = pair.atom2();
            double distance = atom1.distanceFrom(atom2);
            score += potential.getScore(atom1,atom2,distance,protein);
        }

        MeshiInfo out = new MeshiInfo(InfoType.UNKNOWN, score, "Place holder");
        return out;
    }

    public static void test(String dataFilesListName, KbPotential kbPotential) throws IOException,UpdateableException {
        String          dataLine;
        File dataFilesList       = new File(dataFilesListName);
        MeshiLineReader dataFilesListReader = new MeshiLineReader(dataFilesList);
        while ((dataLine = dataFilesListReader.readLine())!= null){
            String pdbFileName = dataLine.substring(0,4).toUpperCase();
            String chainName   = dataLine.substring(4, 5);
            System.out.printf("%s %s \n", pdbFileName, chainName);
            AtomList atoms   = new AtomList(pdbFileName);
            Protein  protein = new Protein(atoms, new ResidueExtendedAtomsCreator());
            System.out.println(protein);
            KbEnergy energy = new KbEnergy(protein,kbPotential);
            System.out.println(energy.score(null));
        }
    }
}
