package meshi.util.proteinReaders;

import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.Loop;
import meshi.parameters.SecondaryStructure;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.Utils;
import meshi.util.dssp.DSSP;
import meshi.util.dssp.DsspReader;

import java.io.File;

/**
 * Created by chen on 02/12/2015.
 */
public class StandardProteinReader implements ProteinReader{
    DsspReader dsspReader;
    CommandList commands;

    public StandardProteinReader(DsspReader dsspReader, String commandsFileName) {
        this.dsspReader = dsspReader;
        commands =  new CommandList(commandsFileName);
    }

    public Protein readProtein(ProteinMetaData proteinMetaData) {
        Command command = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        BondParametersList bondParametersList  = new BondParametersList(parametersDirectory+"/" + Loop.BOND_PARAMETERS);
        AngleParametersList angleParametersList = new AngleParametersList(parametersDirectory+"/" + Loop.ANGLE_PARAMETERS);
        String directory = (String) proteinMetaData.get(ProteinMetaData.MetaDataKey.DIRECTORY);
        if (directory == null) throw new RuntimeException("This is weird.");
        String fileName  = (String) proteinMetaData.get(ProteinMetaData.MetaDataKey.FILE_NAME);
        if (fileName == null) throw new RuntimeException("This is weird.");
        File file = new File(directory,fileName);
        AtomList atomList = new AtomList(file.getAbsolutePath());
        Protein protein = new Protein(atomList, ResidueExtendedAtomsCreator.creator);
        protein.molecularSystem.terminator().reset();
        for (ProteinMetaData.MetaDataKey key : ProteinMetaData.MetaDataKey.values()) {
            Object obj = proteinMetaData.get(key);
            if (obj != null)
                protein.putMetaData(key,obj);
        }

        for (Residue residue: protein.residues()) {
            if (!residue.dummy()) ((ResidueExtendedAtoms) residue).addHydrogens(bondParametersList, angleParametersList);
        }
       if (dsspReader != null) {
           DSSP dssp = dsspReader.read(protein);
           if (!Utils.AssignDSSP(protein, dssp))    {
               System.out.println("DSSP assignment failed. Protein "+protein+" is ignored.");
           }

           for (int iResidue = 0; iResidue< protein.residues().size(); iResidue++) {
               Residue residue = protein.residues().get(iResidue);
               if (residue.getSecondaryStructure() == SecondaryStructure.UNK) {
                   residue.setSecondaryStructure(SecondaryStructure.COIL);
               }
           }
        }
        else throw new RuntimeException("Cannot find DSSP to "+protein);
        return protein;
    }
}
