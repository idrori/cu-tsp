package meshi.energy.solvation;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

import java.io.File;

/**
 * Created by chen on 03/04/2016.
 */
public class AtomEnvironmentCreatorSS extends EnergyCreator{
    private SolvationCreator solvationCreator;
    public AtomEnvironmentCreatorSS(SolvationCreator solvationCreator) {
        super(InfoType.ATOM_ENVIRONMENT_ENERGY_SS);
        this.solvationCreator = solvationCreator;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,CommandList commands) {
        AtomEnvironmentParameters atomEnvironmentParametersEnergyHelix;
        AtomEnvironmentParameters atomEnvironmentParametersPropensityHelix;
        AtomEnvironmentParameters atomEnvironmentParametersEnergySheet;
        AtomEnvironmentParameters atomEnvironmentParametersPropensitySheet;
        AtomEnvironmentParameters atomEnvironmentParametersEnergyCoil;
        AtomEnvironmentParameters atomEnvironmentParametersPropensityCoil;
        SolvationEnergy solvationEnergy = (SolvationEnergy) solvationCreator.term();
        double[] cnc = solvationEnergy.cnc();
        double[] hbc = solvationEnergy.hbc();
        String helixParametersFileName = commands.firstWord("atomEnvironmentParametersHelix").secondWord();
        String sheetParametersFileName = commands.firstWord("atomEnvironmentParametersSheet").secondWord();
        String coilParametersFileName  = commands.firstWord("atomEnvironmentParametersCoil").secondWord();
        File helixParametersFile = new File(helixParametersFileName);
        File sheetParametersFile = new File(sheetParametersFileName);
        File coilParametersFile  = new File(coilParametersFileName);
        atomEnvironmentParametersEnergyHelix     = new AtomEnvironmentParameters(helixParametersFile,"energy");
        atomEnvironmentParametersPropensityHelix = new AtomEnvironmentParameters(helixParametersFile,"propensity");
        atomEnvironmentParametersEnergySheet     = new AtomEnvironmentParameters(sheetParametersFile,"energy");
        atomEnvironmentParametersPropensitySheet = new AtomEnvironmentParameters(sheetParametersFile,"propensity");
        atomEnvironmentParametersEnergyCoil      = new AtomEnvironmentParameters(coilParametersFile,"energy");
        atomEnvironmentParametersPropensityCoil  = new AtomEnvironmentParameters(coilParametersFile,"propensity");
        EnergyInfoElement info = new EnergyInfoElement(InfoType.ATOM_ENVIRONMENT_ENERGY_SS, InfoType.ATOM_ENVIRONMENT_ENERGY_SS.tag );
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_COIL));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_HELIX));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_SHEET));
        info.getChildren().add(new EnergyInfoElement(InfoType.ATOM_ENVIRONMENT_PROPENSITY_SS));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_COIL));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_HELIX));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_SHEET));
        return new AtomEnvironmentEnergySS(protein.atoms(),
                                           atomEnvironmentParametersEnergyHelix,
                                           atomEnvironmentParametersPropensityHelix,
                                           atomEnvironmentParametersEnergySheet,
                                           atomEnvironmentParametersPropensitySheet,
                                           atomEnvironmentParametersEnergyCoil,
                                           atomEnvironmentParametersPropensityCoil,
                                           cnc,
                                           hbc,
                                           info);
    }
}

