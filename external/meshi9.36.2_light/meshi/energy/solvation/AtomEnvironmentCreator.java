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
public class AtomEnvironmentCreator extends EnergyCreator{
    private SolvationCreator solvationCreator;
    public AtomEnvironmentCreator(SolvationCreator solvationCreator) {
        super(InfoType.ATOM_ENVIRONMENT_ENERGY);
        this.solvationCreator = solvationCreator;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,CommandList commands) {
        AtomEnvironmentParameters atomEnvironmentParametersEnergy;
        AtomEnvironmentParameters atomEnvironmentParametersPropensity;
        SolvationEnergy solvationEnergy = (SolvationEnergy) solvationCreator.term();
        double[] cnc = solvationEnergy.cnc();
        double[] hbc = solvationEnergy.hbc();
        String parametersFileName = commands.firstWord("atomEnvironmentParameters").secondWord();
        File parametersFile = new File(parametersFileName);
        atomEnvironmentParametersEnergy     = new AtomEnvironmentParameters(parametersFile,"energy");
        atomEnvironmentParametersPropensity = new AtomEnvironmentParameters(parametersFile,"propensity");
        EnergyInfoElement info = new EnergyInfoElement(InfoType.ATOM_ENVIRONMENT_ENERGY,"atomEnvironmentEnergy");
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_BACKBONE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_POLAR));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_POSITIVE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_NEGATIVE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_AROMATIC));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEE_ALIPHATIC));
        info.getChildren().add(new EnergyInfoElement(InfoType.ATOM_ENVIRONMENT_PROPENSITY,"atomEnvironmentPropensity"));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_BACKBONE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_POLAR));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_POSITIVE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_NEGATIVE));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_AROMATIC));
        info.getChildren().add(new EnergyInfoElement(InfoType.AEP_ALIPHATIC));
        return new AtomEnvironmentEnergy(protein.atoms(),
                                         atomEnvironmentParametersEnergy,
                                         atomEnvironmentParametersPropensity,
                                         cnc,
                                         hbc,
                                         info);

    }
}

