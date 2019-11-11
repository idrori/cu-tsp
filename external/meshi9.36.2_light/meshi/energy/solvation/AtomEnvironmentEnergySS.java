package meshi.energy.solvation;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.AtomType;
import meshi.parameters.SecondaryStructure;
import meshi.util.Utils;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

public class AtomEnvironmentEnergySS extends AbstractEnergy {
    private AtomList atoms;
    private double[] cnc;
    private double[] hbc;
    AtomEnvironmentParameters atomEnvironmentParametersEnergyHelix;
    AtomEnvironmentParameters atomEnvironmentParametersPropensityHelix;
    AtomEnvironmentParameters atomEnvironmentParametersEnergySheet;
    AtomEnvironmentParameters atomEnvironmentParametersPropensitySheet;
    AtomEnvironmentParameters atomEnvironmentParametersEnergyCoil;
    AtomEnvironmentParameters atomEnvironmentParametersPropensityCoil;

    public AtomEnvironmentEnergySS() {
        super();
    }

    /**
     * See the comment at the top of the class for descriptions isOn the weights.
     */

    public AtomEnvironmentEnergySS(AtomList atomList,
                                   AtomEnvironmentParameters atomEnvironmentParametersEnergyHelix,
                                   AtomEnvironmentParameters atomEnvironmentParametersPropensityHelix,
                                   AtomEnvironmentParameters atomEnvironmentParametersEnergySheet,
                                   AtomEnvironmentParameters atomEnvironmentParametersPropensitySheet,
                                   AtomEnvironmentParameters atomEnvironmentParametersEnergyCoil,
                                   AtomEnvironmentParameters atomEnvironmentParametersPropensityCoil,
                                   double[] cnc,
                                   double[] hbc,
                                   EnergyInfoElement info) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        this.atoms = atomList;
        this.atomEnvironmentParametersEnergyHelix     = atomEnvironmentParametersEnergyHelix;
        this.atomEnvironmentParametersPropensityHelix = atomEnvironmentParametersPropensityHelix;
        this.atomEnvironmentParametersEnergySheet     = atomEnvironmentParametersEnergySheet;
        this.atomEnvironmentParametersPropensitySheet = atomEnvironmentParametersPropensitySheet;
        this.atomEnvironmentParametersEnergyCoil      = atomEnvironmentParametersEnergyCoil;
        this.atomEnvironmentParametersPropensityCoil  = atomEnvironmentParametersPropensityCoil;
        this.cnc = cnc;
        this.hbc = hbc;
    }

    public EnergyInfoElement evaluate() {
        double energy = 0;
        double propensity = 0;
        double cncValue, hbcValue, e, p, currentE, currentP;
        for (MeshiInfo meshiInfo : info.getChildren())
            meshiInfo.setValue(new Double(0));

        for (int iAtom = 0; iAtom < atoms.size(); iAtom++) {
            Atom atom = atoms.atomAt(iAtom);
            AtomType atomType = atom.type();
            if (atomType.isHydrogen() || (atomType == AtomType.TRN) || (atomType == AtomType.TRC) || (atomType == AtomType.TRO))
                continue;
            cncValue = cnc[iAtom];
            hbcValue = hbc[iAtom];
            SecondaryStructure secondaryStructure = atom.residue().getSecondaryStructure();
            AtomEnvironmentParameters atomEnvironmentEnergyParameters, atomEnvironmentPropensityParameters;
            MeshiInfo energyInfo, propensityInfo;
            if (secondaryStructure == SecondaryStructure.HELIX) {
                atomEnvironmentEnergyParameters = atomEnvironmentParametersEnergyHelix;
                atomEnvironmentPropensityParameters = atomEnvironmentParametersPropensityHelix;
                energyInfo = getInfo(info, InfoType.AEE_HELIX);
                propensityInfo = getInfo(info, InfoType.AEP_HELIX);
            } else if (secondaryStructure == SecondaryStructure.SHEET) {
                atomEnvironmentEnergyParameters = atomEnvironmentParametersEnergySheet;
                atomEnvironmentPropensityParameters = atomEnvironmentParametersPropensitySheet;
                energyInfo = getInfo(info, InfoType.AEE_SHEET);
                propensityInfo = getInfo(info, InfoType.AEP_SHEET);
            } else if (secondaryStructure == SecondaryStructure.COIL) {
                atomEnvironmentEnergyParameters = atomEnvironmentParametersEnergyCoil;
                atomEnvironmentPropensityParameters = atomEnvironmentParametersPropensityCoil;
                energyInfo = getInfo(info, InfoType.AEE_COIL);
                propensityInfo = getInfo(info, InfoType.AEP_COIL);
            } else throw new RuntimeException("This is weird "+secondaryStructure+" "+atom+" "+atom.residue());
            atomEnvironmentEnergyParameters.getParameters(atomType).calc1(cncValue, hbcValue);
            atomEnvironmentPropensityParameters.getParameters(atomType).calc1(cncValue, hbcValue);
            e = atomEnvironmentEnergyParameters.getParameters(atomType).getValue();
            p = atomEnvironmentPropensityParameters.getParameters(atomType).getValue();
            energy += e;
            propensity += p;
            currentE = ((Double) energyInfo.getValue()).doubleValue();
            currentP = ((Double) propensityInfo.getValue()).doubleValue();
            energyInfo.setValue(new Double(currentE + e));
            propensityInfo.setValue(new Double(currentP + p));
        }
        info.setValue(energy);
        for (MeshiInfo meshiInfo : info.getChildren()) {
            if (meshiInfo.type == InfoType.ATOM_ENVIRONMENT_PROPENSITY_SS)
                meshiInfo.setValue(propensity);
        }
        return info;
    }

    private static MeshiInfo getInfo(MeshiInfo meshiInfo, InfoType infoType) {
        for (MeshiInfo info : meshiInfo.getChildren()) {
            if (info.type == infoType) return info;
        }
        return null;
    }

    public void evaluateAtoms() {}
    public void test(TotalEnergy energy, Atom atom) {}
}
