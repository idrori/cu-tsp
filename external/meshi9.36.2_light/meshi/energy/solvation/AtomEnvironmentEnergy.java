package meshi.energy.solvation;

import meshi.energy.*;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueType;
import meshi.util.ResidueData;
import meshi.util.Utils;
import meshi.util.info.ChainsInfo;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ResidueInfo;

public class AtomEnvironmentEnergy extends AbstractEnergy implements EvaluatesResidues{
    private AtomList atoms;
    private double[] cnc;
    private double[] hbc;
    private AtomEnvironmentParameters energyParameters;
    private AtomEnvironmentParameters propensityParameters;

    public AtomEnvironmentEnergy() {
        super();
    }

    /**
     * See the comment at the top of the class for descriptions isOn the weights.
     */

    public AtomEnvironmentEnergy(AtomList atomList,
                                 AtomEnvironmentParameters energyParameters,
                                 AtomEnvironmentParameters propensityParameters,
                                 double[] cnc,
                                 double[] hbc,
                                 EnergyInfoElement info) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        this.atoms = atomList;
        this.energyParameters = energyParameters;
        this.propensityParameters = propensityParameters;
        this.cnc = cnc;
        this.hbc = hbc;
    }

    public boolean evaluatesResidues() {
        return true;
    }
    public EnergyInfoElement evaluate() {
        evaluateResidues(null);
        return info;
    }
    public void evaluateResidues(ChainsInfo chainsInfo) {
        double energy = 0;
        double propensity = 0;
        double cncValue, hbcValue, e,p;

        ResidueData residueEnergy = null;
        ResidueData residuePropensity = null;
        if (chainsInfo != null) {
            residueEnergy = new ResidueData(chainsInfo);
            residuePropensity = new ResidueData(chainsInfo);
        }

        for (MeshiInfo meshiInfo : info.getChildren())
            meshiInfo.setValue(new Double(0));
        for (int iAtom = 0; iAtom < atoms.size(); iAtom++) {
            AtomType atomType = atoms.atomAt(iAtom).type();
            if(atomType.isHydrogen() || (atomType == AtomType.TRN) || (atomType == AtomType.TRC) || (atomType == AtomType.TRO) ) continue;
            cncValue = cnc[iAtom];
            hbcValue = hbc[iAtom];
            energyParameters.getParameters(atomType).calc1(cncValue,hbcValue);
            e = energyParameters.getParameters(atomType).getValue();
            energy += e;
            propensityParameters.getParameters(atomType).calc1(cncValue,hbcValue);
            p = propensityParameters.getParameters(atomType).getValue();
            addToCategory(e,p,atomType);
            propensity += p;
            if (residueEnergy != null) {
                residueEnergy.add(atoms.atomAt(iAtom), e);
                residuePropensity.add(atoms.atomAt(iAtom), p);
            }
        }
        if (residueEnergy != null) {
            chainsInfo.add(residueEnergy, InfoType.ATOM_ENVIRONMENT_ENERGY);
            chainsInfo.add(residuePropensity, InfoType.ATOM_ENVIRONMENT_PROPENSITY);
        }
        info.setValue(energy);
        for (MeshiInfo meshiInfo : info.getChildren()) {
            if (meshiInfo.type == InfoType.ATOM_ENVIRONMENT_PROPENSITY)
                meshiInfo.setValue(propensity);
        }
   }

    private void addToCategory(double energy, double propensity, AtomType type) {
        for (MeshiInfo meshiInfo : info.getChildren()) {
            double currentValue = ((Double) meshiInfo.getValue()).doubleValue();
            if (type.backbone()) {
                if (meshiInfo.type == InfoType.AEE_BACKBONE) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_BACKBONE) meshiInfo.setValue(new Double(currentValue + propensity));
            }
            else if (ResidueType.isPolar(type)) {
                if (meshiInfo.type == InfoType.AEE_POLAR) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_POLAR) meshiInfo.setValue(new Double(currentValue + propensity));
            }
            else if (ResidueType.isPositive(type)){
                if (meshiInfo.type == InfoType.AEE_POSITIVE) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_POSITIVE) meshiInfo.setValue(new Double(currentValue + propensity));
            }
            else if (ResidueType.isNegative(type)){
                if (meshiInfo.type == InfoType.AEE_NEGATIVE) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_NEGATIVE) meshiInfo.setValue(new Double(currentValue + propensity));
            }
            else if (ResidueType.isAromatic(type)){
                if (meshiInfo.type == InfoType.AEE_AROMATIC) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_AROMATIC) meshiInfo.setValue(new Double(currentValue + propensity));
            }
            else if (ResidueType.isAliphatic(type)){
                if (meshiInfo.type == InfoType.AEE_ALIPHATIC) meshiInfo.setValue(new Double(currentValue + energy));
                if (meshiInfo.type == InfoType.AEP_ALIPHATIC) meshiInfo.setValue(new Double(currentValue + propensity));
            }
        }
    }
    public void evaluateAtoms() {}
    public void test(TotalEnergy energy, Atom atom) {}
}
