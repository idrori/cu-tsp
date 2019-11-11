package meshi.energy.simpleEnergyTerms.torsionConstraints;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.SimpleEnergyCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.geometry.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.atoms.AtomPairList;
import meshi.util.CommandList;
import meshi.util.info.InfoType;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

public class TorsionConstraintsCreator extends SimpleEnergyCreator{
    JSONObject torsions;
    public TorsionConstraintsCreator(JSONObject torsions) {
        super(InfoType.TORSION_CONSTRAINT);
        this.torsions = torsions;
    }
    public TorsionConstraintsCreator(double[][][] torsions) throws ParseException{
        super(InfoType.TORSION_CONSTRAINT);
        this.torsions = new JSONObject();
        for (int i = 0; i < torsions.length; i++) {
            if ((torsions[i][0][0] > -1000) & (torsions[i][1][0] > -1000)) {
                double[][] phiPsiOmega = {{Math.PI * torsions[i][0][0] / 180, torsions[i][0][1]}, {Math.PI * torsions[i][1][0] / 180, torsions[i][1][1]}};
                this.torsions.put("" + i, phiPsiOmega);
            }
        }
        String torsionsString = this.torsions.toJSONString();
        JSONParser parser = new JSONParser();
        this.torsions =  (JSONObject) parser.parse(torsionsString);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands) {

        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.TORSION_CONSTRAINT,
                "Torsion constraints energy", weight);
        AtomPairList allBonds = new AtomPairList(protein.chains(),null);
        AtomPairList bonds = new AtomPairList();
        for (AtomPair bond : allBonds) {
            if ((!bond.atom1().nowhere()) & (!bond.atom2().nowhere()))
                bonds.add(bond);
        }
        QuickAndDirtyAngleList angleList = new QuickAndDirtyAngleList(bonds,dm);
        QuickAndDirtyTorsionList allTorsions = new QuickAndDirtyTorsionList(angleList,dm);

        TorsionList torsionList = new TorsionList();
        TorsionConstraintsParametersList targets = new TorsionConstraintsParametersList();
        Chain chain = protein.chain();
        for (Residue residue : chain) {
            if ((!residue.dummy()) && (!residue.ca().nowhere())) {
                JSONArray residueTorsionsArray = (JSONArray) torsions.get("" + residue.number());
                if (residueTorsionsArray != null) {
                    ResidueTorsions residueTorsions = new ResidueTorsions(residue, allTorsions);
                    Torsion phiTorsion = residueTorsions.phi();
                    Torsion psiTorsion = residueTorsions.psi();
                    Torsion omegaTorsion = residueTorsions.omega();
                    double[] phi = new double[2], psi = new double[2], omega = new double[2];

                    if (ok(phiTorsion) ) {
                        JSONArray phiArray = (JSONArray) residueTorsionsArray.get(0);
                        phi[0] = (double) phiArray.get(0);
                        phi[1] = (double) phiArray.get(1);
                        torsionList.add(phiTorsion);
                        targets.add(new TorsionConstraintsParameters(phi));

                    }
                    if (ok(psiTorsion)) {
                        JSONArray psiArray = (JSONArray) residueTorsionsArray.get(1);
                        psi[0] = (double) psiArray.get(0);
                        psi[1] = (double) psiArray.get(1);
                        torsionList.add(psiTorsion);
                        targets.add(new TorsionConstraintsParameters(psi));
                    }
                    if (ok(omegaTorsion)) {
                        if (residueTorsionsArray.size() > 2) {
                            JSONArray omegaArray = (JSONArray) residueTorsionsArray.get(2);
                            omega[0] = (double) omegaArray.get(0);
                            omega[1] = (double) omegaArray.get(1);
                            torsionList.add(omegaTorsion);
                            targets.add(new TorsionConstraintsParameters(omega));
                        }
                    }
                } //else Utils.println("No torsions for " + residue);
            }
        }
        term = new TorsionConstraints(torsionList, targets, energyInfoElement);
        return term;
    }

    private boolean ok(Torsion torsion) {
        if (torsion == null) return false;
        if (torsion.atom1.nowhere()) return false;
        if (torsion.atom2.nowhere()) return false;
        if (torsion.atom3.nowhere()) return false;
        if (torsion.atom4.nowhere()) return false;
        return true;
    }

    /*
                phiC1 = phiC2 = phiCA = phiN = null;
            psiN1 = psiCA = psiC  = psiN2 = omegaCA1 = omegaC   = omegaN = omegaCA2 = null;
            if (!residue.dummy()) {
                int residueNumber = residue.number();
                if (residueNumber > 0)
                    prevResidue = chain.residueAt(residueNumber - 1);
                else prevResidue = null;
                if (residueNumber < chain.size()-1)
                    nextResidue = chain.residueAt(residueNumber + 1);
                else nextResidue = null;
                if ((prevResidue != null) && (!prevResidue.dummy())) {
                    phiC1 = prevResidue.carbonylC();
                    phiN  = residue.amideN();
                    phiCA = residue.ca();
                    phiC2 = residue.carbonylC();
                }
                if ((nextResidue != null) && (!nextResidue.dummy())) {
                    psiN1 = residue.amideN();
                    psiCA = residue.ca();
                    psiC  = residue.carbonylC();
                    psiN2 = nextResidue.amideN();
                    omegaCA1 = residue.ca();
                    omegaC   = residue.carbonylC();
                    omegaN   = nextResidue.amideN();
                    omegaCA2 = nextResidue.ca();
                }
     */
}
