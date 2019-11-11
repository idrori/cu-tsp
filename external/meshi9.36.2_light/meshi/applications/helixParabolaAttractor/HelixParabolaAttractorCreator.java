/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.helixParabolaAttractor;

import meshi.energy.*;
import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.util.CommandList;
import meshi.util.filters.Filter;
import meshi.parameters.SecondaryStructure;
import meshi.util.info.InfoType;

/**
 * Creates a parabole around the optimal torsion angles for helices. Each
 * residue in helix secondary structure is encouraged to converge to phi=-40,
 * psi=-60, using the following formula:
 * <p/>
 * e1phi = (phi- phiTarget).^4./ ((phi- phiTarget).^4+a)-1;
 * e2phi = (phi- max).^4./( (phi- max).^4+epsilon );
 * e3phi = (phi-min).^4./((phi-min).^4+epsilon);
 * ephi = e2phi.*e1phi.*e3phi;
 * <p/>
 * the same for psi
 *
 * @author Ami Levy-Moonshine
 *         31/8/09
 */
public class HelixParabolaAttractorCreator extends EnergyCreator {

    public HelixParabolaAttractorCreator(double weight) {
        super(weight,InfoType.HELIX_ATTRACTOR);
    }

    public HelixParabolaAttractorCreator() {
        super(InfoType.HELIX_ATTRACTOR);
    }

    public AbstractEnergy createEnergyTerm(
            Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {

        /* create a list of all (phi,psi) torsion pairs */
        TorsionPairList torsionPairList = TorsionPairList.createQuickAndDirtyTorsionPairList(protein, distanceMatrix, new PhiPsiFilter());

        /* isolate helix only */
        TorsionPairList helixPairList = torsionPairList.filter(new HelixOnlyFilter());
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.HELIX_ATTRACTOR, "Forces regions with declared helix secondary structure to actually become helicall", weight);
        return term = new HelixParabolaAttractorEnergy(helixPairList, distanceMatrix, energyInfoElement);
    }

    /**
     * Filters phi,psi torsions only.
     */
    private static class PhiPsiFilter implements Filter {
        public boolean accept(Object obj) {
            return (((Torsion) obj).getTorsionName().equals("PHI") ||
                    ((Torsion) obj).getTorsionName().equals("PSI"));
        }
    }

    /**
     * filters torsions of residues marked as helix
     */
    private static class HelixOnlyFilter implements Filter {
        public boolean accept(Object obj) {
            /* accept torsion pair objects only */
            if (obj instanceof TorsionPair) {
                TorsionPair pair = (TorsionPair) obj;
                /* and verify they fit helix secondary structure */
                return pair.torsion1().atom2.residue().getSecondaryStructure().equals(SecondaryStructure.HELIX);
            }
            return false;
        }
    }

}
