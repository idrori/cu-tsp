/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvateRot1;

import meshi.energy.*;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

public class SolvateRot1Creator extends EnergyCreator implements KeyWords {

    private double[][] pp = null;
    private double maximalZ = 1e10;
    private static SolvateRot1ParametersList parametersList = null;

/*
    public SolvateRot1Creator(double weight,double maximalZ) {
        super(weight);
        this.maximalZ = maximalZ;
    }

    public SolvateRot1Creator(double maximalZ) {
        super(SOLVATE_ENERGY);
        this.maximalZ = maximalZ;
    }
*/

    public SolvateRot1Creator(double weight, double[][] pp, double maximalZ) {
        super(weight,InfoType.SOLVATION);
        this.pp = pp;
        this.maximalZ = maximalZ;
    }

    public SolvateRot1Creator(double[][] pp, double maximalZ) {
        super(InfoType.SOLVATION);
        this.pp = pp;
        this.maximalZ = maximalZ;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
/*
if (pp==null) {
        DunbrackLib lib = new DunbrackLib(commands, 1.0 , 2);
        pp = RotamericTools.putIntoRot1(protein, distanceMatrix, lib);
}
*/
        if (parametersList == null) {
            // Appanding the path to the list of parameter filenames.
            String[] strlist = new String[SOLVATE_NOHB_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < SOLVATE_NOHB_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(SOLVATE_NOHB_PARAMETERS[cc]);
            parametersList = new SolvateRot1ParametersList(strlist);
        }
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.SOLVATION, "Rot1 Solvation", weight);
        SolvateRot1Energy solvateEnergy = new SolvateRot1Energy(protein.atoms(),
                distanceMatrix,
                (SolvateRot1ParametersList) parametersList,
                pp,
                maximalZ,
                energyInfoElement);
        solvateEnergy.setComment("Rot1 Solvate");
        return solvateEnergy;
    }

}
