package meshi.energy.cooperativeDistanceMapConstraints;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.SimpleEnergyCreator;
import meshi.geometry.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomBuilder;
import meshi.parameters.ResidueType;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.info.InfoType;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.util.ArrayList;


public class CooperativeDistanceMapConstraintsCreator extends SimpleEnergyCreator{
    JSONObject distances;
    String mode;
    Chain chain;
    public CooperativeDistanceMapConstraintsCreator(JSONObject distances, String mode) {
        super(InfoType.COOPERATIVE_DISTANCE_MAP_CONSTRAINT);
        this.distances = distances;
        this.mode = mode;
        this.chain = null;
    }
    public CooperativeDistanceMapConstraintsCreator(Chain chain, String mode) {
        super(InfoType.COOPERATIVE_DISTANCE_MAP_CONSTRAINT);
        this.distances = null;
        this.mode = mode;
        this.chain = chain;
    }


    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands) {
        if (chain == null)
            return createEnergyTerm1(protein, dm, commands);
        else
            return createEnergyTerm2(protein);
    }

    public AbstractEnergy createEnergyTerm2(Protein protein) {
        FreeDistanceList distanceList = new FreeDistanceList();
        ArrayList<double[]> targets = new ArrayList<>();
        double[] targetAndRange;
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.COOPERATIVE_DISTANCE_MAP_CONSTRAINT,
                "cooperative distance map constrintnts energy", weight);
        for (int i = 1; i < chain.size(); i++) {
            Residue residue = chain.residueAt(i);
            Residue nextResidue = null;
            if (i < chain.size()-1) {
                nextResidue = chain.residueAt(i + 1);
                FreeDistance distance = new FreeDistance(residue.ca(), nextResidue.ca());
                distanceList.add(distance);
                targetAndRange = new double[2];
                targetAndRange[0] = 3.8;
                targetAndRange[1] = 0;
                targets.add(targetAndRange);
            }
            if (residue.type != ResidueType.GLY) {
                FreeDistance distance = new FreeDistance(residue.ca(), residue.cb());
                distanceList.add(distance);
                targetAndRange = new double[2];
                targetAndRange[0] = AtomBuilder.CA_CB_len;
                targetAndRange[1] = 0;
                targets.add(targetAndRange);
            }
        }
        term = new CooperativeDistanceMapConstraints(protein, distanceList, targets, energyInfoElement, mode);
        return term;
    }

    public AbstractEnergy createEnergyTerm1(Protein protein,
                                           DistanceMatrix dm, CommandList commands) {
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.COOPERATIVE_DISTANCE_MAP_CONSTRAINT,
                "cooperative distance map constrintnts energy", weight);
        CommandList myCommands = null;
        FreeDistanceList distanceList = new FreeDistanceList();
        ArrayList<double[]> targets = new ArrayList<>();
        Atom atomI, atomJ;
        for (Residue residueI : protein.chain()) {
            if ((! residueI.dummy()) && (! residueI.ca().nowhere())) {
                JSONObject residueDistances = (JSONObject) distances.get("" + residueI.number());
                if (residueDistances != null) {
                    for (Object jResidueString : residueDistances.keySet()) {
                        Residue residueJ = protein.chain().residueAt(new Integer((String) jResidueString));
                        if ((!residueJ.dummy()) && (!residueJ.ca().nowhere())) {
                            if (residueI.number() + 1 < residueJ.number()) {
                                if (mode.equals("ca")) {
                                    atomI = residueI.ca();
                                    atomJ = residueJ.ca();
                                }
                                else {
                                    atomI = residueI.cb();
                                    if (atomI == null)
                                        atomI = residueI.ca();
                                    atomJ = residueJ.cb();
                                    if (atomJ == null)
                                        atomJ = residueJ.ca();
                                }
                                FreeDistance distance = new FreeDistance(atomI, atomJ);
                                distanceList.add(distance);
                                double[] targetAndRange = new double[2];
                                targetAndRange[0] = (double) ((JSONArray) residueDistances.get(jResidueString)).get(0);
                                targetAndRange[1] = (double) ((JSONArray) residueDistances.get(jResidueString)).get(1);
                                targets.add(targetAndRange);
                            }
                        }
                    }
                } //else Utils.println("No residueDistances for " + residueI);
            }
        }
        term = new CooperativeDistanceMapConstraints(protein, distanceList, targets, energyInfoElement, mode);
        return term;
    }
}
