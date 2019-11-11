package meshi.energy.simpleEnergyTerms.distanceMapConstraints;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.cooperativeDistanceMapConstraints.CooperativeDistanceMapConstraints;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraints;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParametersList;
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


public class DistanceMapConstraintsCreator extends SimpleEnergyCreator{
    private JSONObject distances;
    private String mode;
    private Chain chain;
    public DistanceMapConstraintsCreator(JSONObject distances, String mode) {
        super(InfoType.DISTANCE_MAP_CONSTRAINT);
        this.distances = distances;
        this.mode = mode;
        chain = null;
    }

    public DistanceMapConstraintsCreator(Chain chain, String mode) {
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
        DistanceConstraintsParametersList targets = new DistanceConstraintsParametersList();
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
                targets.add(new DistanceConstraintsParameters(targetAndRange));
            }
            if (residue.type != ResidueType.GLY) {
                FreeDistance distance = new FreeDistance(residue.ca(), residue.cb());
                distanceList.add(distance);
                targetAndRange = new double[2];
                targetAndRange[0] = AtomBuilder.CA_CB_len;
                targetAndRange[1] = 0;
                targets.add(new DistanceConstraintsParameters(targetAndRange));
            }
        }
        term = new DistanceMapConstraints(distanceList, targets, energyInfoElement);
        return term;
    }

    public AbstractEnergy createEnergyTerm1(Protein protein,
                                           DistanceMatrix dm, CommandList commands) {
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.DISTANCE_MAP_CONSTRAINT,
                "distance map constrintnts energy", weight);
        CommandList myCommands = null;
        FreeDistanceList distanceList = new FreeDistanceList();
        DistanceConstraintsParametersList targets = new DistanceConstraintsParametersList();
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
                                double[] target = new double[2];
                                target[0] = (double) ((JSONArray) residueDistances.get(jResidueString)).get(0);
                                target[1] = (double) ((JSONArray) residueDistances.get(jResidueString)).get(1);
                                targets.add(new DistanceConstraintsParameters(target));
                                if ((distance.distance > 50) && (target[0] > 50) && (distance.distance - target[0] > 10)) {
                                    Atom atom = residueJ.amideN();
                                    distanceList.add(new FreeDistance(atomI, atom));
                                    targets.add(new DistanceConstraintsParameters(target));
                                    atom = residueJ.carbonylC();
                                    distanceList.add(new FreeDistance(atomI, atom));
                                    targets.add(new DistanceConstraintsParameters(target));
                                    atom = residueJ.amideH();
                                    if ((atom != null) && (!atom.nowhere())) {
                                        distanceList.add(new FreeDistance(atomI, atom));
                                        targets.add(new DistanceConstraintsParameters(target));
                                    }
                                    atom = residueJ.carbonylO();
                                    if((atom != null) &&  (!atom.nowhere())){
                                        distanceList.add(new FreeDistance(atomI, atom));
                                        targets.add(new DistanceConstraintsParameters(target));
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    int i = 0;
                    for (Residue residueJ : protein.chain()) {
                        if ((!residueJ.dummy()) && ((Math.abs(residueI.number() - residueJ.number()) >= 4) & (!residueJ.ca().nowhere()))) {
                            FreeDistance distance = new FreeDistance(residueI.ca(), residueJ.ca());
                            distanceList.add(distance);
                            double[] target = {-10,0};
                            targets.add(new DistanceConstraintsParameters(target));
                        }
                        i++;
                    }
                }
            }
        }
        term = new DistanceMapConstraints(distanceList, targets, energyInfoElement);
        return term;
    }
}
