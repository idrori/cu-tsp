package meshi.energy.simpleEnergyTerms.repulsion;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParametersList;
import meshi.energy.simpleEnergyTerms.SimpleEnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;
import meshi.util.CommandList;
import meshi.util.info.InfoType;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;


public class RepulsionCreator extends SimpleEnergyCreator{
    JSONArray repulsions;
    public enum LocalPairs {O_CB, O_C, O_O, O_N, C_CB, C_C, C_O, C_N, N_O, N_N, CB_O, CB_N, CB_CB} //Todo - this is a copy of the relevant enum in the program ComntactConstraints
    public enum GPTypes {GLY, PRO, OTHER}//Todo - this is a copy of the relevant enum in the program ComntactConstraints

    public RepulsionCreator(JSONArray repulsions) {
        super(InfoType.REPULSION);
        if (repulsions == null)
            throw new RuntimeException("This is weird");
        this.repulsions = repulsions;
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands) {

        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.REPULSION,
                "Repultion energy", weight);
        Object[] allDistances;
        allDistances = repulsions.toArray();
        FreeDistanceList distanceList = new FreeDistanceList();
        DistanceConstraintsParametersList targets = new DistanceConstraintsParametersList();
        Residue prevResidue = null;
        for (Residue residue : protein.chain()) {
            if ((prevResidue != null) && (!prevResidue.dummy()) && (!prevResidue.ca().nowhere()) &
                    (! residue.dummy()) && (! residue.ca().nowhere())) {
                Object[] distances;
                if (residue.type == ResidueType.GLY)
                    distances = ((JSONArray) allDistances[GPTypes.GLY.ordinal()]).toArray();
                else if (residue.type == ResidueType.PRO)
                    distances = ((JSONArray) allDistances[GPTypes.PRO.ordinal()]).toArray();
                else distances = ((JSONArray) allDistances[GPTypes.OTHER.ordinal()]).toArray();

                FreeDistance distance;
                Atom cb1 = prevResidue.cb();
                Atom o1 = prevResidue.carbonylO();
                Atom c1 = prevResidue.carbonylC();
                Atom cb2 = residue.cb();
                Atom c2 = residue.carbonylC();
                Atom o2 = residue.carbonylO();
                Atom n2 = residue.amideN();

//                distance = new FreeDistance(o1, c2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.O_C.ordinal()]));

                if (cb2 != null) {
                    distance = new FreeDistance(o1, cb2);
                    distanceList.add(distance);
                    targets.add(new DistanceConstraintsParameters((Double)distances[LocalPairs.O_CB.ordinal()]));
                    if (cb1 != null) {
                        distance = new FreeDistance(cb1, cb2);
                        distanceList.add(distance);
                        targets.add(new DistanceConstraintsParameters((Double)distances[LocalPairs.CB_CB.ordinal()]));

                    }
                }

//                distance = new FreeDistance(o1, n2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.O_N.ordinal()]));
//
//                distance = new FreeDistance(o1, o2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.O_O.ordinal()]));
//
//                distance = new FreeDistance(c1, c2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.C_C.ordinal()]));

                if (cb2 != null) {
                    distance = new FreeDistance(c1, cb2);
                    distanceList.add(distance);
                    targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.C_CB.ordinal()]));
                }

//                distance = new FreeDistance(c1, n2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.C_N.ordinal()]));
//
//                distance = new FreeDistance(c1, o2);
//                distanceList.add(distance);
//                targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.C_O.ordinal()]));

                if (cb2 != null) {
                    distance = new FreeDistance(cb2, n2);
                    distanceList.add(distance);
                    targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.CB_N.ordinal()]));
                    distance = new FreeDistance(cb2, o2);
                    distanceList.add(distance);
                    targets.add(new DistanceConstraintsParameters((Double) distances[LocalPairs.CB_O.ordinal()]));
                }
            }
            prevResidue = residue;
           }
        term = new Repulsion(distanceList, targets, energyInfoElement);
        return term;
    }
}
