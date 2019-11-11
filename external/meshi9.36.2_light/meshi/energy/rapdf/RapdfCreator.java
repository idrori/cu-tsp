package meshi.energy.rapdf;
import meshi.energy.simpleEnergyTerms.*;

import meshi.energy.EnergyCreator;
import meshi.energy.rapdf.*;
import java.util.*;
import meshi.util.*;
import meshi.util.filters.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.*;

import meshi.geometry.DistanceMatrix;

/**
 * Created by N. Vanetik
 */

import meshi.energy.*;

import meshi.molecularElements.*;
import meshi.util.info.InfoType;


public class RapdfCreator extends SimpleEnergyCreator  implements KeyWords {
    private String comment, rapdfFileNameKey;
    private final InfoType infoType;
    public RapdfCreator() {
        super(InfoType.RAPDF);
        this.comment = "Rapdf";
        this.rapdfFileNameKey = "rapdfoutputfile";
        infoType = InfoType.RAPDF;
    }
        public RapdfCreator(String comment, String rapdfFileNameKey) {
        super(InfoType.RAPDF1);
            this.comment = comment;
            this.rapdfFileNameKey = rapdfFileNameKey;
            infoType = InfoType.RAPDF1;
    }

  public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                               CommandList commands) {
                 Utils.println("Creating RapdfCreator energy term!");
        
        if (term != null) return term;

        AtomList atomList = protein.atoms();
        double[][] initialLocationsMatrix=new double[3][atomList.size()];
        for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (! atom.nowhere()) {
                int number = atom.number();
                try {
                    initialLocationsMatrix[0][number]=atom.x();
                    initialLocationsMatrix[1][number]=atom.y();
                    initialLocationsMatrix[2][number]=atom.z();
                }
                catch (RuntimeException ex) {
                    Utils.println("Problem in RapdfCreator.createEnergyTerm while processing atom number :"+number+"\n"+atom);
                    throw ex;
                }
            }
        }
        EnergyInfoElement info = new EnergyInfoElement(infoType, comment+"  energy",weight);
        term = new  RapdfEnergy(atomList, initialLocationsMatrix,info,comment, rapdfFileNameKey);
        ((RapdfEnergy)term).set(commands);
        //term.off();
        return term;
    }
}

