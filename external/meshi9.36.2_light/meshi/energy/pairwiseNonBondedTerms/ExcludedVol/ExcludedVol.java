/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

/**
 *Excluded Volume potential.
 *
 *A potential that creates a strong repulsion between atoms when they are getting
 *near their VDW radius. There is no attraction part like in the VDW.
 *The functional form of the term is:  
 *
 *dis =                EV
 *
 *[0,sigma]  C*(dis-sigma)^4
 *[sigma,Inf]          0
 *
 *width is setResidue in ExcludedVolParameters.java. Currently it is 0.2 Ang.
 *width is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
 *
 **/
package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceList;
import meshi.geometry.Distance;
import meshi.util.Classes;
import meshi.util.Utils;
import meshi.util.filters.Filter;

public class ExcludedVol extends NonBondedEnergyTerm implements Classes {
    Filter filter;


    public ExcludedVol(){}
    public ExcludedVol(DistanceMatrix distanceMatrix,
                       ExcludedVolParametersList parametersList,
                       int type,
                       EnergyInfoElement info,
                       double Rfac,
                       Filter filter) {
        super(toArray(distanceMatrix), info, distanceMatrix);
        comment = "ExcludedVol";
        energyElement = new ExcludedVolEnergyElement(parametersList, distanceMatrix, type, weight, Rfac);
        this.filter = filter;
    }


    public void scaleWeight(double factor) {energyElement.scaleWeight(factor);}
    public EnergyInfoElement evaluate() throws EvaluationException{
        if (!on) {
            info.setValue(0);
            return info;
        }
        double energy = 0;
        DistanceLists nonBondedLists;
        if (filter == null)
            nonBondedLists = distanceMatrix.nonBondedList();
        else {
            nonBondedLists = new DistanceLists(100);
            Utils.filter(distanceMatrix.nonBondedList(), filter, nonBondedLists);
        }
        if (nonBondedLists == null)
            System.out.println("!!!!!!");

        for (DistanceList distanceList : nonBondedLists) {
            for (Distance distance : distanceList) {
                if (distance != null && !distance.mode().frozen) {
                    energyElement.set(distance);
                    double e = energyElement.evaluate();
//                if (e > 1000)
//                    throw new RuntimeException("*************** debug ************  "+distance.distance+"\n"+
//                                               distance.atom1()+"\n"+ distance.atom2()+"\n");
                    energy += e;
                }
            }
        }
        info.setValue(energy);
        return info;
    }//evaluate
}

	
 
