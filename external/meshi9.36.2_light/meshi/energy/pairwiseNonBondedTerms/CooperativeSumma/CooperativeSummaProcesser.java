/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomCore;
import meshi.util.MeshiAttribute;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 15/02/2010
 * Time: 12:16:41
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeSummaProcesser {

    private final DistanceMatrix distanceMatrix;
    //private DistanceList nonBondedList;
//for cooperative energy terms
    protected  double[] atomEnergies;
    protected  double[] atomEnergiesPOLARS;
    protected  double[] atomEnergiesNonPOLARS;
    protected  double[] atomEnergiesNEUTRALS;
    //protected static double[] atomEnergiesPOLARSH;
    protected  double[] atomEnergiesPOLARSnoH;
    protected  double[] atomEnergiesPOLARSnoH_NNorOO;

    public double[] atomEnergies() {
        return atomEnergies;
    }

    public double[] atomEnergiesPOLARS() {
        return atomEnergiesPOLARS;
    }

    public double[] atomEnergiesNonPOLARS() {
        return atomEnergiesNonPOLARS;
    }

    public double[] atomEnergiesNEUTRALS() {
        return atomEnergiesNEUTRALS;
    }
    //public static double [] atomEnergiesPOLARSH(){return atomEnergiesPOLARSH;}

    public double[] atomEnergiesPOLARSnoH() {
        return atomEnergiesPOLARSnoH;
    }

    public double[] atomEnergiesPOLARSnoH_NNorOO() {
        return atomEnergiesPOLARSnoH_NNorOO;
    }


    public CooperativeSummaProcesser(DistanceMatrix distanceMatrix) {
        this.distanceMatrix = distanceMatrix;
        atomEnergies = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARS = new double[distanceMatrix.molecularSystem.size()];
        //atomEnergiesPOLARSH = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARSnoH = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARSnoH_NNorOO = new double[distanceMatrix.molecularSystem.size()];

        atomEnergiesNonPOLARS = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesNEUTRALS = new double[distanceMatrix.molecularSystem.size()];

    }

    public void reset() {
        for (int i = 0; i < atomEnergiesPOLARS.length; i++) {
            atomEnergies[i] = 0;
            atomEnergiesPOLARS[i] = 0;
            //       atomEnergiesPOLARSH[i] = 0;
            atomEnergiesPOLARSnoH[i] = 0;
            atomEnergiesPOLARSnoH_NNorOO[i] = 0;

            atomEnergiesNonPOLARS[i] = 0;
            atomEnergiesNEUTRALS[i] = 0;
        }
        /*
DistanceList nonBondedList = distanceMatrix.nonBondedList();
for (Distance distance:nonBondedList) {
if (distance.mode().frozen) continue;
summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
if (summaAttribute != null)
                 summaAttribute.restart();
}
        */
    }

    public void setCooperativeStatistic(Distance distance, double fx, double fy, double fz, double halfEnergy0) {
        SummaAttribute summaAttribute;
        summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
        if (summaAttribute == null) {
            summaAttribute = new SummaAttribute();
            distance.addAttribute(summaAttribute);
        }

        /* apply force to atoms */
        summaAttribute.fx = fx;
        summaAttribute.fy = fy;
        summaAttribute.fz = fz;

/*
//For statistic

        atomEnergies[atom1.number] += halfEnergy0;
        atomEnergies[atom2.number] += halfEnergy0;
     //for atom1
     //  if (atom2.type().isPolar() || atom2.type().isPolarSideChains()) {
          if (atom2.type().isPolar() ) {
           atomEnergiesPOLARS[atom1.number] += halfEnergy0;
       }
       else
         if (atom2.type().isNonPolar()) atomEnergiesNonPOLARS[atom1.number] += halfEnergy0;
         else
             if (atom2.type().isNeutral()) atomEnergiesNEUTRALS[atom1.number] += halfEnergy0;

     //for atom2
       if (atom1.type().isPolar() || atom1.type().isPolarSideChains()) {
         atomEnergiesPOLARS[atom2.number] += halfEnergy0;
     }
     else
       if (atom1.type().isNonPolar()) atomEnergiesNonPOLARS[atom2.number] += halfEnergy0;
       else
           if (atom1.type().isNeutral()) atomEnergiesNEUTRALS[atom2.number] += halfEnergy0;
          //*/

//--------------------------------------------------For statistic with H bond--------------------------------------------------
        AtomCore atom1, atom2;

        atom1 = distance.atom1;
        atom2 = distance.atom2;
        if (atom1.type().isOther() || atom2.type().isOther())
            return;

        atomEnergies[atom1.number] += halfEnergy0;
        atomEnergies[atom2.number] += halfEnergy0;

        if ((!atom1.type().isPolarBackbone()) || (!atom2.type().isPolarBackbone())) {
            //for atom1
            if (atom2.type().isPolar())
                atomEnergiesPOLARS[atom1.number] += halfEnergy0;
            else if (atom2.type().isNonPolar())
                atomEnergiesNonPOLARS[atom1.number] += halfEnergy0;
            else if (atom2.type().isNeutral())
                atomEnergiesNEUTRALS[atom1.number] += halfEnergy0;
            else throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this);
            //for atom2
            if (atom1.type().isPolar())
                atomEnergiesPOLARS[atom2.number] += halfEnergy0;
            else if (atom1.type().isNonPolar())
                atomEnergiesNonPOLARS[atom2.number] += halfEnergy0;
            else if (atom1.type().isNeutral())
                atomEnergiesNEUTRALS[atom2.number] += halfEnergy0;
            else throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this);
        } else {
            if (!(atom1.type().isPolarBackbone() && atom2.type().isPolarBackbone()))
                throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this);
//only two Polars N or O are stayed

            if ((atom1.type().isOxygen() && atom2.type().isOxygen())
                    ||
                    (atom1.type().isNitrogen() && atom2.type().isNitrogen())) {
                atomEnergiesPOLARSnoH_NNorOO[atom1.number] += halfEnergy0;
                atomEnergiesPOLARSnoH_NNorOO[atom2.number] += halfEnergy0;
            } else {  //hBond
//                atomEnergiesPOLARS[atom1.number] += halfEnergy0;
                //              atomEnergiesPOLARS[atom2.number] += halfEnergy0;

                if (Math.abs(atom1.atom.residue().number() - atom2.atom.residue().number()) < 3)
                    return;
                atomEnergiesPOLARSnoH[atom1.number] += halfEnergy0;
                atomEnergiesPOLARSnoH[atom2.number] += halfEnergy0;

                if (!((atom1.type().isNitrogen() && atom2.type().isOxygen()) ||
                        (atom2.type().isNitrogen() && atom1.type().isOxygen())))
                    throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this);
            }
        }  // end of two Backbone Polars
    }

}
