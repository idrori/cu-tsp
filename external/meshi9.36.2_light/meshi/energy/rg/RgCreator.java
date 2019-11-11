/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.energy.*;
import meshi.energy.rg.filters.*;
import meshi.molecularElements.*;
import meshi.geometry.*;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.*;
import meshi.util.filters.*;
import meshi.util.info.InfoType;

/**
 * An implicit solvation energy term for all-atom models, modeling a 4.0 angs solvation shell around
 * each atom.
 */

public class RgCreator extends EnergyCreator {
    private static Filter hydrophobicSideChains = new HydrophobicSideChains(true);
    private static Filter backboneFilter = new BackboneCAs();
    private static Filter polarFilter = new PolarSideChains();
    private static Filter heavyAtomsFilter = new HeavyAtomsFilter();
    private static Filter secondaryStructureFilter = new SecondaryStructureFilter();
    private static Filter coilFilter = new CoilFilter();


    public RgCreator() {
        super(InfoType.RG);
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        if (term != null) return term;

        RgInfoElement info = new RgInfoElement(weight,
                commands.firstWordFilter(RG).secondWord(RG_LOGARITMIC).thirdWordDouble(),
                commands.firstWordFilter(RG).secondWord(RG_RATIO).thirdWordDouble(),
                commands.firstWordFilter(RG).secondWord(RG_LINEAR_POLAR).thirdWordDouble(),
                commands.firstWordFilter(RG).secondWord(RG_LINEAR_NON_POLAR).thirdWordDouble(),
                commands.firstWordFilter(RG).secondWord(RG_LINEAR_BACKBONE).thirdWordDouble());


        // atom lists
        AtomList atomList = protein.atoms();
        AtomList heavyAtoms = filter(atomList, heavyAtomsFilter);
        AtomList nonPolarSS = filter(atomList, hydrophobicSideChains, secondaryStructureFilter);
        AtomList nonPolarCoil = filter(atomList,hydrophobicSideChains, coilFilter);
        AtomList backboneSS = filter(atomList,backboneFilter,secondaryStructureFilter);
        AtomList backboneCoil = filter(atomList,backboneFilter,coilFilter);
        AtomList polarSS = filter(atomList,polarFilter,secondaryStructureFilter);
        AtomList polarCoil = filter(atomList,polarFilter,coilFilter);
        if (anyIsNull(atomList,heavyAtoms,nonPolarSS,nonPolarCoil,backboneSS,backboneCoil,polarSS,polarCoil)) {
            term = new RgEnergy(info);
            return term;
        }

        //
        CenterOfMass centerOfMass = new CenterOfMass(nonPolarSS);
        CenterOfMass backboneCenterOfMass = new CenterOfMass(backboneSS);
        //
        // LogRGs
        LogRadiusOfGyration logNonPolarSsRG = new LogRadiusOfGyration(new RadiusOfGyration(nonPolarSS, centerOfMass));
        LogRadiusOfGyration logNonPolarCoilRG = new LogRadiusOfGyration(new RadiusOfGyration(nonPolarCoil, centerOfMass));
        LogRadiusOfGyration logBackboneSsRG = new LogRadiusOfGyration(new RadiusOfGyration(backboneSS, centerOfMass));
        LogRadiusOfGyration logBackboneCoilRG = new LogRadiusOfGyration(new RadiusOfGyration(backboneCoil, centerOfMass));
        LogRadiusOfGyration logPolarSsRG = new LogRadiusOfGyration(new RadiusOfGyration(polarSS, centerOfMass));
        LogRadiusOfGyration logPolarCoilRG = new LogRadiusOfGyration(new RadiusOfGyration(polarCoil, centerOfMass));
        // pack them
        LogRadiusOfGyration[] logRgArray = {logNonPolarSsRG, logNonPolarCoilRG, logBackboneSsRG, logBackboneCoilRG, logPolarSsRG, logPolarCoilRG};
        Updateable[] updateables = {centerOfMass,
                backboneCenterOfMass,
                logNonPolarSsRG,
                logNonPolarCoilRG,
                logBackboneSsRG,
                logBackboneCoilRG,
                logPolarSsRG,
                logPolarCoilRG};
        //
        //
        term = new RgEnergy(updateables, logRgArray, Math.log(heavyAtoms.size()), info);
        return term;
    }

    private AtomList filter(AtomList atoms,Filter filter){
        return filter(atoms,filter, new KolDichfin());
    }
    private AtomList filter(AtomList atoms,Filter filter1,Filter filter2){
            AtomList out = atoms.filter(filter1).filter(filter2);
            //if (out.isEmpty()) return null;
            return out;
        }

    private boolean anyIsNull(Object o1,Object o2,Object o3,Object o4,Object o5,Object o6,Object o7,Object o8){
        return (o1 == null)||(o2 == null)||(o3 == null)||(o4 == null)||(o5 == null)||(o6 == null)||(o7 == null)||(o8 == null);
    }

}