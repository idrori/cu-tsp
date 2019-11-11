/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomCore;
import meshi.geometry.DistanceMatrix;
import meshi.parameters.AtomType;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 22/12/2009
 * Time: 13:24:41
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZ5typesSummaCreator extends EnergyCreator implements KeyWords {
    static CooperativeZ5typesParameters parameters = null;
    AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator;
    private static CooperativeSummaProcesser cooperativeSummaProcesser;
    private CooperativeZSummaPolarElement polarElement;
    private CooperativeZSummaPolarElement nonPolarElement;
    private CooperativeZSummaPolarElement neutralElement;
    private CooperativeZSummaPolarElement polarNN_OOElement;
    private CooperativeZSummaPolarElement polarBBElement;


    public CooperativeZ5typesSummaCreator(AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator) {
        super(InfoType.COOPERATIVE_Z_SUMMA);
        this.atomicPairwisePMFSummaCreator = atomicPairwisePMFSummaCreator;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        /* load parameters if this is the first time the function is called */
        if (term != null) return term;
        if ((parameters == null) | (weight() == 0))
            parameters = new CooperativeZ5typesParameters(commands);

        cooperativeSummaProcesser = AtomicPairwisePMFSumma.cooperativeSummaProcesser;
        if (cooperativeSummaProcesser == null)
            throw new RuntimeException("The AtomicPairwisePMFSumma data have not been processed for " + this);

        int nAtoms;
        int nAtomsPolar;
        int nAtomsNonPolar;
        int nAtomsNeutral;
        int nPolarBackboneNO;
        double nonPolarWeight = commands.firstWordFilter(COOPERATIVE_Z_SUMMA_ENERGY).secondWord(COOPERATIVE_Z_SUMMA_ENERGY_NON_POLAR).thirdWordDouble();
        double polarWeight = commands.firstWordFilter(COOPERATIVE_Z_SUMMA_ENERGY).secondWord(COOPERATIVE_Z_SUMMA_ENERGY_POLAR).thirdWordDouble();
        double neutralWeight = commands.firstWordFilter(COOPERATIVE_Z_SUMMA_ENERGY).secondWord(COOPERATIVE_Z_SUMMA_ENERGY_NEUTRAL).thirdWordDouble();
        double backboneWeightNNOO = commands.firstWordFilter(COOPERATIVE_Z_SUMMA_ENERGY).secondWord(COOPERATIVE_Z_SUMMA_ENERGY_BACKBONE_OO_NN).thirdWordDouble();
        double backboneWeight = commands.firstWordFilter(COOPERATIVE_Z_SUMMA_ENERGY).secondWord(COOPERATIVE_Z_SUMMA_ENERGY_BACKBONE).thirdWordDouble();


        nAtoms = distanceMatrix.molecularSystem.size();
        AtomType type;
        nAtomsPolar = nAtomsNonPolar = nAtomsNeutral = nPolarBackboneNO = 0;
        int nOther = 0;
        for (AtomCore atom : distanceMatrix.molecularSystem) {
            type = atom.type();
            if (type.isPolar() || type.isPolarSideChains()) {
                nAtomsPolar++;
                if (type.isPolarBackbone() && (type.isOxygen() || type.isNitrogen()))
                    nPolarBackboneNO++;
            } else if (type.isNonPolar())
                nAtomsNonPolar++;
            else if (type.isNeutral())
                nAtomsNeutral++;
            else nOther++;

        }
        nAtoms = distanceMatrix.molecularSystem.size();
        if ((nAtomsPolar + nAtomsNonPolar + nAtomsNeutral + nOther) != nAtoms)
            throw new RuntimeException("Something weird with atom.type() Polarity assignment in " + this);

        nAtomsPolar = nAtoms - nOther;
        nAtomsNonPolar = nAtoms - nOther;
        nAtomsNeutral = nAtoms - nOther;

        polarElement = new CooperativeZSummaPolarElement(parameters.meanPolar, parameters.stdPolar, parameters.proportionPolar, nAtomsPolar, cooperativeSummaProcesser.atomEnergiesPOLARS, weight, distanceMatrix, "Polar", polarWeight);
        nonPolarElement = new CooperativeZSummaPolarElement(parameters.meanNonPolar, parameters.stdNonPolar, parameters.proportionNonPolar, nAtomsNonPolar, cooperativeSummaProcesser.atomEnergiesNonPOLARS, weight, distanceMatrix, "NonPolar", nonPolarWeight);
        neutralElement = new CooperativeZSummaPolarElement(parameters.meanNeutral, parameters.stdNeutral, parameters.proportionNeutral, nAtomsNeutral, cooperativeSummaProcesser.atomEnergiesNEUTRALS, weight, distanceMatrix, "Neutrals", neutralWeight);
// Polar backbone
        polarNN_OOElement = new CooperativeZSummaPolarElement(parameters.meanPOLARSnoH_NNorOO, parameters.stdPOLARSnoH_NNorOO, parameters.proportionPOLARSnoH_NNorOO, nPolarBackboneNO, cooperativeSummaProcesser.atomEnergiesPOLARSnoH_NNorOO, weight, distanceMatrix, "POLARS_BB_NNorOO", backboneWeightNNOO);
        polarBBElement = new CooperativeZSummaPolarElement(parameters.meanPOLARS_BB, parameters.stdPOLARS_BB, parameters.proportionPOLARS_BB, nPolarBackboneNO, cooperativeSummaProcesser.atomEnergiesPOLARSnoH, weight, distanceMatrix, "POLARS_BB", backboneWeight);
        //polarBBElement = null;
        EnergyInfoElement info = new CooperativeSummaInfoElement(InfoType.COOPERATIVE_Z_SUMMA, "Cooperative term that affects the pairwise term of Summan and Levitt.",
                new EnergyInfoElement(InfoType.COOPERATIVE_SUMMA_POLAR, "Polar element of Cooperative Summa term", polarWeight),
                new EnergyInfoElement(InfoType.COOPERATIVE_SUMMA_NON_POLAR, "Polar element of Cooperative Summa term", nonPolarWeight),
                new EnergyInfoElement(InfoType.COOPERATIVE_SUMMA_NEUTRAL, "Neutral element of Cooperative Summa term", neutralWeight),
                new EnergyInfoElement(InfoType.COOPERATIVE_SUMMA_POLAR_NN_OO, "Backbone (no H-bond) element of Cooperative Summa term", backboneWeightNNOO),
                new EnergyInfoElement(InfoType.COOPERATIVE_SUMMA_POLAR_BB, "Backbone H-bonds element of Cooperative Summa term", backboneWeight), weight);

        return term = new CooperativeSumma5PolarTypesEnergy(distanceMatrix, info,
                (AtomicPairwisePMFSumma) atomicPairwisePMFSummaCreator.term(),
                polarElement,
                nonPolarElement,
                neutralElement,
                polarNN_OOElement,
                polarBBElement,
                "CmmZ5t_Summa"
        );
    }


}

