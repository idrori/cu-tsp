/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.beta.peptideBonds.PeptideBond;
import meshi.molecularElements.BetaSegmentList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.Segment;
import meshi.molecularElements.atoms.Atom;

import java.util.ArrayList;

/**
 *
 */
public class BetaEnergy extends AbstractEnergy {
    public static final int NUMBER_OF_DISTANCES = 5;
    public static final double TARGET_DISTANCE = 2;
    private ArrayList<BetaEnergyElement> elements;
    private BetaSegmentList segments;
    private Protein protein;

    public BetaEnergy() {
    }

    public BetaEnergy(Protein protein, EnergyInfoElement info) {
        super(toArray(), info);
        comment = "Beta";
        this.protein = protein;
        PeptideBond.attachPeptideBonds(protein);
        segments = new BetaSegmentList(protein);
        reset();
    }

    public void reset() {
        elements = new ArrayList<BetaEnergyElement>();
        for (Segment segment : segments)
            for (Residue residue : segment)
                if (residue.peptideBond() != null)
                    residue.peptideBond().neighbors(segments);
        for (Segment segment : segments)
            for (Residue residue : segment)
                elements.add(new BetaEnergyElement(residue, segments, weight));
    }

    public void off() {
        for (BetaEnergyElement element : elements)
            element.off();
    }

    public void on() {
        for (BetaEnergyElement element : elements)
            element.on();
    }

    public EnergyInfoElement evaluate() {

        double energy = 0;
        for (BetaEnergyElement element : elements)
            energy += element.evaluate();

        info.setValue(energy);      //*11
        return info;
    }

    public void evaluateAtoms() {
        for (BetaEnergyElement element : elements)
            element.evaluateAtoms();
    }

    public void test(TotalEnergy energy, Atom atom)  {
        for (BetaEnergyElement element : elements)
            if (element.isOn()) element.test(energy, atom);

    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        reset();
        super.setNumberOfUpdates(numberOfUpdates);
    }

    public void scale(double factor) {
        weight = weight * factor;
    }


}




