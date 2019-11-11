/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms;

import meshi.energy.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.ResidueData;
import meshi.util.Updateable;
import meshi.util.info.ChainsInfo;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A super class for those energy terms that operate isOn a fixed list of elements.
 * Typical cases are the bonded terms (bond angle etc.), where the list of elements
 * (say covalent bonds)
 * is fixed by the chemical structure of the molecule.
 */
public abstract class SimpleEnergyTerm extends AbstractEnergy implements EvaluatesResidues{

    protected ArrayList<EnergyElement> elementsList;

    /**
     * A list of all the parameters needed by the energy term.
     */
    protected ParametersList parametersList;

    public SimpleEnergyTerm() {
    }

    public SimpleEnergyTerm(Updateable[] updateableResources,
                            ParametersList parametersList,
                            EnergyInfoElement info, EnergyType type) {
        super(updateableResources, info,type);
        this.parametersList = parametersList;
    }
     public SimpleEnergyTerm(Updateable[] updateableResources,
                            ParametersList parametersList,
                            EnergyInfoElement info) {
        this(updateableResources, parametersList, info, EnergyType.DIFFERENTIAL);
    }



    /**
     * Creates the fixed list of energy elements in the initiation phase of a simulation.
     * These elements typically warp some object such as AtomPair, Angle or Torsion.
     * The warped elements are provided by the parameter "baseList".
     * Note that Different energy terms treat this list differently.
     * BondEnergy for example tries to warp every single element of "baseList" while PlaneEnergy
     * assumes that the list is redundant and warps only those Torsion elements for which it finds parameters.
     * The binding of "baseElement" to the parameters is done by the method parameters of ParametersList.
     */
    public void createElementsList(ArrayList baseList) {
        elementsList = new ArrayList();

        Parameters parameters;
        for (Object baseElement : baseList) {
            parameters = parametersList.parameters(baseElement);
            if (parameters == null) handleMissingParameters(baseElement);
            else {
                EnergyElement newElement = createElement(baseElement, parameters);
                if ((newElement != null) && (!newElement.frozen()))
                    elementsList.add(newElement);
            }
        }
    }


    public boolean evaluatesResidues() { return true; }

    public ArrayList elementsList() {
        return elementsList;
    }

    public abstract EnergyElement createElement(Object baseElement, Parameters parameters);

    /**
     * Testing of one atom in all energy elements
     *
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom        an criminal <code>Atom</code> value
     */
    public void test(TotalEnergy totalEnergy, Atom atom) throws EvaluationException{
        if (!on) System.out.println("" + this + " is off");
        for (Iterator eei = elementsList.iterator(); eei.hasNext();) {
            EnergyElement energyElement = (EnergyElement) eei.next();
            if (!energyElement.frozen()) {
                energyElement.test(totalEnergy, atom);
            }
        }
    }


    /**
     * Evaluates the current value of the energy function and <b> update </b> the derivatives.
     */
    public EnergyInfoElement evaluate() throws EvaluationException {
        evaluateResidues(null);
        return info;
    }


        //for debug
//        EnergyElement hotest=null;
//        double highestEnergy = -999999;

    public void evaluateResidues(ChainsInfo chainsInfo) {
        if (!on) {
            info.setValue(0);
            return;
        }

        double e, energy = 0;
        ResidueData residueData = null;
        if (chainsInfo != null) {
            residueData = new ResidueData(chainsInfo);
        }
        if (elementsList == null)
            throw new RuntimeException("this is very weird");
        for (EnergyElement energyElement : elementsList) {
            e = energyElement.evaluate();
            energy += e;
            if (residueData != null) {
                int nAtoms = energyElement.atoms().size();
                for (Atom atom :energyElement.atoms()) {
                    residueData.add(atom,e/nAtoms);
                }
            }
//            if (highestEnergy < e) {
//                highestEnergy = e;
//                hotest = energyElement;
//            }
        }
        if (chainsInfo != null) {
            chainsInfo.add(residueData,info.type);
        }
        if (info == null) throw new RuntimeException("null info in " + this);
        info.setValue(energy);            //*36
//        if (highestEnergy > 100)
//            System.out.println("\n"+hotest+"\n"+highestEnergy+"\n---- debug -----\n");
    }

    public void evaluateAtoms() throws EvaluationException{
        if (on) {
            for (EnergyElement energyElement : elementsList) {
                energyElement.evaluateAtoms();
            }
        }
    }
}


  
