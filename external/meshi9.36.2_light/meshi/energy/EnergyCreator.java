/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.parameters.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;

import java.util.*;

import meshi.util.*;
import meshi.util.filters.*;
import meshi.energy.simpleEnergyTerms.ParametersList;

/**
 * Factory classes that generate energy terms. This is the carpet under which we hide all the boring details of
 * what is needed for the creation of a given energy term. Each energy term requires a specific subclass.
 * See meshi.energy.bond.BondCreator for a simple example.
 */
public abstract class EnergyCreator implements MeshiPotential, KeyWords {
    /**
     * The weight of the energy term within the total energy function.
     */
    protected double weight;
    protected AbstractEnergy term = null;

    public AbstractEnergy term() {
        return term;
    }

    /**
     * The keyword for the energy term. All commands relevant to this energy terms will sart with this word.
     */
    final protected InfoType infoType;
    protected boolean weightWasSet = false;


    private static String parametersDirectory = null;


    /**
     * Constructs an energy creator object. The key parameter will serve as a
     * keyword that identify relevant commands in the commands list.
     */
    public EnergyCreator(InfoType infoType) {
        this.infoType = infoType;
    }

    /**
     * Constructs an energy creator object. The key parameter will serve as a
     * keyword that identify relevant commands in the commands list.
     */
    public EnergyCreator(InfoType infoType, CommandList commands) {
        this.infoType = infoType;
        getWeight(commands);
    }

    //------------------

    /**
     * Construct a somewhat degenerate creator object that cannot read commands from the commands list.
     */
    public EnergyCreator(double weight, InfoType infoType) {
        this.weight = weight;
        weightWasSet = true;
        this.infoType = infoType;
    }

    //------------------

    /**
     * The weight of the energy term.
     */
    public double weight() {
        if (!weightWasSet)
            throw new RuntimeException(infoType + " was not setResidue");
        return weight;
    }

    //------------------

    /**
     * Extract the path of the parameters directory from the commands list.
     */
    public static String parametersDirectory(CommandList commands) {
        if (parametersDirectory == null)
            getParametersDirectory(commands);
        return parametersDirectory;
    }

    public static void getParametersDirectory(CommandList commands) {
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        parametersDirectory = command.secondWord();
    }

    //------------------

    /**
     * Where each sub class gets the chance to show what it knows. Offers a standard interface to many energy
     * functions that have different requirements.
     */
    public abstract AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                                    CommandList commands);

    /**
     * Get the weight of the energy term.
     */
    public void getWeight(CommandList commands) {
        if (!weightWasSet) {
            try {
                weight = commands.getWeight(infoType.tag);
            } catch (WeightNotFoundException ex) {
                throw new RuntimeException("getWeight of " + this + " failed\n" + ex.getMessage());
            }
            weightWasSet = true;
            Utils.println("weight = " + weight);
        }
    }

    public boolean weightWasSet() {
        return weightWasSet;
    }

    public String toString() {
        return "Creator of " + infoType.tag;
    }

    /**
     * Picking the relevant objects.
     */
    protected static class HaveParametersFilter implements Filter {
        ParametersList parametersList;

        public HaveParametersFilter(ParametersList parametersList) {
            this.parametersList = parametersList;
        }

        public boolean accept(Object obj) {
            return (parametersList.parameters(obj) != null);
        }
    }
    public void resetTerm() {
        term = null;
    }
}
