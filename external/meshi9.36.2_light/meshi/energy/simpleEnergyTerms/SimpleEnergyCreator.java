/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms;

import meshi.energy.*;
import meshi.parameters.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;

import java.util.*;

import meshi.util.*;
import meshi.util.filters.*;

/**
 * Factory classes that generate energy terms. This is the carpet under which we hide all the boring details of
 * what is needed for the creation of a given energy term. Each energy term requires a specific subclass.
 * See meshi.energy.bond.BondCreator for a simple example.
 */
public abstract class SimpleEnergyCreator extends EnergyCreator {
    /**
     * A list of the parameters needed for this energy term.
     */
    protected ParametersList parametersList = null;

    /**
     * get method for parametersList.
     */
    public ParametersList parametersList() {
        return parametersList;
    }

    private static String parametersDirectory = null;


    /**
     * Constructs an energy creator object. The key parameter will serve as a
     * keyword that identify relevant commands in the commands list.
     */
    public SimpleEnergyCreator(InfoType infoType) {
        super(infoType);
    }

    /**
     * Constructs an energy creator object. The key parameter will serve as a
     * keyword that identify relevant commands in the commands list.
     */
    public SimpleEnergyCreator(InfoType infoType, CommandList commands) {
        super(infoType, commands);
    }

    //------------------

    /**
     * Construct a somewhat degenerate creator object that cannot read commands from the commands list.
     */



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
}
