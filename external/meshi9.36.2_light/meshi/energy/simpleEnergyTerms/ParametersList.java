/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms;

import meshi.energy.*;
import meshi.util.Utils;
import meshi.util.file.*;

import java.util.*;

/**
 * A list of parameters for an energy term.
 * Energy terms (extensions of AbstractEnergy) typically need a large number of parameters. These
 * parameters are provided by the user in a text file. ParametersList objects read, parse and stores
 * the contennts of these files. They provide the parameters to the energy term with the
 * getParameters(Parameters key) method. In general, each energy term requires a specific
 * ParametersList class. See for example meshi.energy.bond.BondParametersList.
 */
public abstract class ParametersList extends ArrayList<Parameters> {
    /**
     * True if the list is sortable.
     */
    boolean sortable;

    /**
     * True if the list is sorted.
     */
    boolean sorted = false;


    public ParametersList() {
        super();
    }

    /**
     * Construct a ParametersList object from parameters file. It is assumed that the file is arranged
     * in raws where each raw includes the parameters for a specific interaction. the # sign may appear
     * anywhere in a line and indicates the beginning of a comment. If the parameters are sortable the list
     * is sorted. See for example the parameters for
     * bond-energy in meshi/parameters/meshiPotential/bondParameters.dat .
     */
    public ParametersList(String parametersFileName, boolean sortable) {
        super();
        this.sortable = sortable;
        if (parametersFileName != null) {
            Utils.println("Loading parameters from " + parametersFileName);
            try {
                MeshiLineReader lines = new MeshiLineReader(parametersFileName);
                String line;
                while ((line = lines.readLine("#")) != null) {
                    add(createParameters(line));
                }
                if (sortable) sort();
            }
            catch (RuntimeException e) {
                System.out.println("A problem while reading parameters file " +
                        parametersFileName);
                throw e;
            }
        }
    }


    /**
     * Construct a ParametersList object from multiple files.
     * This constructor is useful when the parameters are distributed over many files,
     * where each file holds a single parameter. An example for this case are the TwoTorsions
     * energies.
     */
    public ParametersList(String[] parametersFileName, boolean sortable) {
        super();

        this.sortable = sortable;
        if (parametersFileName != null) {
            if (parametersFileName.length > 0) {
                Utils.println("Loading " + this + " parameters");
                for (int cc = 0; cc < parametersFileName.length; cc++) {
                    add(createParameters(parametersFileName[cc]));
                    if (sortable) sort();
                }
            }
        }
    }


    /**
     * Adds an element (must be an instance of Parameters) to the list.
     */
    public boolean add(Parameters element) {
        if (sorted) throw new RuntimeException("Cannot add to ParametersList after sorting");
        return super.add(element);
    }

    /**
     * Sort the list. Note that if the List is not sortable this method simply does nothing.
     */
    public void sort() {
        if (!sortable) return;
        Parameters[] array = toArray(new Parameters[size()]);
        Arrays.sort(array);
        clear();
        for (Parameters p : array)
            add(p);
        sorted = true;
    }

    public Iterator iterator() {
        if (sortable & (!sorted))
            throw new RuntimeException("Sortable ParametersList not sorted yet.\n" +
                    "Cannot provide an iterator");
        return super.iterator();
    }

    /**
     * <pre>
     * Fetches a parameter from the list. A sorted list is assumed and a binary search is performed.
     * If the list is not sortable this method needs to be overrun. Example (from BondParametersList):
     * <p/>
     * 	Parameters key = new BondParameters(pair.largeType(), pair.smallType());
     *  return getParameters(key);
     * </pre>
     */
    public Parameters getParameters(Parameters key) {
        if (sortable) {
            int index = Arrays.binarySearch(toArray(new Parameters[size()]), key);
            if (index < 0) return null;
            return (Parameters) get(index);
        }
        throw new RuntimeException("The generic getParameters(Parameters key, \n" +
                "                          ParametersList parametersList)\n" +
                "method of ParametersList assume sorable parameters and \n" +
                "uses binary search. Apparently this approach is not \n" +
                "applicable here.");
    }

    /**
     * Energy term specific method to fetch parameters for the interactions between a st of atoms.
     */
    public abstract Parameters parameters(Object Obj);


    /**
     * Energy term specific method to create a Parameters object from a line of the parameters file.
     */
    public abstract Parameters createParameters(String line);

    /**
     * Get an element from the list.
     **/

}

    
    
    
