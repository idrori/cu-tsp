/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.util.*;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * A super class for all meshi energy terms.
 * Currently oriented towards drivable energy functions. See meshi.energy.bond.BondEnergy.
 */
public abstract class AbstractEnergy implements Updateable, Attributable {
    /**
     * Energy function/term name.
     **/
    /**
     * A short name for the energy term.
     */
    protected String comment = "Abstract energy";

    protected Object debug = null;
    public void setDebug(Object obj) {debug = obj;}
    public void resetDebug() {debug = null;}

    public boolean hasTheSameDistanceMatrix(DistanceMatrix distanceMatrix) {
        for (Updateable updateable : updateableResources)
            if (updateable == distanceMatrix) return true;
        return false;
    }
    /**
     * The weight of the energy term within the total energy.
     */
    protected double weight;
    public double weight() { return weight;}
    public void setWeight(double weight) {this.weight = weight;}
    public void scaleWeight(double factor) {
        weight *= factor;
        Utils.print(toString()+" weight was scaled by "+factor+"  to "+weight+"\n");
    }

    /**
     * The energy term is evaluated only if it is isOn.
     */
    protected boolean on = true;


    public static final Filter filter = new IsEnergy();

    /**
     * 1/0 good for dubugging.
     */
    public static final double INFINITY = 1 / 0.0;
    /**
     * sqrt(-1) good for dubugging.
     */
    public static final double NaN = Math.sqrt(-1.0);


    /*
    * A list of all the updateable resources needed by the energy term. For a detailed description
    * of MESHI's treatment of updatable resources see the documentation of TotalEnergy.
    **/
    protected UpdateableList updateableResources;



    /**
     * construct An energy term.
     */
    protected EnergyInfoElement info;

    public EnergyInfoElement info() {
        return info;
    }

    public final EnergyType type;

    public AbstractEnergy(Updateable[] updateableResources,
                                                     EnergyInfoElement info, EnergyType type) {
        this.updateableResources = new UpdateableList(updateableResources);
        this.info = info;
        weight = info.weight();
        this.type =  type;
    }

public AbstractEnergy(Updateable[] updateableResources,
                          EnergyInfoElement info) {
        this(updateableResources, info, EnergyType.DIFFERENTIAL);
    }

    /**
       * Construct a dummy object that cannot be evaluated. Such an object is useful as a key
       * during searches.
       */
      public AbstractEnergy() {
          type = EnergyType.UNKNOWN;
      }


    /**
     * Updates the updatable resources. For a detailed description
     * of MESHI's treatment of updatable resources see the documentation of TotalEnergy.
     */
    public void update(int numberOfUpdates) throws UpdateableException {
        updateableResources.update(numberOfUpdates);
    }

    /**
     * Evaluates the energy term and <b> update </b> the derivatives.
     */
    public abstract EnergyInfoElement evaluate() throws EvaluationException;

    public EnergyInfoElement genericEvaluate() throws EvaluationException{
        if (info == null) throw new RuntimeException(this + " has null info");
        info.reset();
        if (!on) return info;
        return evaluate();
    }

    public boolean differential() {
        return type == EnergyType.DIFFERENTIAL;
    }

    /**
     * Evaluates the energy term and devides the energy between the atoms. The energy field of
     * each atom is assigned a value - its contribution to the total energy sum.
     */
    public abstract void evaluateAtoms() throws EvaluationException;

    /**
     * Turnes the energyTerm ON.
     */
    public void on() {
        on = true;
    }

    /**
     * Turnes the energyTerm OFF.
     */
    public void off() {
        on = false;
    }

    public boolean isOn() {
        return on;
    }

    /**
     * Looking for one "criminal" atom whose derivation is wrong.
     */
    public abstract void test(TotalEnergy totalEnergy, Atom atom) throws EvaluationException;

    public boolean evaluatesResidues() {
        return false;
    }


    //----------------------------------------- housekeeping ---------------------------------------------

    public void handleMissingParameters(Object obj) {
        throw new RuntimeException("Missing parameters for:\n" + obj);
    }

    public String comment() {
        return comment;
    }

    public String toString() {
        return comment;
    }


    //--------------------------------------- auxiliary methods--------------------------------------------------

    /**
     * Generates an empty array.
     */
    protected static Updateable[] toArray() {
        Updateable[] out = {};
        return out;
    }

    /**
     * Generates an aray with one element - the parameter.
     */
    protected static Updateable[] toArray(Updateable o1) {
        Updateable[] out = {o1};
        return out;
    }

    protected static Updateable[] toArray(Updateable o1, Updateable o2) {
        Updateable[] out = {o1, o2};
        return out;
    }

    protected static Updateable[] toArray(Updateable o1, Updateable o2, Updateable o3) {
        Updateable[] out = {o1, o2, o3};
        return out;
    }
    // ------------------------------------- internal auxiliary classes -----------------------------------------

    /**
     * A list of updateable elements (implementing the meshi.util.Updateable interface).
     */
    protected class UpdateableList extends ArrayList<Updateable> {
        public UpdateableList(Updateable[] array) {
            super();
            if (array.length == 0) return;
            for (Updateable updatable : array)
                add(updatable);
        }

        public UpdateableList() {
            super();
        }

        /**
         * Iterates over the list updates each of its elements.
         */
        public void update(int numberOfUpdates) throws UpdateableException {
            for (Updateable resource : this) {
                resource.update(numberOfUpdates);
            }
        }
    }


    private static class IsEnergy implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof AbstractEnergy);
        }
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        if (updateableResources == null) return;
        if (updateableResources.size() == 0) return;
        for (Updateable updateable : updateableResources) {
            if (updateable == null)
                throw new RuntimeException("Null element in updateableResources " + updateableResources);
            updateable.setNumberOfUpdates(numberOfUpdates);
        }
    }

    protected EnergyInfoElement energy() {
        return info;
    }


    //----------------------------------- Attributes ---------------------------------------     
    private HashMap attributes = new HashMap();

    public final void addAttribute(MeshiAttribute attribute) {
        attributes.put(attribute.key(), attribute);
    }

    public final MeshiAttribute getAttribute(int key) {
        return (MeshiAttribute) attributes.get(key);
    }
}


  
