/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.tether;

import meshi.energy.simpleEnergyTerms.*;

/**
 * Created by IntelliJ IDEA.
 * User: guykarl
 * Date: 23/06/2005
 * Time: 12:22:59
 * To change this template use File | Settings | File Templates.
 */


import meshi.molecularElements.atoms.*;

import meshi.energy.*;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;

import java.util.*;

public class TetherEnergy extends SimpleEnergyTerm {
    /**
     * The constructor associates any bond with its parameters.
     */
    protected double[][] initialLocationsMatrix;

    public TetherEnergy() {
    }

    public TetherEnergy(AtomList atomList, double[][] inputMatrix, EnergyInfoElement info, String comment) {
        super(toArray(), null, info);
        initialLocationsMatrix = inputMatrix;
        elementsList = new ArrayList();
        this.weight = weight;
        for (int i = 0; i < atomList.size(); i++) {
            //for (Iterator baseElements = atomList.iterator(); baseElements.hasNext();) {
            Atom atom = atomList.get(i);
            if (!atom.nowhere()) {
                EnergyElement newElement = createElement(atom);
                if (!newElement.frozen()) {
                    elementsList.add(newElement);
                }
            }
        }
        if (comment == null) this.comment = "Tether";
        else this.comment = comment;
    }


    public boolean evaluatesResidues() {return false;}

    public void scaleWeight(double factor) {
        super.scaleWeight(factor);
        for (EnergyElement element : elementsList)
            ((TetherEnergyElement) element).force *= factor;
    }


    public void  restoreWeight() {
             for (EnergyElement element : elementsList)
                ((TetherEnergyElement) element).restoreWeight();
    }
    public EnergyElement createElement(Object baseElement) {
        Atom atom = (Atom) baseElement;
        EnergyElement out = new TetherEnergyElement(atom, initialLocationsMatrix[0][atom.number()],
                initialLocationsMatrix[1][atom.number()], initialLocationsMatrix[2][atom.number()], weight);
        return out;
    }

    public boolean allAtomsTether() {
        for (Object element : elementsList) {
            if (((TetherEnergyElement) element).reliability() != 1) return false;
        }
        return true;
    }

    public void setResidue(CommandList commands) {
        for (EnergyElement element : elementsList) {
            ((TetherEnergyElement) element).setReliability(0);
        }
        if (commands.keyExists("residueTether")) {
            String line = commands.firstWord("residueTether").secondWord();
            if (line.startsWith("allCa")) {
                for (EnergyElement element : elementsList) {
                    TetherEnergyElement tetherElement = (TetherEnergyElement) element;
                    if (tetherElement.atom.name().equals("CA")) {
                        tetherElement.setReliability(1.0);
                        if (Utils.verbose()) System.out.println("Tethering " + tetherElement.atom);
                    }
                }
            }
            else {
                String[] words = line.split(",");
                for (String word : words) {
                    int residueNumber = Integer.valueOf(word);
                    for (EnergyElement element : elementsList) {
                        TetherEnergyElement tetherElement = (TetherEnergyElement) element;
                        if (tetherElement.atom.residueNumber() == residueNumber) {
                            tetherElement.setReliability(1.0);
                            if (Utils.verbose()) System.out.println("Tethering " + tetherElement.atom);
                        }
                    }
                }
            }
        } else {
            on = false;
            if (Utils.verbose()) System.out.println("No residue was specifically tethered");
        }
    }
    public void setBackbone(CommandList commands) {
           for (EnergyElement element : elementsList) {
               ((TetherEnergyElement) element).setReliability(0);
           }
           if (commands.keyExists("tetherBackboneOn")) {
                   for (EnergyElement element : elementsList) {
                       TetherEnergyElement tetherElement = (TetherEnergyElement) element;
                       if (tetherElement.atom.isBackbone()) {
                           tetherElement.setReliability(1.0);
                           if (Utils.verbose()) System.out.println("Tethering " + tetherElement.atom);
                       }
                   }
           } else {
               on = false;
               if (Utils.verbose()) System.out.println("No residue was specifically tethered");
           }
       }


    public void update() {
    }

    public void update(int i) {
    }

    public EnergyElement createElement(Object baseElement, Parameters p) {
        return null;
    }

    public void reset() {
        for (EnergyElement e : elementsList) {
            TetherEnergyElement element = (TetherEnergyElement) e;
            element.reset();
        }
        on();
    }

    public void remove(Atom atom) {
        for (EnergyElement e : elementsList) {
            TetherEnergyElement element = (TetherEnergyElement) e;
            if (element.atom == atom) {
                element.remove();
            }
        }
    }

    /*public void restart(){

            if (atomList == null)
                throw new RuntimeException ("Can not make TetherCreator.restart() while TetherCreator.atomList is empty");
            if (initialLocationsMatrix == null)
                throw new RuntimeException ("Can not make TetherCreator.restart() while TetherCreator.initialLocationsMatrix is empty");

            for (Iterator getAtoms = atomList.iterator(); getAtoms.hasNext();) {
                Atom atom = (Atom) getAtoms.next();
                if (! atom.nowhere()) {
                int number = atom.number();
                try {
                initialLocationsMatrix[0][number]=atom.x();
                initialLocationsMatrix[1][number]=atom.y();
                initialLocationsMatrix[2][number]=atom.z();
                }
                catch (RuntimeException ex) {
                    System.out.println("Problem in TetherCreator. createEnergyTerm while processing atom number :"+number+"\n"+atom);
                    throw ex;
                }
                    if (filter != null) {
                        if (filter.accept(atom)) atom.setReliability(1);
                        else atom.setReliability(0);
                    }
                    else atom.setReliability(1);
                }
            }
            if (term != null)
                   ((TetherEnergy)term).restart(atomList,initialLocationsMatrix);
        }

    public void restart(AtomList atomList, double[][] inputMatrix){
            initialLocationsMatrix=inputMatrix;
            elementsList.clearRows();

            for (Iterator baseElements = atomList.iterator(); baseElements.hasNext();) {
                Atom atom = (Atom) baseElements.next();
                if ((! atom.nowhere()) & (atom.reliability() > 0))  {
                EnergyElement newElement = createElement(atom);
                if (! newElement.frozen()) {
                    elementsList.add(newElement);
                }
                }
            }
        }                                            */

}


