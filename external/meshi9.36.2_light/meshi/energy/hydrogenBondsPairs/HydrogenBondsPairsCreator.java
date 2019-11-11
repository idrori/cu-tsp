/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

/*
 * Created reset 27/12/2004
 */
package meshi.energy.hydrogenBondsPairs;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * @author amilev
 */
public class HydrogenBondsPairsCreator extends EnergyCreator {

    //----------------------------------- data ----------------------------
    
	HydrogenBondsCreator hydrogenBondsCreator;
    private BetaParametersList betaParametersList;
    private HelixParametersList helixParametersList;
    private BetaParametersList betaBackgroundParametersList;
    private HelixParametersList helixBackgroundParametersList;
    private HydrogenBondsEnergy HBE = null; //chen 26.9.11
    private PairsOfHBEElementsList elements = null;
    //---------------------------------- constructors ---------------------
    
	/**
	 * @param hydrogenBondsCreator used to extract the HB list
	 **/
	public HydrogenBondsPairsCreator(HydrogenBondsCreator hydrogenBondsCreator) {
		super(InfoType.HYDROGEN_BONDS_PAIRS);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}

	/**
	 * @param weight of the energy term
	 * @param hydrogenBondsCreator used to extract the HB list
	 **/
	public HydrogenBondsPairsCreator(double weight,HydrogenBondsCreator hydrogenBondsCreator) {
		super(weight,InfoType.HYDROGEN_BONDS_PAIRS);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}

    //---------------------------------- methods ---------------------------

    /*
    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands,HydrogenBondsEnergy HBE,double punish,double hpunish) {
        if (parametersList== null)
			{
				parametersList = new HydrogenBondsPairsParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_PARAMETERS);
			}
        HydrogenBondsParametersList HBparamsList = HBE.getParametersList();
        return new HydrogenBondsPairsEnergy (distanceMatrix,
                                             ( HydrogenBondsPairsParametersList) parametersList,
                                             HBparamsList,
                                             new PairsOfHBEElementsList(HBE),
                                             weight(),
                                             punish,
                                             hpunish);
	}
	
    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands,HydrogenBondsEnergy HBE,double punish) {
        if (parametersList== null)
			{
				parametersList = new HydrogenBondsPairsParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_PARAMETERS);
			}
        HydrogenBondsParametersList HBparamsList = HBE.getParametersList();
        return new HydrogenBondsPairsEnergy (distanceMatrix,
                                             ( HydrogenBondsPairsParametersList) parametersList,
                                             HBparamsList,
                                             new PairsOfHBEElementsList(HBE),
                                             weight(),
                                             punish,
                                             -50);
	}
	
    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands,HydrogenBondsEnergy HBE) {
        if (parametersList== null)
			{
				parametersList = new HydrogenBondsPairsParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_PARAMETERS);
			}
        HydrogenBondsParametersList HBparamsList = HBE.getParametersList();
        return new HydrogenBondsPairsEnergy (distanceMatrix,
                                             ( HydrogenBondsPairsParametersList) parametersList,
                                             HBparamsList,
                                             new PairsOfHBEElementsList(HBE),
                                             weight(),
                                             -10,
                                             -50);
	}
    */
    
	/** 
	 * 
	 **/
	public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands) {

        //get the HydrogenBondsEnergy  which will be used to  instruct the HB list (than generate the pairs list based reset it)
        HBE = hydrogenBondsCreator.getHydrogenBondsEnergy();
        elements = new PairsOfHBEElementsList(HBE,distanceMatrix);
        if (helixParametersList==null |  betaParametersList== null | betaBackgroundParametersList == null | helixBackgroundParametersList == null)
			{
                helixParametersList = new HelixParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_HELIX_PARAMETERS);
                betaParametersList = new BetaParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_BETA_PARAMETERS);
                helixBackgroundParametersList = new HelixParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_HELIX_BACKGROUND_PARAMETERS);
                betaBackgroundParametersList = new BetaParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_BETA_BACKGROUND_PARAMETERS);


            }
        double width = commands.firstWord(HYDROGEN_BONDS).thirdWordDouble();
               EnergyInfoElement energyInfoElement = new EnergyInfoElement(infoType, "Cooperative HB term", weight);
        return term = new HydrogenBondsPairsEnergy (distanceMatrix,
                                             helixParametersList,
                                             betaParametersList,
                                             helixBackgroundParametersList,
                                             betaBackgroundParametersList,
                                             elements,
                                             energyInfoElement,
                                             hydrogenBondsCreator .getSpecialDisArray(),
                                             hydrogenBondsCreator .getAntiParalel(),width );
	}

}
