/*
 * Created on 27/12/2004
 */
package meshi.energy.hydrogenBondsPairs;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.*;
import meshi.util.CommandList;

/**
 * @author amilev
 */
public class HydrogenBondsPairsCreator extends EnergyCreator {

    //----------------------------------- data ----------------------------
    
	HydrogenBondsCreator hydrogenBondsCreator;
    private BetaParametersList betaParametersList;
    private HelixParametersList helixParametersList;

    //---------------------------------- constructors ---------------------
    
	/**
	 * @param hydrogenBondsCreator
	 **/
	public HydrogenBondsPairsCreator(HydrogenBondsCreator hydrogenBondsCreator) {
		super(HYDROGEN_BONDS_PAIRS);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}

	/**
	 * @param weight
	 * @param hydrogenBondsCreator
	 **/
	public HydrogenBondsPairsCreator(double weight,HydrogenBondsCreator hydrogenBondsCreator) {
		super(weight);
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
        //get an HydrogenBondsEnergy element witch will bw used to generate the pairs list from its HB list
        HydrogenBondsEnergy HBE = hydrogenBondsCreator.getHydrogenBondsEnergy();//(HydrogenBondsEnergy)(new HydrogenBondsCreator(0).createEnergyTerm(protein,distanceMatrix,commands));
        if (helixParametersList==null |  betaParametersList== null)
			{
                helixParametersList = new HelixParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_HELIX_PARAMETERS);
                betaParametersList = new BetaParametersList(parametersDirectory(commands)+
                                                                      "/"+HYDROGEN_BONDS_PAIRS_BETA_PARAMETERS);
			}
        return new HydrogenBondsPairsEnergy (distanceMatrix,
                                             helixParametersList,
                                             betaParametersList,
                                             new PairsOfHBEElementsList(HBE),
                                             weight(),
                                             hydrogenBondsCreator .getSpecialDisArray(),
                                            hydrogenBondsCreator .getAntiParalel() );
	}

}
