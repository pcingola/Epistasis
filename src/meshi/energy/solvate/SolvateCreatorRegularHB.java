package meshi.energy.solvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The Creator class for setting up an instance of SolvateEnergy as described in Kalisman & Keasar (2008).
 * Unlike the paper it is possible to give a different weight to the solvation energies of certain atom types in the final
 * summation described in Eq. 5. These atom types are side-chain polars, side-chain carbons, and backbone 
 * polars. The class also include a regular hydrogen bond term. The functional form is therefore: 
 *
 * Esolv = weightSCPolarSolvate*Eside_chain_polars + 
 * 		weightSCCarbonSolvate*Eside_chain_carbons + 
 *      weightBBPolarSolvate*Ebackbone_polars + 
 *      weightHB*Ehb
 *      
 * Where Ehb is the negative of the HBC summation over all atoms. The weights are given as parameters to the constructor. 
 * Default values of 1.0 for the weights are as parameters to the constructors with the shorter signatures. 
 *
 * CREATOR TYPES
 * -------------
 * This creator (SolvateCreatorRegularHB), define hydrogen bonds in a strict manner, as defined by the hydrogen bond atlas
 * of McDonald and Thoronton (1994). The values of the hydrogen-bond distance sigmoid are chosen so that the hydrogen-bond 
 * is nearly 0 in value for donor-acceptor distance above 3.2 Ang. 
 *
 **/

public class SolvateCreatorRegularHB extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;
    private static SolvateParametersList parametersList = null;

    public SolvateCreatorRegularHB(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
    					double weightBBPolarSolvate, double weightHB) {
		super(1.0);
		this.weightSCPolarSolvate = weightSCPolarSolvate;
		this.weightSCCarbonSolvate = weightSCCarbonSolvate;
		this.weightBBPolarSolvate = weightBBPolarSolvate;
		this.weightHB = weightHB;
    }
   
    public SolvateCreatorRegularHB(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorRegularHB() {
		super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_PARAMETERS[cc]);
	    	parametersList = new SolvateParametersList(strlist);
	 	}

		return new SolvateEnergy(protein.atoms(), 
					distanceMatrix, 
					(SolvateParametersList) parametersList,
				    weightSCPolarSolvate,
				    weightBBPolarSolvate,
				    weightSCCarbonSolvate,
				    weightHB);
    }

}
