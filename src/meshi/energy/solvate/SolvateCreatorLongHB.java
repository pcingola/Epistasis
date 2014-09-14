package meshi.energy.solvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The different between this creator and the "SolvateCreatorRegularHB" (see the latter main comment for a lot of details), is in the 
 * way the hydrogen-bond distance sigmoid is set. Here parameters are chosen so that the hydrogen-bond drops to 0 for donor-acceptor 
 * distance of 3.6 Ang. This is good in cases where the atom placement is approximate, e.g. when using discreet rotamers.  
 **/

public class SolvateCreatorLongHB extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;
    private static SolvateParametersList parametersList = null;

    public SolvateCreatorLongHB(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
    					double weightBBPolarSolvate, double weightHB) {
		super(1.0);
		this.weightSCPolarSolvate = weightSCPolarSolvate;
		this.weightSCCarbonSolvate = weightSCCarbonSolvate;
		this.weightBBPolarSolvate = weightBBPolarSolvate;
		this.weightHB = weightHB;
    }
   
    public SolvateCreatorLongHB(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorLongHB() {
		super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_LONG_HB_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_LONG_HB_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_LONG_HB_PARAMETERS[cc]);
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
