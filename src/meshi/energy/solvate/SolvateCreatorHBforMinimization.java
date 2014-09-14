package meshi.energy.solvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The different between this creator and the "SolvateCreatorRegularHB" (see the latter main comment for a lot of details), is in the 
 * way the hydrogen-bond distance sigmoid is set. Here parameters are chosen so that the hydrogen-bond starts its drop to 0 at a 
 * donor-acceptor distance of 3.1 Ang, and is 0 only for donor-acceptor distance of 4.0 Ang. This is good to increase the basin of 
 * attraction of direct minimizations protocols.   
 **/

public class SolvateCreatorHBforMinimization extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;
    private static SolvateParametersList parametersList = null;
    private String comment = null;
    public SolvateCreatorHBforMinimization(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
					   double weightBBPolarSolvate, double weightHB, CommandList commands) {
	this( weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate, weightHB, commands, "Solvation");
    }

    public SolvateCreatorHBforMinimization(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
					   double weightBBPolarSolvate, double weightHB, CommandList commands, String comment) {
	super(SOLVATE_ENERGY,commands);
	this.weightSCPolarSolvate  = weightSCPolarSolvate  * weight();
	this.weightSCCarbonSolvate = weightSCCarbonSolvate * weight();
	this.weightBBPolarSolvate  = weightBBPolarSolvate  * weight();
	this.weightHB              = weightHB              * weight();
	this.comment               = comment;
    }
   
    public SolvateCreatorHBforMinimization(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorHBforMinimization() {
		super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_MINIMIZE_HB_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_MINIMIZE_HB_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_MINIMIZE_HB_PARAMETERS[cc]);
	    	parametersList = new SolvateParametersList(strlist);
	 	}

		if (comment == null) return new SolvateEnergy(protein.atoms(), 
							      distanceMatrix, 
							      (SolvateParametersList) parametersList,
							      weightSCPolarSolvate,
							      weightBBPolarSolvate,
							      weightSCCarbonSolvate,
							      weightHB);
		else  return new SolvateEnergy(protein.atoms(), 
							      distanceMatrix, 
							      (SolvateParametersList) parametersList,
							      weightSCPolarSolvate,
							      weightBBPolarSolvate,
							      weightSCCarbonSolvate,
					       weightHB,comment);
    }

}
