package meshi.symmetryComplex.energy.edmEnergy;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.ParametersList;


public class EDMEnergyParametersList extends ParametersList {

	private EDMEnergyParameters ep = null;
	
	public EDMEnergyParametersList( String edmDataFileName ) {
		super();
		
		EDMSplinedPolynomial esp = null;
		
		try {
			esp = new EDMSplinedPolynomial( edmDataFileName );
		}
		catch( Exception e ) {
			e.printStackTrace();
			throw new RuntimeException( "unable to read polynomials parameters file" );
		}
		
		/* parameters are always the same, prepare them in advance */
		ep = new EDMEnergyParameters( esp );
	}

	/** the only relevant "parameters" are the splined polynomial,
	 * so always return it through a dummy Parameters class.
	 */
	public Parameters parameters(Object Obj) {
		return ep;
	}

	public Parameters createParameters(String line) {
		throw new RuntimeException( "EDMEnergy term uses " +
			"EDMSplinedPolynomial in order to create parameters list" );
	}

}
