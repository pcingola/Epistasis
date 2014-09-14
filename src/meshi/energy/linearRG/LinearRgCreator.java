package meshi.energy.linearRG;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.string.*;
import java.util.*;
import meshi.util.*;

/**
 *An implicit solvation energy term for all-atom models, modeling a 4.0 angs solvation shell around
 *each atom.
 **/

public class LinearRgCreator extends EnergyCreator  implements KeyWords {
    double targetRG = 0;

    public LinearRgCreator() {
	super(LINEAR_RG);
    }
    public LinearRgCreator(double targetRG) {
	super(LINEAR_RG);
	this.targetRG = targetRG;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		return new LinearRgEnergy(protein.atoms(), weight(), targetRG, new RgCalculator(protein.atoms()));
    }

}
