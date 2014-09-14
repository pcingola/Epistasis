package meshi.energy.simpleEnergyTerms.inflate.simpleInflate;
import  meshi.energy.simpleEnergyTerms.inflate.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.*;

public class SimpleInflateCreator extends EnergyCreator  implements KeyWords {
    public SimpleInflateCreator() {
	super(INFLATE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	Command command = commands.firstWordFilter(key).secondWord(RMS_TARGET);
	double rmsTarget = command.thirdWordDouble();
	term = new SimpleInflate(distanceMatrix, rmsTarget, weight());
	return term;
    }
}
