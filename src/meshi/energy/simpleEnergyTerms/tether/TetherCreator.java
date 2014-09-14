package meshi.energy.simpleEnergyTerms.tether;
import meshi.energy.simpleEnergyTerms.*;

import meshi.energy.EnergyCreator;

import java.util.*;
import meshi.util.*;
import meshi.util.filters.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.*;

import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: guykarl
 * Date: 23/06/2005
 * Time: 12:32:27
 * To change this template use File | Settings | File Templates.
 */

import meshi.energy.*;

import meshi.molecularElements.*;


public class TetherCreator extends SimpleEnergyCreator  implements KeyWords {
    Filter filter = null;
    String comment = null;
    public TetherCreator() {
	super( TETHER_ENERGY);
    }

    public TetherCreator(Key key, Filter filter, CommandList commands,String comment) {
	super(key,commands);
	this.filter = filter;
	this.comment = comment;
    }

    public TetherCreator(Key key) {
	super(key);
	
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
					   CommandList commands) {
 	AtomList atomList = protein.atoms();
	double[][] initialLocationsMatrix=new double[3][atomList.size()];
	for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
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
	
	term = new  TetherEnergy(atomList, initialLocationsMatrix,weight(),comment);
	return term;
    }
}
