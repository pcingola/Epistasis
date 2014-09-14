package meshi.energy.simpleEnergyTerms.tether;
import meshi.energy.simpleEnergyTerms.*;

/**
 * Created by IntelliJ IDEA.
 * User: guykarl
 * Date: 23/06/2005
 * Time: 12:22:59
 * To change this template use File | Settings | File Templates.
 */


import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

import meshi.energy.*;

import java.util.*;

public class TetherEnergy  extends SimpleEnergyTerm{
    /**
     * The constructor associates any bond with its parameters.
     **/
    protected double[][] initialLocationsMatrix;
    private static int counter;
    public TetherEnergy() {}
    public TetherEnergy(AtomList atomList, double[][] inputMatrix,double weight,String comment) {
	initialLocationsMatrix=inputMatrix;
	elementsList = new ArrayList();
	this.weight = weight;
	Object baseElement;

	for (Iterator baseElements = atomList.iterator(); baseElements.hasNext();) {
	    Atom atom = (Atom) baseElements.next();
	    if ((! atom.nowhere()) & (atom.reliability() > 0))  {
		EnergyElement newElement = createElement(atom);
		if (! newElement.frozen()) {
		    elementsList.add(newElement);
		}
	    }
        }
	if (comment == null) this.comment = "Tether";
	else this.comment = comment;
    }


    public EnergyElement createElement(Object baseElement) {
	Atom atom = (Atom) baseElement;
	EnergyElement out = new TetherEnergyElement(atom,initialLocationsMatrix[0][atom.number()],
				       initialLocationsMatrix[1][atom.number()] ,initialLocationsMatrix[2][atom.number()], weight);
	return out;
    }

     public void update(){}
     public void update(int i ){}

        public EnergyElement createElement(Object baseElement,Parameters p){return null;}

}


