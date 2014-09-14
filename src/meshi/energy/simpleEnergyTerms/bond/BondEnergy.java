package meshi.energy.simpleEnergyTerms.bond;
import meshi.util.*;
import meshi.molecularElements.*;   
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import java.util.*;
/**
 * Bond energy term. 
 * Has the general form <b> Eb = SIGMAi(Ki(Di-D0i)^2)</b>
 * where <b>Di</b> 
 * is the distance between the two bonded atoms, <b> D0i </b> 
 * is their expected average distance (depends on their types) and <b>Ki</b>
 * is a force constant that again, depends on the atom types.<br>
 * This class is used for both calculating the bond-energy term of an energy function 
 * and for updating the forces on each atom accordingly.<b>
 * It is assumed that the list of bonds is constant during the simulation. That is 
 * no bonds are made or broken.
 */
public class BondEnergy extends SimpleEnergyTerm {
    /**
     * The constructor associates any bond with its parameters.
     **/
    protected DistanceMatrix distanceMatrix;

    public BondEnergy() {}

    public BondEnergy(AtomPairList bondList, 
		      DistanceMatrix distanceMatrix,
		      BondParametersList  parametersList, 
		      double weight) {
	super(toArray(distanceMatrix), parametersList, weight);
	comment = "Bond";
	this.distanceMatrix = distanceMatrix;
	createElementsList(bondList);
    }

 
    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	AtomPair atomPair = (AtomPair)baseElement;
	if (atomPair.atom1().nowhere() || atomPair.atom2().nowhere()) return null; 
	if ((distanceMatrix.distance(atomPair.atom1(),atomPair.atom2()) == null) &&
	    (distanceMatrix.distance(atomPair.atom2(),atomPair.atom1()) == null)) return null; 
	return new BondEnergyElement(atomPair, parameters, distanceMatrix, weight);
    }
    
    public void removeBadBonds(double farAway) {
	for (Iterator elements = elementsList.iterator(); elements.hasNext();) {
	    BondEnergyElement element = (BondEnergyElement) elements.next();
	    Atom              atom1   = element.atom1();
	    Atom              atom2   = element.atom2();
	    if (atom1.distanceFrom(atom2) > farAway) 
		element.turnOff();
	    else element.turnOn();
	}
    }
}
	
	
