package meshi.energy.simpleEnergyTerms.alphaTorsion;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import java.util.*;

/**
 *This energy limits the alpha torsion (4 consecutive CAs) to be in a specific range.
 *It is also secondary structure sensitive, and the range depends on the residue SS.
 *It operates currently only on HELIX,SHEET secondary structure states, since the COIL,ALL 
 *states practically don't have any limitation on the alpha torsion.
 *
 *Important Note: This energy term must be accompanied by an ALPHA-angle energy term. This is 
 *because it has a non-continous point at torsion values of -Pi or Pi , and also when 
 *one of the 2 angles that make up the torsion is close to 0 or Pi. These discontinuites should 
 *not affect normal operation if the ALPHA-angle term is working. On very rare starting condition
 *these problems might never the less be encountered.  
 **/
 
             
public class AlphaTorsionEnergy extends SimpleEnergyTerm{
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public AlphaTorsionEnergy() {}

    public AlphaTorsionEnergy(TorsionList torsionList, DistanceMatrix distanceMatrix, 
		       AlphaTorsionParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	this.torsionList = torsionList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionList);
	comment = "alphaTorsion";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new AlphaTorsionEnergyElement((Torsion)baseElement,
	                                    (AlphaTorsionParameters) parameters, weight);
    }

    public void handleMissingParameters(Object obj) {}

    public void removeBadBonds(double farAway) {
	for (Iterator elements = elementsList.iterator(); elements.hasNext();) {
	    AlphaTorsionEnergyElement element = (AlphaTorsionEnergyElement) elements.next();
	    Atom                       atom1   = element.atom1();
	    Atom                       atom2   = element.atom2();
	    Atom                       atom3   = element.atom3();
	    Atom                       atom4   = element.atom4();
	    double                     dis12   = atom1.distanceFrom(atom2);
	    double                     dis23   = atom2.distanceFrom(atom3);
	    double                     dis34   = atom3.distanceFrom(atom4);
	    if ((dis12>farAway) ||(dis23>farAway) ||(dis34>farAway))  
		element.turnOff();
	    else element.turnOn();
	}
    }

}    
