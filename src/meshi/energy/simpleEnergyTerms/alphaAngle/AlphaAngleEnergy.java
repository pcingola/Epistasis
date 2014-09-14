package meshi.energy.simpleEnergyTerms.alphaAngle;
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
 *This class limits the alpha angle (3 consecutive CAs) to be in a specific range.
 *It is also secondary structure sensitive, and the range depends on the residue SS.
 *The range limitiation is done by putting very steep parabolas on both side of the 
 *range. Within the range the energy is a constant 0.0
 *
 *Important Note: This energy term has a non-continous point at angle values of 0 or Pi.
 *These discontinuites should not normally affect normal operation if the weight set to this
 *term is sufficiently high (of the order of 1-10).On very rare starting condition, however,
 *these problems might never the less be encountered.  
 **/
 
public class AlphaAngleEnergy extends SimpleEnergyTerm{
    protected AngleList angleList;
    protected DistanceMatrix distanceMatrix;

    public AlphaAngleEnergy() {}

    public AlphaAngleEnergy(AngleList angleList, DistanceMatrix distanceMatrix, 
		       AlphaAngleParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, angleList), parametersList, weight);
	this.angleList = angleList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(angleList);
	comment = "alphaAngle";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new AlphaAngleEnergyElement((Angle)baseElement,
	                                    (AlphaAngleParameters) parameters, weight);
    }

    public void handleMissingParameters(Object obj) {}

    public void removeBadBonds(double farAway) {
	for (Iterator elements = elementsList.iterator(); elements.hasNext();) {
	    AlphaAngleEnergyElement element = (AlphaAngleEnergyElement) elements.next();
	    Atom                    atom1   = element.atom1();
	    Atom                    atom2   = element.atom2();
	    Atom                    atom3   = element.atom3();
	    double                  dis12   = atom1.distanceFrom(atom2);
	    double                  dis23   = atom2.distanceFrom(atom3);
	    if ((dis12>farAway) ||(dis23>farAway))  
		element.turnOff();
	    else element.turnOn();
	}
    }

}    
