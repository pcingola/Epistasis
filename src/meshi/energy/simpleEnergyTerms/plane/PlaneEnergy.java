package meshi.energy.simpleEnergyTerms.plane;
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
 * Plane energy term. 
 **/
public class PlaneEnergy extends SimpleEnergyTerm{
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public PlaneEnergy() {}

    public PlaneEnergy(TorsionList torsionList, DistanceMatrix distanceMatrix, 
		       PlaneParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	this.torsionList = torsionList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionList);
	comment = "Plane";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new PlaneEnergyElement(((Torsion)baseElement), parameters, weight);
    }
}    
