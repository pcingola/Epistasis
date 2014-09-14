package meshi.energy.simpleEnergyTerms.inflate.simpleInflate;
import  meshi.energy.simpleEnergyTerms.inflate.*;
import meshi.energy.simpleEnergyTerms.inflate.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import java.util.*;
import meshi.optimizers.*;

public class SimpleInflate extends SimpleEnergyTerm{
    private static ArrayList dummy = null;
    //private static Parameters parameters = new InflateParameters();
    private static AtomList atomList, atomListCopy;
    private DistanceMatrix distanceMatrix;
    private Filter filter;
    private double rmsTarget; 

    public SimpleInflate() {}

    public SimpleInflate(DistanceMatrix distanceMatrix, double rmsTarget, double weight) {
	super(toArray(distanceMatrix),null, weight);
	comment = "Inflate ;)";
	this.distanceMatrix = distanceMatrix;
	off();
	this.rmsTarget = rmsTarget;
    }
    public double evaluate() {
	if (! on) return 0;
	boolean targetReached;
	try {
	    targetReached = atomList.getRms(atomListCopy) > rmsTarget;
	}
	catch (Exception ex) {
	    System.out.println("Simple inflate failed due to:\n"+ex+"\n"+"Trying to continue.");
	    targetReached = true;
	}
   if (targetReached) {
	    System.out.println("inflate reached RMS of "+rmsTarget);
	    Minimizer.terminator.kill("Inflate reached RMS of "+rmsTarget);
	}
	return super.evaluate();
    }

    public void on() {
	atomList = new AtomList();
	for (AtomCore atom:distanceMatrix.molecularSystem)
	    atomList.add(atom.atom);
	atomListCopy = Utils.duplicateInAnewMolecularSystem(atomList);
       	createElementsList(distanceMatrix.nonBondedList());
	super.on();
    }


    public void createElementsList(DistanceList baseList) {
	elementsList = new ArrayList();
	
	Parameters parameters;
	for (Distance baseElement:baseList) {
	    EnergyElement newElement = createElement(baseElement, null);
	    if (! newElement.frozen())
		elementsList.add(newElement);
	}
    }
    
    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new SimpleInflateEnergyElement(((Distance)baseElement), distanceMatrix, weight);
    }

    public static class TargetFilter implements Filter {
	private double target;
	public TargetFilter(double target) {
	    this.target = target;
	}
	public boolean accept(Object obj) {
	    AtomPair ap = (AtomPair) obj;
	    if (ap.atom1().distanceFrom(ap.atom2()) < target) return true;
	    return false;
	}
    }
    
    public double weight() {return weight;}

}

	
    
