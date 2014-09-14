package meshi.energy.simpleEnergyTerms;
import meshi.energy.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.util.filters.*;
import meshi.util.file.*;
import meshi.molecularElements.*;
import meshi.geometry.*;
import java.util.*;
/** 
 * A super class for those energy terms that operate on a fixed list of elements. 
 * Typical cases are the bonded terms (bond angle etc.), where the list of elements 
 * (say covalent bonds)
 * is fixed by the chemical structure of the molecule. 
 **/
public abstract class SimpleEnergyTerm extends AbstractEnergy  {

    protected ArrayList<EnergyElement> elementsList;

    /**
     * A list of all the parameters needed by the energy term.
     **/
    protected ParametersList parametersList;

    public SimpleEnergyTerm() {}

    public SimpleEnergyTerm(Updateable[] updateableResources,
			    ParametersList  parametersList, 
			    double weight) {
	super(updateableResources, weight);
	this.parametersList = parametersList;
    }

    /**
     * Creates the fixed list of energy elements in the initiation phase of a simulation. 
     * These elements typically warp some object such as AtomPair, Angle or Torsion.
     * The warped elements are provided by the parameter "baseList". 
     * Note that Different energy terms treat this list differently. 
     * BondEnergy for example tries to warp every single element of "baseList" while PlaneEnergy 
     * assumes that the list is redundant and warps only those Torsion elements for which it finds parameters.
     * The binding of "baseElement" to the parameters is done by the method parameters of ParametersList.
     **/    
    public void createElementsList(ArrayList baseList) {
	elementsList = new ArrayList();

	Parameters parameters;
	for(Object baseElement:baseList) {
	    parameters = parametersList.parameters(baseElement);
	    if (parameters == null) handleMissingParameters(baseElement);
	    else {
		EnergyElement newElement = createElement(baseElement, parameters);
		if ((newElement != null) &&(! newElement.frozen()))
		    elementsList.add(newElement);
	    }
	}
    }


    public ArrayList elementsList() {return elementsList;} 
    
    public abstract EnergyElement createElement(Object baseElement, Parameters parameters);

    /**
     * Testing of one atom in all energy elements
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom an criminal <code>Atom</code> value
     */
    public void test(TotalEnergy totalEnergy,Atom atom){
        if (! on) System.out.println(""+this +" is off");
        for(Iterator eei = elementsList.iterator();eei.hasNext();){
            EnergyElement energyElement = (EnergyElement)eei.next();
	    if(!energyElement.frozen()){
                energyElement.test(totalEnergy,atom);
	    }
	}
    }
 
   /**
    * Evaluates the current value of the energy function and <b> update </b> the derivatives.
    **/ 
    public double evaluate() {
	if (! on) return 0.0;
	double e, energy = 0;
	for (EnergyElement energyElement:elementsList) {
	    e = energyElement.evaluate();
	    energy += e;
	}
	return energy;
    }
    
    public void evaluateAtoms() {
	if (on) {
		for (EnergyElement energyElement:elementsList) {
		energyElement.evaluateAtoms();
	    }
	}
    }
    

}


  
