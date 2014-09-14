/**
 *Excluded Volume potential.
 *
 *A potential that creates a strong repulsion between atoms when they are getting
 *near their VDW radius. There is no attraction part like in the VDW.
 *The functional form of the term is:  
 *
 *dis =                EV
 *
 *[0,sigma]  C*(dis-sigma)^4
 *[sigma,Inf]          0
 *
 *ALPHA is set in ExcludedVolParameters.java. Currently it is 0.2 Ang.
 *ALPHA is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
 *
 **/
package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceList;
import meshi.geometry.Distance;
import meshi.util.Classes;
import meshi.util.Utils;
import meshi.util.filters.Filter;

public class ExcludedVol extends NonBondedEnergyTerm implements Classes{
    Filter filter;

    public ExcludedVol(){super();}

    public ExcludedVol(DistanceMatrix distanceMatrix,
                           ExcludedVolParametersList parametersList,
                           int type,
                           double weight,
                          double Rfac,
                           Filter filter){
	super(toArray(distanceMatrix), weight, distanceMatrix);
	comment = "ExcludedVol";
	energyElement = new ExcludedVolEnergyElement(parametersList, distanceMatrix, type, weight,Rfac);
    this.filter = filter;    
    }


    public double evaluate() {
	if (! on) return 0.0;
	double energy = 0;
	DistanceList nonBondedList;
	if (filter == null)
		nonBondedList = distanceMatrix.nonBondedList();
	else {
	    nonBondedList = new DistanceList(100);
	    Utils.filter(distanceMatrix.nonBondedList(),filter,nonBondedList);
	}
      if(nonBondedList  == null)
                System.out.println("!!!!!!");

    for (Distance distance:nonBondedList) {
            if (distance != null &&  !distance.mode().frozen) {
                energyElement.set(distance);
	            energy += energyElement.evaluate();
            }
        }
	return energy;
    }//evaluate
}

	
 
