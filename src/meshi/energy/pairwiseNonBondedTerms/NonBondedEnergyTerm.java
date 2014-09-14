package meshi.energy.pairwiseNonBondedTerms;
import meshi.energy.*;
import meshi.util.*;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;

import java.util.Iterator;

public abstract class NonBondedEnergyTerm extends AbstractEnergy{
    protected DistanceMatrix distanceMatrix;
    /**
     * Must be initialized on subclass construction.
     **/
    protected NonBondedEnergyElement energyElement=null;

    public NonBondedEnergyTerm(){super();}

    /**
     * Creates a new <code>NonBondedEnergyTerm</code> instance.
     *
     * @param updateableResources an <code>Object[]</code> value
     * @param weight a <code>double</code> value
     * @param distanceMatrix a <code>DistanceMatrix</code> value
     */
    public NonBondedEnergyTerm(Updateable[] updateableResources,
                               double weight,
                               DistanceMatrix distanceMatrix){
        super(updateableResources, weight);
        this.distanceMatrix = distanceMatrix;
    }



    /**
     * Evaluates energy for each distance
     * @return a sum of all energy elements
     */
    public double evaluate() {
	if (! on) return 0.0;
	double energy = 0;
	DistanceList nonBondedList = distanceMatrix.nonBondedList();
	for (Distance distance:nonBondedList) {
            if (!distance.mode().frozen) {
	       energyElement.set(distance);
	       energy += energyElement.evaluate();
            }
        }
	return energy;
    }
    /**
     * Describe <code>evaluateAtoms</code> method here.
     */
    public void evaluateAtoms() {
	if (on) {
        DistanceList nonBondedList = distanceMatrix.nonBondedList();
        for(Distance nonBonded : nonBondedList){
                if (!nonBonded.mode().frozen) {
                    energyElement.set(nonBonded);
                    energyElement.evaluateAtoms();
                }
	    }
	}
    }

    /**
     * Testing of one atom in all energy elements
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom an criminal <code>Atom</code> value
     */
    public void test(TotalEnergy totalEnergy,Atom atom) {
	if (! on) System.out.println(""+this +" is off");
    DistanceList nonBondedList = distanceMatrix.nonBondedList();
        for(Distance nonBonded : nonBondedList){
            if (nonBonded == null) {
                    //System.out.println("Null distance in NonBondedList");
            System.out.println("Weird nonBondedList:");
            for(Distance d : nonBondedList) System.out.println(d);
            throw new RuntimeException("Null distance in NonBondedList");
            }
            if ((nonBonded.atom1() == null) || (nonBonded.atom2() == null)){
                    //System.out.println("Null distance in NonBondedList"+nonBonded);
                System.out.println("Weird nonBondedList:");
                for(Distance d : nonBondedList) System.out.println(d);
                throw new RuntimeException("Null atom in distance "+nonBonded);
            }

            energyElement.set(nonBonded);
            energyElement.test(totalEnergy,atom);
    }
    }
}
