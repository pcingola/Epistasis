package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.*;
import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.parameters.*;

public class CooperativeAtomicPairwisePMFSumma extends AbstractEnergy{
    private static final double[] targetSigma       = new double[AtomType.values().length];
    private static final double[] diffSigma         = new double[AtomType.values().length];
    public  static final int      numberOfAtomTypes = AtomType.values().length;
    private  static DistanceMatrix distanceMatrix;

    private static SummaAttribute summaAttribute;
    private static AtomCore atom1, atom2;
    private static AtomType type1, type2;
    private static double   diff1, diff2, factor;
    private AtomicPairwisePMFSumma atomicPairwisePMFSumma;

    public CooperativeAtomicPairwisePMFSumma(DistanceMatrix distanceMatrix, double weight, CooperativeAtomicPairwisePMFSummaParameters parameters,
					     AtomicPairwisePMFSumma  atomicPairwisePMFSumma) {
	super(toArray(distanceMatrix), weight);
	if (weight != 0) {
	    this.atomicPairwisePMFSumma = atomicPairwisePMFSumma; 
	    if (atomicPairwisePMFSumma  == null) throw new RuntimeException("CooperativeAtomicPairwisePMFSumma term must not be created "+
									    "before the non-cooperative term is created.");
	    comment = "CooperativeAtomicPairwisePMFSumma";
	    for (AtomType type:AtomType.values()) {
		int n = atomicPairwisePMFSumma.numberOfEachAtomType[type.ordinal()];
		targetSigma[type.ordinal()] = parameters.minAvgSigma[type.ordinal()] * n;
	    }
	    this.distanceMatrix = distanceMatrix; 
	}
    }
    
    public double evaluate() {
	if (weight == 0) return 0;
	return evaluate(false);
    }
    public void evaluateAtoms() {
	evaluate(true);
    }
    
    public double evaluate(boolean evaluateAtoms) {
	if (! on) return 0.0;
	DistanceList nonBondedList = distanceMatrix.nonBondedList();
	double energy = 0;

	for (int i = 0; i < numberOfAtomTypes; i++) {
	    diffSigma[i] = atomicPairwisePMFSumma.sumPerAtomType[i]-targetSigma[i];
	    if (diffSigma[i] > 0 ) diffSigma[i] = 0;
	    energy += weight * diffSigma[i] * diffSigma[i];
	    if (evaluateAtoms) {
		AtomType type = AtomType.values()[i];
		int      n    = AtomicPairwisePMFSumma.numberOfEachAtomType[i];
		for (AtomCore atom:distanceMatrix.molecularSystem)
		    if (atom.type() == type) atom.atom.addEnergy((weight * diffSigma[i]* diffSigma[i])/(n*n));
	    }
	}
	    

	for (Distance distance:nonBondedList) {
            if (!distance.mode().frozen) {
		summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
		if (summaAttribute == null) continue; 
		atom1     = distance.atom1;
		atom2     = distance.atom2;
		type1 = atom1.type();
		type2 = atom2.type();
        //int      n1    = AtomicPairwisePMFSumma.numberOfEachAtomType[type1.ordinal()];
        //int      n2    = AtomicPairwisePMFSumma.numberOfEachAtomType[type2.ordinal()];
		//diff1     = diffSigma[type1.ordinal()]/(n1*n1);
		//diff2     = diffSigma[type2.ordinal()]/(n2*n2);
         diff1     = diffSigma[type1.ordinal()];
         diff2     = diffSigma[type2.ordinal()];
        factor    = (diff1 + diff2)*weight;            // (2* diff1 + 2*diff2)/2 * weight                // factor    = (diff1 + diff2)*weight/n           - WHY NOT SO?
		if (factor != 0) {
		    if( !atom1.status().frozen() ) atom1.addForce(factor*summaAttribute.fx,  factor*summaAttribute.fy,  factor*summaAttribute.fz);
		    if( !atom2.status().frozen() ) atom2.addForce(-factor*summaAttribute.fx, -factor*summaAttribute.fy, -factor*summaAttribute.fz);
		}
	    }
	}
	return energy;
    }

    public void test(TotalEnergy energy, Atom atom) {
//	throw new RuntimeException("Cannot test");
        	System.out.println("Cannot test "+this);
    }
}


		
	