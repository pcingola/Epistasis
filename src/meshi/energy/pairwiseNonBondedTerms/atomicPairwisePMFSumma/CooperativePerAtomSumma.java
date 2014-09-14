package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceList;
import meshi.geometry.Distance;
import meshi.parameters.AtomType;
import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.util.MeshiAttribute;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 01/06/2009
 * Time: 12:46:56
 * To change this template use File | Settings | File Templates.
 */
public class CooperativePerAtomSumma  extends AbstractEnergy {
//    private static final double DEFAULT_KOEFF_OUTSTD = 1;
 //   private static final double DEFAULT_KOEFF_INSTD = 1;
    private AtomicPairwisePMFSumma atomicPairwisePMFSumma;
    private static final double[] mean       = new double[AtomType.values().length];
    private static final double[] std         = new double[AtomType.values().length];
    private  static DistanceMatrix distanceMatrix;

    private static SummaAttribute summaAttribute;
    private static AtomCore atom1, atom2;
    private static AtomType type1, type2;
    private static double   diff1, diff2, factor;

   // private double koeffOutStd = DEFAULT_KOEFF_INSTD;
    

    public CooperativePerAtomSumma(DistanceMatrix distanceMatrix, double weight, CooperativePerAtomSummaParameters parameters,
					     AtomicPairwisePMFSumma  atomicPairwisePMFSumma) {
	super(toArray(distanceMatrix), weight);
	if (weight != 0) {
	    this.atomicPairwisePMFSumma = atomicPairwisePMFSumma;
	    if (atomicPairwisePMFSumma  == null) throw new RuntimeException("CooperativePerAtomSumma term must not be created "+
									    "before the non-cooperative term is created.");
	    comment = "CooperativePerAtomSumma";
	    for (AtomType type:AtomType.values()) {
		    mean[type.ordinal()] = parameters.mean[type.ordinal()];
            std[type.ordinal()] = parameters.std[type.ordinal()];
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
	double energy = 0, d, e, eAt;


        for (AtomCore atom:distanceMatrix.molecularSystem){
            AtomType type = atom.type();
            //d = atomicPairwisePMFSumma.atomEnergies[atom.number] - mean[type.ordinal()];
          eAt =   atomicPairwisePMFSumma.atomEnergies[atom.number];
  /*          if (Math.abs(mean[type.ordinal()] - eAt) < std[type.ordinal()] ){
                        koeffOutStd = 1;
                    }
                    else  koeffOutStd = DEFAULT_KOEFF_OUTSTD;
    */
            d = zScore(eAt, mean[type.ordinal()], std[type.ordinal()]);
            //e = weight * d*d;
      //      e = weight * koeffOutStd * (d*d - 1);
           // e = weight * koeffOutStd * d*d;
             e = weight * d*d;
            energy += e;
		    atom.atom.addEnergy(e);
	    }

	for (Distance distance:nonBondedList) {
      if (!distance.mode().frozen) {
        summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
		if (summaAttribute == null) continue;
		atom1     = distance.atom1;
		atom2     = distance.atom2;
		type1 = atom1.type();
		type2 = atom2.type();
//        diff1 = atomicPairwisePMFSumma.atomEnergies[atom1.number] - mean[type1.ordinal()];
//		diff2 = atomicPairwisePMFSumma.atomEnergies[atom2.number] - mean[type2.ordinal()];
    /*      eAt =   atomicPairwisePMFSumma.atomEnergies[atom1.number];
          if (Math.abs(mean[type1.ordinal()] - eAt) < std[type1.ordinal()] ){
                      koeffOutStd = 1;
                  }
                  else  koeffOutStd = DEFAULT_KOEFF_OUTSTD;
      */
          diff1 = zScore(atomicPairwisePMFSumma.atomEnergies[atom1.number], mean[type1.ordinal()], std[type1.ordinal()])/std[type1.ordinal()];
       //   diff1 *= koeffOutStd;
       /*
          eAt =   atomicPairwisePMFSumma.atomEnergies[atom2.number];
          if (Math.abs(mean[type2.ordinal()] - eAt) < std[type2.ordinal()] ){
                      koeffOutStd = 1;
                  }
                  else  koeffOutStd = DEFAULT_KOEFF_OUTSTD;
         */
          diff2 = zScore(atomicPairwisePMFSumma.atomEnergies[atom2.number], mean[type2.ordinal()], std[type2.ordinal()])/std[type2.ordinal()];
     //     diff2 *= koeffOutStd;
        factor    = (diff1 + diff2)*weight;            // (2*diff1 + 2*diff2)/2*weight
		if (factor != 0) {
		    if( !atom1.status().frozen() ) atom1.addForce(factor*summaAttribute.fx,  factor*summaAttribute.fy,  factor*summaAttribute.fz);
		    if( !atom2.status().frozen() ) atom2.addForce(-factor*summaAttribute.fx, -factor*summaAttribute.fy, -factor*summaAttribute.fz);
		}
	    }
	}
	return energy;
    }

    private double zScore(double curE, double mean, double std){
       // double koeffOutStd = 1;
        if (std == 0)
                throw new RuntimeException("Something weird in Parameters "+this);
 /*       if (Math.abs(mean - curE) < std ){
            koeffOutStd = DEFAULT_KOEFF_INSTD;
        }
        else
            koeffOutStd = DEFAULT_KOEFF_OUTSTD;
   */
        return (curE-mean)/std;
    }

    public void test(TotalEnergy energy, Atom atom) {
//	throw new RuntimeException("Cannot test");
        	System.out.println("Cannot test "+this);
    }
}
