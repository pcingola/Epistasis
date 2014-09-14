package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.energy.EnergyElement;
import meshi.parameters.AtomType;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceList;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.Atom;
import meshi.util.MeshiAttribute;
import meshi.util.UpdateableException;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 04/06/2009
 * Time: 12:00:16
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZSumma extends AbstractEnergy{

    private AtomicPairwisePMFSumma atomicPairwisePMFSumma;
    private static final double[] mean       = new double[AtomType.values().length];
    private static final double[] std         = new double[AtomType.values().length];
    
    private  static DistanceMatrix distanceMatrix;

    private static SummaAttribute summaAttribute;
    private static AtomCore atom1, atom2;
    private static AtomType type1, type2;
    private static double   diff1, diff2, factor;


    public CooperativeZSumma(DistanceMatrix distanceMatrix, double weight, CooperativeZParameters parameters,
					     AtomicPairwisePMFSumma  atomicPairwisePMFSumma) {
	super(toArray(distanceMatrix), weight);
	if (weight != 0) {
	    this.atomicPairwisePMFSumma = atomicPairwisePMFSumma;
	    if (atomicPairwisePMFSumma  == null) throw new RuntimeException("CooperativeAtomicPairwisePMFSumma term must not be created "+
									    "before the non-cooperative term is created.");
        comment = "CooperativeZSumma";
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
	double energy = 0, e = 0;
        int nAtoms;
        nAtoms = distanceMatrix.molecularSystem.size();
        AtomType type;
        for (AtomCore atom:distanceMatrix.molecularSystem){
            type = atom.type();
   /*          double val = (atomicPairwisePMFSumma.atomEnergies[atom.number]-mean[type.ordinal()])/std[type.ordinal()];
            int a;
             if (atom.number == 1131)
             a = 1;
           if (val > 100)
                    a = 100;
            else
            if (val > 80)
                    a = 80;
            else
            if (val > 60)
                    a = 60;
            else
            if (val > 50)
                    a = 50;
        */


            e += zScore(atomicPairwisePMFSumma.atomEnergies[atom.number], mean[type.ordinal()], std[type.ordinal()]);
        }
        if (e >0) return 0;
        energy =e*e*weight;

        if (evaluateAtoms) {
            for (AtomCore atom:distanceMatrix.molecularSystem){
                atom.atom.addEnergy(energy/nAtoms);
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
                      if ((std[type1.ordinal()] == 0.0) || (std[type2.ordinal()] == 0.0))
                                        throw new RuntimeException("Something is weird in Parameters "+this);
                int a;
                 if ((atom1.number == 1131) || (atom2.number == 1131))
                 a = 1;

        //int      n1    = AtomicPairwisePMFSumma.numberOfEachAtomType[type1.ordinal()];
        //int      n2    = AtomicPairwisePMFSumma.numberOfEachAtomType[type2.ordinal()];
		//diff1     = diffSigma[type1.ordinal()]/(n1*n1);
		//diff2     = diffSigma[type2.ordinal()]/(n2*n2);

            diff1     =e/std[type1.ordinal()];
            diff2     =e/std[type2.ordinal()];
            factor    = (diff1 + diff2)*weight;            // (2* diff1 + 2*diff2)/2 * weight
		    if (factor != 0) {
		     if(!atom1.status().frozen() ) atom1.addForce(factor*summaAttribute.fx,  factor*summaAttribute.fy,  factor*summaAttribute.fz);
		     if( !atom2.status().frozen() ) atom2.addForce(-factor*summaAttribute.fx, -factor*summaAttribute.fy, -factor*summaAttribute.fz);
		    }
	    }
	}
	return energy;
    }


    private double zScore(double curE, double mean, double std){
       // double koeffOutStd = 1;
        if (std == 0)
                return 0;
      //      throw new RuntimeException("Something weird in Parameters "+this);
 /*       if (Math.abs(mean - curE) < std ){
            koeffOutStd = DEFAULT_KOEFF_INSTD;
        }
        else
            koeffOutStd = DEFAULT_KOEFF_OUTSTD;
   */

        return (curE-mean)/std;

    }


    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (! on) {System.out.println(""+this +" is off"); return;}
	System.out.println("Testing "+this+" using "+atom);
        if (atom == null)
	    throw new RuntimeException("Cannot test "+this);

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for(int i = 0; i< 3; i++) {
            try{totalEnergy.update();}catch(UpdateableException ue){}
	        atomicPairwisePMFSumma.evaluate();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
            try{totalEnergy.update();}catch(UpdateableException ue){}
	        atomicPairwisePMFSumma.evaluate();
            double e2 = evaluate();
            double de = e2-e1;
            double numericalForce = - de/ EnergyElement.DX;
            coordinates[i][0] -=  EnergyElement.DX;
            try{totalEnergy.update();}catch(UpdateableException ue){}

            double diff = Math.abs(analiticalForce - numericalForce);

            if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+EnergyElement.VERY_SMALL)) > EnergyElement.relativeDiffTolerance){
                System.out.println("Testing "+this);
                System.out.println("Atom["+atom.number()+"]."+EnergyElement.XYZ.charAt(i)+" = "+x);
                System.out.println("Analytical force = "+analiticalForce);
                System.out.println("Numerical force  = "+numericalForce);

                System.out.println("diff = "+diff+"\n"+
                                   "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+EnergyElement.VERY_SMALL) = "+
                                   2*diff/(Math.abs(analiticalForce) + Math.abs(numericalForce)+EnergyElement.VERY_SMALL));
                System.out.println();
            }
            if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\ne1 = "+e1);
            if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\ne2 = "+e2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\nanaliticalForce = "+analiticalForce);
        }

    }
}

