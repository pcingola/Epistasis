package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.util.*;

import java.util.Iterator;

/** Atomic pairwise potential of mean force, original energy term by Chris
 * Summa.
 *
 * Reference:
 * Summa CM and Levitt M, Near-native structure refinement using in vacuo
 * energy minimization. PNAS 2007;104(9):3177-82.
 *
 * @author El-ad David Amir
 */
public class AtomicPairwisePMFSumma extends AbstractEnergy {
    public final DistanceMatrix  distanceMatrix;
    private boolean frozenFlag = false;
    private static AtomCore atom1, atom2;
    private static int bin;
    private static double d, energy0, energy1,x, x2, x3,fx,fy,fz, halfEnergy0, e;
    /* polynomial spline for this pair */
    private static CoefficientsMatrixForAtomPairSpline coefs;

    private static int    numberOfTypes = AtomType.values().length;
    private AtomicPairwisePMFSummaParameters parameters;
    protected static final double[] sumPerAtomType      = new double[numberOfTypes];
    protected static final double[] sm2PerAtomType      = new double[numberOfTypes];
    protected static       int[]    numberOfEachAtomType = null;
    public    static final double[] sumPerAtomType() { return sumPerAtomType;}
    public    static final double[] sm2PerAtomType() { return sm2PerAtomType;}
    public    static       int[]    numberOfEachAtomType() {return numberOfEachAtomType;} 
    public    static       int      atom1TypeOrdinal, atom2TypeOrdinal, atomTypeOrdinal;
    private static SummaAttribute summaAttribute;
    public static double[] atomEnergies;

    public AtomicPairwisePMFSumma( DistanceMatrix distanceMatrix, double weight, AtomicPairwisePMFSummaParameters parameters ) {
	    super( toArray( distanceMatrix ), weight);
	    this.distanceMatrix = distanceMatrix;
	    comment = "AtomicPairwisePMFSumma";
	    atomEnergies = new double[distanceMatrix.molecularSystem.size()];
	    for (AtomCore atom:distanceMatrix.molecularSystem) 
		if (atom.status() == AtomStatus.FROZEN) frozenFlag = true; 
	    this.parameters = parameters;
	    numberOfEachAtomType = new int[numberOfTypes];
	    for (int i = 0; i < sumPerAtomType.length; i++) {
		numberOfEachAtomType[i] = 0;
	    }
	    for (AtomCore atom:distanceMatrix.molecularSystem)
		numberOfEachAtomType[atom.type().ordinal()]++;
    }
   /**
     * Evaluates energy for each distance
     * @return a sum of all energy elements
     */
    public double evaluate() {
	return evaluate(false);
    }


    /**
     * Describe <code>evaluateAtoms</code> method here.
     */
    public void evaluateAtoms() {
	evaluate(true);
    }
    
    public double evaluate(boolean evaluateAtoms) {
	if (! on) return 0.0;
	double energy = 0;
	DistanceList nonBondedList = distanceMatrix.nonBondedList();
	if (nonBondedList == null) throw new RuntimeException("nonBondedList == null");
	for (int i = 0; i < atomEnergies.length; i++) atomEnergies[i] = 0;
	for (Distance distance:nonBondedList) {
            if (!distance.mode().frozen) {
	                    energy += evaluateElement(distance, evaluateAtoms);
            }
        }

	for (int i = 0; i < numberOfTypes; i++) {
	    sumPerAtomType[i] = 0;
	    sm2PerAtomType[i] = 0;
	}
	for (AtomCore atom:distanceMatrix.molecularSystem) {
	    atomTypeOrdinal = atom.type().ordinal();
	    e = atomEnergies[atom.number];
	    sumPerAtomType[atomTypeOrdinal] += e;
	    sm2PerAtomType[atomTypeOrdinal] += e*e;
	}

	return energy;
  }

    public double evaluateElement(Distance distance, boolean evaluateAtoms){
        atom1 = distance.atom1;
		atom2 = distance.atom2;
		atom1TypeOrdinal = atom1.type().ordinal();
		atom2TypeOrdinal = atom2.type().ordinal();
		coefs = parameters.splineForAtomPair( atom1TypeOrdinal, atom2TypeOrdinal);
		/* turn off term if any of the atoms was not found in parameters file */
		if (coefs == null) return 0;
		summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
		if (summaAttribute == null) {
		    summaAttribute = new SummaAttribute();
		    distance.addAttribute(summaAttribute);
		}
		/* calculate energy and derivative */
		d = distance.distance();
		bin = (int) (d*10);
// 		bin = coefs.bins.length-1;
// 		while( (bin > 0) && coefs.bins[bin] >= d )
// 			bin--;

		/* correct x according to its distance from start of bin */
		x = d - coefs.bins[bin];
		/* calculate powers of x */
		x2 = x*x;
		x3 = x2*x;

		/* calculate energy */
		try {
		    energy0 = coefs.coefs[bin][0]*x3 + coefs.coefs[bin][1]*x2 + coefs.coefs[bin][2]*x + coefs.coefs[bin][3];
		}
		catch (RuntimeException ex) {
		    System.out.println("Problem in evaluate\n"+
				       atom1+" "+atom1.type()+"\n"+
				       atom2+" "+atom2.type()+"\n"+
				       "bin = "+bin+"\n"+
				       "distance = "+d);
		    System.out.println("coefs = "+coefs);
		    System.out.println("coefs.coefs = "+coefs.coefs);
		    System.out.println("coefs.coefs = "+coefs.coefs);
		    throw ex;
		}
		/* calculate first derivative */
		energy1 = 3*coefs.coefs[bin][0]*x2 + 2*coefs.coefs[bin][1]*x + coefs.coefs[bin][2];

		/* apply weight */
		energy0     = energy0 * weight;
		energy1     = energy1 * weight;
		halfEnergy0 = energy0*0.5;
		/* apply force to atoms */
		summaAttribute.fx = fx = -energy1 * distance.dx*distance.invDistance;
		summaAttribute.fy = fy = -energy1 * distance.dy*distance.invDistance;
		summaAttribute.fz = fz = -energy1 * distance.dz*distance.invDistance;

		if (frozenFlag) {
		    if( !atom1.status().frozen() ) atom1.addForce( fx,  fy,  fz);
		    if( !atom2.status().frozen() ) atom2.addForce(-fx, -fy, -fz);
		}
		else {
		    atom1.addForce( fx,  fy,  fz);
		    atom2.addForce(-fx, -fy, -fz);
		}
	       if (evaluateAtoms) {
		   atom1.atom.addEnergy(halfEnergy0);
		   atom2.atom.addEnergy(halfEnergy0);
	       }
	       atomEnergies[atom1.number] += halfEnergy0;
	       atomEnergies[atom2.number] += halfEnergy0;
       return energy0;
    }

     public void test(TotalEnergy totalEnergy,Atom atom) {
//	 System.out.println("Sorry no test for  AtomicPairwisePMFSumma yet");
     if (! on) {System.out.println(""+this +" is off"); return;}
	System.out.println("Testing "+this+" using "+atom);
        if (atom == null)
	    throw new RuntimeException("Cannot test "+this);

         DistanceList nonBondedList = distanceMatrix.nonBondedList();
         if (nonBondedList == null) throw new RuntimeException("nonBondedList == null");
         for (int i = 0; i < atomEnergies.length; i++) atomEnergies[i] = 0;

         for (Distance distance:nonBondedList) {
             if (distance.mode().frozen) continue;
             if ((distance.atom1() != atom) &&(distance.atom2() != atom)) continue;
             
             //energyElement test
             double[][] coordinates = new double[3][];
             coordinates[0] = atom.X();
             coordinates[1] = atom.Y();
             coordinates[2] = atom.Z();
             for(int i = 0; i< 3; i++) {
                 try{totalEnergy.update();}catch(UpdateableException ue){}
                 double x = coordinates[i][0];
                 coordinates[i][1] = 0;
                 double e1 = evaluateElement(distance, false);
                 double analiticalForce = coordinates[i][1];
                 coordinates[i][0] += EnergyElement.DX;
                 // Whatever should be updated ( such as distance matrix torsion list etc. )
                 try{totalEnergy.update();}catch(UpdateableException ue){}
                 double e2 = evaluateElement(distance, false);
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

}
