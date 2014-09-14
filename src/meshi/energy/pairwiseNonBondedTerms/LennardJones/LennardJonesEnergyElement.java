package meshi.energy.pairwiseNonBondedTerms.LennardJones;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.parameters.*;
import meshi.util.filters.Filter;
import java.util.*;

public  class LennardJonesEnergyElement extends NonBondedEnergyElement {
    public static final double MAX_ENERGY = 10000;
    public static final double ALPHA = 0.5;
    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected int atom1Number, atom2Number;    
    protected double epsilon, sigma, sigma6, sigma6EpsilonFour, minusTwelveSigma6;
    protected boolean frozen;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected double maxEnergy;
    protected double breakEnergy;
    protected double breakEnergy4;
    protected double breakEnergySquare4;
    protected double rMax;
    protected double contact, dCdD, dCdX, dCdY, dCdZ;
    private final int FIRST = 0,SECOND=1;
    private Distance distance;
 
    protected static double[][][] parameters;

    public  LennardJonesEnergyElement() {}
    public  LennardJonesEnergyElement(DistanceMatrix distanceMatrix, double weight, double[][][] parameters) {
	this.parameters = parameters;
	this.weight = weight;
	this.distanceMatrix = distanceMatrix;
	maxEnergy = MAX_ENERGY/weight;
	breakEnergy = maxEnergy/3;
	breakEnergy4 = breakEnergy*4;
	breakEnergySquare4 =  breakEnergy4*breakEnergy;
	rMax = distanceMatrix.rMax();
        //----------------------
    }
    protected void setAtoms(){
	throw new RuntimeException("setAtoms() may not be used by LennardJonesEnergyElement for "+
				   "efficiency.");
    }
    
    public double distance() {
	return (distanceMatrix.distance(atom1Number, atom2Number)).distance();
    }
	
    public void set(Object obj) {
	distance = (Distance) obj;
	atoms = new AtomList(2);
	atom1 = distance.atom1();
	atom2 = distance.atom2();
	atoms.add(atom1);
	atoms.add(atom2);
	atom1Number = atom1.number();
	atom2Number = atom2.number();
	epsilon = parameters[atom1.type().ordinal()][atom2.type().ordinal()][0];
	sigma = parameters[atom1.type().ordinal()][atom2.type().ordinal()][1];
	if ((sigma == Double.MAX_VALUE) || (epsilon == Double.MAX_VALUE)) 
	    throw new RuntimeException("No parameters for "+atom1+" and "+atom2); 
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
	sigma6EpsilonFour = 4*epsilon*sigma6;
	minusTwelveSigma6 = -12*sigma6;
    }        
    
    public double evaluate() {
        updateEnergy();
	updateAtoms();
	return energy*weight;
    }
    
    public double updateEnergy() {
	double LJ;
	double LJPlus;
	double dLJdD;
	double invD = 0;
	double invD2;
	double invD6 = 0;
	double invD7;
	double invD12 = 0;
	double invD13;
	double dis = -1;
	double  rMaxMinusDis, rMaxMinusDisPlus;
	double energy1, dE1dD;
            dis = distance.distance();
            rMaxMinusDis = rMax - dis;
	    invD = distance.invDistance();		
	    invD2 = invD*invD;
	    invD6 = invD2*invD2*invD2;
	    invD7 = invD6*invD;
	    invD12 = invD6*invD6;
	    invD13 = invD12*invD;
	    LJ = sigma6EpsilonFour * (sigma6*invD12 - invD6);
	    dLJdD =  sigma6EpsilonFour * ( minusTwelveSigma6*invD13 + 6*invD7);

	    //smooth high energies
	    if (LJ >= breakEnergy) {
		LJPlus = LJ+breakEnergy;
		energy1 = breakEnergy4*LJ/LJPlus-breakEnergy;
		dE1dD =breakEnergySquare4/(LJPlus*LJPlus)*dLJdD;
	    }
	    else {
		energy1 = LJ;
		dE1dD = dLJdD;
	    }

	    //quench to zero in rMax
	    double rMaxMinusDisSquare = rMaxMinusDis*rMaxMinusDis;
	    double rMaxMinusDisSquarePlusAlpha = rMaxMinusDisSquare+ALPHA;
	    double rMaxMinusDisTimesAlpha = rMaxMinusDis*ALPHA;
	    double rMaxMinusDisSquarePlusAlphaSquare = rMaxMinusDisSquarePlusAlpha*rMaxMinusDisSquarePlusAlpha;
	    contact = rMaxMinusDisSquare/rMaxMinusDisSquarePlusAlpha;
	    dCdD = -2*rMaxMinusDisTimesAlpha/rMaxMinusDisSquarePlusAlphaSquare; //The minus comes from rMaxMinusDis
	    energy = energy1*contact;
	    dEdD = dE1dD*contact + dCdD*energy1;	    

	    dEdX = dEdD*distance.dDistanceDx();
	    dEdY = dEdD*distance.dDistanceDy();
	    dEdZ = dEdD*distance.dDistanceDz();
	    //	if ((! (energy < 0) ) & (! ( energy == 0)) & (!(energy > 0)))System.out.println("weird energy "+energy+" "+distance);
	return energy;
    }        
    
    public void updateAtoms(){
	if (! atom1.frozen()) {
	    atom1.addToFx(-1*dEdX*weight); // force = -derivative   
	    atom1.addToFy(-1*dEdY*weight);
            atom1.addToFz(-1*dEdZ*weight);
	}
	if (! atom2.frozen()) {
	    atom2.addToFx(dEdX*weight);
	    atom2.addToFy(dEdY*weight);
	    atom2.addToFz(dEdZ*weight);
	}
    }

     public double contactAtoms(double weight) {
       double e = contact(weight);
       atom1.addEnergy(e/2);
       atom2.addEnergy(e/2);
       return e;
     }
	
    
    public double dEdXAtom(int index){
        if (index == FIRST)
            return dEdX;
        if (index == SECOND)
            return -1*dEdX;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }
    public double dEdYAtom(int index){
        if (index == FIRST)
            return dEdY;
        if (index == SECOND)
            return -1*dEdY;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }
    public double dEdZAtom(int index){
        if (index == FIRST)
            return dEdZ;
        if (index == SECOND)
            return -1*dEdZ;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }
           
     public double contact(double weight) {
       if (! atom1.frozen()) {
           atom1.addToFx(dCdX*weight); // force = -derivative
           atom1.addToFy(dCdY*weight);
           atom1.addToFz(dCdZ*weight);
       }
       if (! atom2.frozen()) {
           atom2.addToFx(-1.0*dCdX*weight);
           atom2.addToFy(-1.0*dCdY*weight);
           atom2.addToFz(-1.0*dCdZ*weight);
       }
       return -1.0*weight*contact;
     }
 
   public double contact() {return contact;}



    public String toString() {
	if ((atom1 == null) & (atom2 == null)) return "LennardJonesEnergyElement - atoms not set yet";
	if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
									  "atom1 = "+atom1+"\n"+
									  "atom2 = "+atom2); 
	Distance distance = distanceMatrix.distance(atom1Number, atom2Number);
	double dis = distance.distance();
	return ("LennardJonesEnergyElement sigma = "+sigma+" epsilon = "+epsilon+" Distance = "+
		dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
    }

   
    /**
     * @return Returns the dEdD.
     */
    public final double getDEdD() {
	return dEdD;
    }
    /**
     * @return Returns the dEdX.
     */
    public final double getDEdX() {
	return dEdX;
    }
    /**
     * @return Returns the dEdY.
     */
    public final double getDEdY() {
	return dEdY;
    }
    /**
     * @return Returns the dEdZ.
     */
    public final double getDEdZ() {
	return dEdZ;
    }

}
