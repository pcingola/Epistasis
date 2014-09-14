package meshi.energy.pairwiseNonBondedTerms.excludedVolOLD;

import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;


public  class ExcludedVolEnergyElementOLD extends NonBondedEnergyElement {
    /** 
    * ALPHA is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
    **/    
    public final double  ALPHA = 0.2; 
    /**
    * C is the parameter of  EV = C*(d-sigma)^4 in the range [0,sigma]
    **/    
    public final double  C = 1/(ALPHA*ALPHA*ALPHA*ALPHA);

    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected Distance distance;
    protected double sigma,Rfac;
    protected boolean frozen;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected double rMax;
    protected double[][] parameters;

    public ExcludedVolEnergyElementOLD() {}
    public ExcludedVolEnergyElementOLD(DistanceMatrix distanceMatrix,
                                     double Rfac,
                                     double weight,
                                     double[][] parameters) {
	this.parameters = parameters;
	this.weight = weight;
	this.distanceMatrix = distanceMatrix;
	this.Rfac = Rfac;
	this.rMax = distanceMatrix.rMax();
        //----------------------
    }
    protected void setAtoms(){
	throw new RuntimeException("setAtoms() may not be used by ExcludedVolEnergyElementOLD for "+
				   "efficiency.");
    }
    
	
    public void set(Object obj) {
	distance = (Distance) obj;
	//atoms = distance.atoms();
	atom1 = distance.atom1();
	atom2 = distance.atom2();
	sigma = parameters[atom1.type().ordinal()][atom2.type().ordinal()];
	if (sigma == Double.MAX_VALUE) 
	    throw new RuntimeException("No parameters for "+atom1+" and "+atom2); 
	if (sigma>rMax) 
	    throw new RuntimeException("Excluded Volume: sigma="+sigma+" and it is larger " +
				       "than rMax="+rMax);
    }

        
    public double evaluate() {
	double dEdD;
	double DMinusSig;
	double D3;
	double localSigma=sigma;
	double dis = distance.distance();
	if (Rfac<0.99) {
	   if (atom1.residueNumber() == atom2.residueNumber())
	      localSigma = 0;
	   else
	      localSigma = sigma*Rfac;
    }
	if (dis>localSigma) {
	    energy = dEdD = 0;
	}
	else {
		DMinusSig = (dis-localSigma);
		D3 = C*DMinusSig*DMinusSig*DMinusSig*weight;
	    energy = D3*(dis-localSigma);
	    dEdD = 4*D3;
		if (! atom1.frozen()) {
			atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
			atom1.addToFy(-1*dEdD*distance.dDistanceDy());
			atom1.addToFz(-1*dEdD*distance.dDistanceDz());
		}
		if (! atom2.frozen()) {
			atom2.addToFx(dEdD*distance.dDistanceDx());
			atom2.addToFy(dEdD*distance.dDistanceDy());
			atom2.addToFz(dEdD*distance.dDistanceDz());
		}
	 }	
//if (energy>10.0)
//       System.out.println(atom1.residueNumber+" "+atom2.residueNumber+" "+atom1.type+" "+atom2.type+" "+sigma+" "+dis+" "+A+" "+B+" "+C+" "+ALPHA+" "+sigmaMinusALPHA);
	return energy;
    }
    

    public String toString() {
	if ((atom1 == null) & (atom2 == null)) return "ExcludedVolumeEnergyElement - atoms not set yet";
	if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
									  "atom1 = "+atom1+"\n"+
									  "atom2 = "+atom2); 
	double dis = distance.distance();
	return ("ExcludedVolumeEnergyElement sigma = "+sigma+" Distance = "+
		dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
    }
}
