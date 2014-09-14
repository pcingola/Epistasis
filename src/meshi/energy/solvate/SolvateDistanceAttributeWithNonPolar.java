package meshi.energy.solvate;
import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.util.mathTools.*;


 
public class SolvateDistanceAttributeWithNonPolar  extends SolvateDistanceAttribute {   
    public final Distance dis;
    public final AtomCore atom1, atom2;
    public final int atom1number, atom2number;

    public final int TsaiAtomicType1;
    public final int TsaiAtomicType2;
    public final double Cend_12, Cp1_12, Cp2_12, CvalAtp1_12, CvalAtp2_12;
    public final double Cend_21, Cp1_21, Cp2_21, CvalAtp1_21, CvalAtp2_21;


     // Carbon sigmoids and related derivatives
    public double sigmCa1;      /*1*/
    public double dsigmCa1dx1;  /*2*/
     public double dsigmCa1dy1;/*2*/
     public double dsigmCa1dz1;/*2*/
     public double dsigmCa1dx2;/*2*/
     public double dsigmCa1dy2;/*2*/
     public double dsigmCa1dz2;/*2*/
    public double sigmCa2;     /*1*/
     public double dsigmCa2dx1;/*2*/
     public double dsigmCa2dy1;/*2*/
     public double dsigmCa2dz1;/*2*/
     public double dsigmCa2dx2;/*2*/
     public double dsigmCa2dy2;/*2*/
     public double dsigmCa2dz2;/*2*/

    public SolvateDistanceAttributeWithNonPolar(Distance dis, AtomCore atom1, AtomCore atom2, int atom1number, int atom2number, 
				    boolean isDisInvolvesHydrogens, boolean isDisBetween2Polars,
				    int TsaiAtomicType1, int TsaiAtomicType2,
				    double Cend_12, double Cp1_12, double Cp2_12, double CvalAtp1_12, double CvalAtp2_12,
				    double Cend_21, double Cp1_21, double Cp2_21, double CvalAtp1_21, double CvalAtp2_21) {
	super(false, false);
	this.dis         = dis;
	this.atom1       = atom1;
	this.atom2       = atom2;
	this.atom1number = atom1number;
	this.atom2number = atom2number;
	this.TsaiAtomicType1 = TsaiAtomicType1;
	this.TsaiAtomicType2 = TsaiAtomicType2;
	this.Cend_12 = Cend_12;
	this.Cp1_12 = Cp1_12;
	this.Cp2_12 = Cp2_12;
	this.CvalAtp1_12 = CvalAtp1_12;
	this.CvalAtp2_12 = CvalAtp2_12;
 	this.Cend_21 = Cend_21;
	this.Cp1_21 = Cp1_21;
	this.Cp2_21 = Cp2_21;
	this.CvalAtp1_21 = CvalAtp1_21;
	this.CvalAtp2_21 = CvalAtp2_21;
	resetSigmVals();
   }
    public final void update() {
       	double sigmas_tag;

     	// Handling Carbons
    	// ----------------
	// Is atom2 non-polar? If so it should contribute to atom1's CNC
    	if (atom2.type().isCarbon() || (atom2.type() == AtomType.MSD)) {
	    Sigma.sigma(dis.distance(),Cend_12, Cp1_12, Cp2_12, CvalAtp1_12, CvalAtp2_12);	    
	    sigmCa1    = Sigma.s;
	    sigmas_tag = Sigma.s_tag;
	    dsigmCa1dx1 = sigmas_tag*dis.dDistanceDx();
	    dsigmCa1dy1 = sigmas_tag*dis.dDistanceDy();
	    dsigmCa1dz1 = sigmas_tag*dis.dDistanceDz();
	    dsigmCa1dx2 = -dsigmCa1dx1;
	    dsigmCa1dy2 = -dsigmCa1dy1;
	    dsigmCa1dz2 = -dsigmCa1dz1;
	}
	// Is atom1 non-polar? If so it should contribute to atom2's CNC
    	if (dis.atom1().type().isCarbon() || (dis.atom1().type() == AtomType.MSD)) {	
	    Sigma.sigma(dis.distance(), Cend_21, Cp1_21, Cp2_21, CvalAtp1_21, CvalAtp2_21);
	    sigmCa2 = Sigma.s;
	    sigmas_tag = Sigma.s_tag;	    
	    dsigmCa2dx1 = sigmas_tag*dis.dDistanceDx();
	    dsigmCa2dy1 = sigmas_tag*dis.dDistanceDy();
	    dsigmCa2dz1 = sigmas_tag*dis.dDistanceDz();
	    dsigmCa2dx2 = -dsigmCa2dx1;
	    dsigmCa2dy2 = -dsigmCa2dy1;
	    dsigmCa2dz2 = -dsigmCa2dz1;
	}
    } // of update





	public final void resetSigmVals() {
	    sigmCa1 = 0.0;
	    dsigmCa1dx1 = 0.0;
	    dsigmCa1dy1 = 0.0;
	    dsigmCa1dz1 = 0.0;
	    dsigmCa1dx2 = 0.0;
	    dsigmCa1dy2 = 0.0;
	    dsigmCa1dz2 = 0.0;
	    sigmCa2 = 0.0;
	    dsigmCa2dx1 = 0.0;
	    dsigmCa2dy1 = 0.0;
	    dsigmCa2dz1 = 0.0;
	    dsigmCa2dx2 = 0.0;
	    dsigmCa2dy2 = 0.0;
	    dsigmCa2dz2 = 0.0;
     	}
}
