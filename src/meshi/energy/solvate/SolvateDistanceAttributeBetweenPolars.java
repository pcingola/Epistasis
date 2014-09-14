package meshi.energy.solvate;
import meshi.molecularElements.hydrogenBonds.*;
import meshi.geometry.*;
public class SolvateDistanceAttributeBetweenPolars extends SolvateDistanceAttribute {
    public  double sigmHB;  // The hydrogen bond value (including distance and angular dependence  /*1*/ /*2*/
    public final double saltBridgeFactorA1;     // Salt bridges vs. hydrogen bond strengths /*1*/ /*2*/
    public final double saltBridgeFactorA2;             /*1*//*2*/
    public final double saltBridgeFactorForHBenergyA1;  /*1*//*2*/   
    public final double saltBridgeFactorForHBenergyA2;  /*1*//*2*/  
    public final int atom1number, atom2number;

    public AbstractHydrogenBondList solvateHB;
    public SolvateDistanceAttributeBetweenPolars(Distance distance, AbstractHydrogenBondList solvateHB,
						 double saltBridgeFactorForHBenergyA1, double saltBridgeFactorForHBenergyA2, double saltBridgeFactorA1, double saltBridgeFactorA2) {
	super(false, true);
	this.solvateHB   = solvateHB;
	this.saltBridgeFactorForHBenergyA1 = saltBridgeFactorForHBenergyA1;
	this.saltBridgeFactorForHBenergyA2 = saltBridgeFactorForHBenergyA2;
	this.saltBridgeFactorA1            = saltBridgeFactorA1;
	this.saltBridgeFactorA2            = saltBridgeFactorA2;
	this.atom1number = distance.atom1.number;
	this.atom2number = distance.atom2.number;
   }

    public final void update() {
	// Is this a hydrogen bond or salt bridge?  Possible HB Ahoy!!
    	// -----------------------------------------------------------
	AbstractHydrogenBond hb = solvateHB.findBondByPolars(atom1number, atom2number);
	if (hb != null) 
	    sigmHB = hb.hbVal();
	else
	    sigmHB = 0.0;
	return;
    }	
}
