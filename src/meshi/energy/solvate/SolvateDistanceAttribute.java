package meshi.energy.solvate;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.geometry.*;
import meshi.parameters.*;
import meshi.molecularElements.hydrogenBonds.*;
import meshi.util.mathTools.*;
/** 
 * In the SolvateEnergy evaluation function these fields are to be recalculated for EVERY 
 * distance in the non-bonded list. They are used several times in the solvate evaluation
 * process, and we would like to calculate them once. To this end we attach this special 
 * class as an attribute on each Distance instance. 
 * For each distance we calculate two sigmoid values, and an hydrogen-bond strength value: 
 * sigmCa1 - The carbon index of atom 2 on atom 1. If atom 2 is not a carbon then this value
 *           should be zero. If atom 2 is a carbon then this value should be ~1.0 if atom 2 
 *           is spatially near atom 1. This index drops sigmoidally to zero the farther 
 *           atom 2 is.  
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atom 1 on atom 2.
 * sigmHBa1 - The hydrogen bond (HB) strength between atom 1 and 2. If atom 2 can not form HB 
 *           with atom 1, because of its chemical type then this value should be zero.
 *           If atom 2 can form HB with atom 1 then this value should be ~1.0 if atom 1 and 2
 *           are sufficiently close to create a HB, and if their orientation (defined
 *           also by their base atoms - see HydrogenBondDahiyat) permit hydrogen bonding.
 *           This value drops steeply to zero if the conditions to hydrogen bonding 
 *           are violated.
 *              
 * Also provided are the carbon sigmoid values derivative relatives to the atom coordinates. They 
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 *   
 **/


public class SolvateDistanceAttribute  implements MeshiAttribute {

    // Hydogen involvement
    public final boolean isDisInvolvesHydrogens; /*2*/
    public final boolean isDisBetween2Polars;   // Is this a polar pair? /*1*/	/*2*/    
    public SolvateDistanceAttribute(boolean isDisInvolvesHydrogens,  boolean isDisBetween2Polars) {
	this.isDisInvolvesHydrogens = isDisInvolvesHydrogens;
	this.isDisBetween2Polars    = isDisBetween2Polars;
	
    }
    public int key() {return SOLVATE_ALL_ATOM_ATTRIBUTE;}

}












