/*
 * Created on 22/8/2005
 *
 */
package meshi.energy.hydrogenBond;

import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.*;
import meshi.energy.pairwiseNonBondedTerms.*;

/**
 * @author amilev
 *
 * This class is an energy element for an Hydrogen Bond. Two atoms are involved: Hydrogen atom (hAtom,
 * and Oxygen Atom (oAtom).
 * The Hydrogen Bond forces energy between hAtom and oAtom is similar to
 * Lennard-Jones in the range [minimum:infinity]
 * where the x value of the minimum is sqrSixOf2*Sigma.
 * in the range of [0:minimum] the energy function is decreasing very slowly.
 * E(Distance) =
 *   0.0001*Epsilon*(Distance - sqrSix(2)*Sigma)*(Distance - sqrSix(2)*Sigma) - Epsilon  //Distance < sqrSix(2)*Sigma
 *   [4*sigma^6*Epsilon*(sigma^6/Distance^12 - 1/Distance^6)]*[(rMax-Distance)^2/(rMax-Distance)^2+Alpha]
 *   (a graph shown this function can be found in file: UnderZeroLJEnergyElement.gif)
 */
public class HydrogenBondsEnergyElement extends NonBondedEnergyElement {

    //---------------------------------- data filds -----------------------------------

    /*
     * weight given to this energy in the TotalEnergy element
     */
    private double weight;  
    private double energy;
    protected Atom oAtom, hAtom;
    private double deDxOAtom, deDyOAtom, deDzOAtom;
    private double deDxHAtom, deDyHAtom, deDzHAtom;
    private double dEdD, dEdX, dEdY, dEdZ;
    
    private int hFactor,oFactor;
    
    //--------- fields needed for energy calculationm -------

    public static final double ALPHA = 0.00001;
    private static final double ONE_DIV_SIX = 1.0/6.0;
    private static final double SQR_SIX_OF_TWO = Math.pow(2,ONE_DIV_SIX);
    
    /*
     * when dis > Rmax the energy quench to zero.
     */
    private double rMax;
    ///*
    // * ParameterList determine the epsilon and sigma according to the atoms type 
    // */
    //private LennardJonesParametersList parametersList;
    private double epsilon, sigma, sigma6, sigma6EpsilonFour, minusTwelveSigma6;
    private Atom atom1, atom2;
    /*
     * contact is used to make sure the function get to 0 in rMax
     */
    private double contact;
    private double dCdD;

    //true if no set was done or explisit call was done to freeElement
    private boolean free = true;
    /*
    * minimum of this energy function achived when distance = sqrSixOfTwoSigma
    */
    private double sqrSixOfTwoSigma;
    /*
     * An Element is built according to distance between two atoms
     */
    private Distance distance;
    /*
     * Attribute of distance contains:
     * is this relevant hydrogenBond;
     * its parameters;
     * its energy and derivative
     */
    public HB_DistanceAttribute hb_Attribute;
    //private LennardJonesParameters parameters;

    
    public HydrogenBondsEnergyElement() {
        rMax = DistanceMatrix.rMax();
    }

     public HydrogenBondsEnergyElement(double weight) {
        this.weight = weight;
        rMax = DistanceMatrix.rMax();
    }

    /*
    * @parm obj should be Distance
    */
    public void set(Object obj){
        free = false;
        distance = (Distance) obj;
        atoms = new AtomList(AtomPair.ATOM_PAIR_CAPACITY);
        atom1 = distance.atom1();
        atom2 = distance.atom2();
	atoms.add(atom1);
	atoms.add(atom2);
        HB_AtomAttribute atom1_attribute = (HB_AtomAttribute) atom1.getAttribute(HB_AtomAttribute.key);
        if (atom1_attribute.isH) {
            hFactor = oFactor = -1;    
            hAtom = atom1;
            oAtom = atom2;
        }
        else {
           hFactor = oFactor = 1;
            hAtom = atom2;
            oAtom = atom1;
        }
        //nAtom = hAtom.residue().amideN();    

        hb_Attribute = (HB_DistanceAttribute) distance.getAttribute(HB_DistanceAttribute.key);
        //update the parameters of the element:
        //parameters = hb_Attribute.parameters; 
        //parameters = (LennardJonesParameters) parametersList.parameters(distance);
        epsilon = 1;  //TODO hb_Attribute.epsilon;
        sigma = 1.7;//TODO hb_Attribute.sigma;
        sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;   //TODO //sigma6 = hb_Attribute.sigma6;
	    sigma6EpsilonFour = 4*epsilon*sigma6;   //TODO   // sigma6EpsilonFour = hb_Attribute.sigma6EpsilonFour;
	    minusTwelveSigma6 = -12*sigma6;     //TODO  // minusTwelveSigma6 = hb_Attribute.minusTwelveSigma6;
        sqrSixOfTwoSigma = SQR_SIX_OF_TWO*sigma;
       
        
    }    


    public double evaluate(double weight) {
        energy = updateEnergy();
        updateAtoms(weight);
        return energy * weight;
    }


    public double evaluate(){
        energy = updateEnergy();
        updateAtoms();
        return energy * weight;
     }


    /*
    * free all the global vars to avoid missuse
    */
    public void freeElement(){
        free = true;
        distance = null;
        hb_Attribute = null;
        atoms = null;
        oAtom = hAtom = atom1 =atom2 =null;
        contact = dCdD =dEdD =dEdX =deDxHAtom =deDxOAtom =dEdY = 0;
        deDyOAtom =deDyHAtom =dEdZ = deDzHAtom =deDzOAtom =0;
        energy =0;
        hFactor =oFactor = 0;
        epsilon = 0;
        sigma = 0;
        sigma6 = 0;
        sigma6EpsilonFour = 0;
        minusTwelveSigma6 = 0;
        sqrSixOfTwoSigma = 0;
    }
    
   /**
    * energy and dirivarives calculation.
    **/
   public double updateEnergy() {
       double dis = distance.distance();          	
       double rMaxMinusDis = rMax - dis;        
       double energy1, dE1dD;
       
       //smooth high energies: the function line goes to x=0 is almost a constant line 
       if (dis < sqrSixOfTwoSigma){
           
           energy1 = 0.0001*epsilon*(dis - sqrSixOfTwoSigma)*(dis - sqrSixOfTwoSigma) - epsilon;
           dE1dD = 0.0002*epsilon*(dis - sqrSixOfTwoSigma);
            if(! (energy1 >0 | energy1 <= 0) )
            System.out.print("eee");
       }    
       else {           
           double invD = distance.invDistance();		
           double invD2 = invD*invD;
           double invD6 = invD2*invD2*invD2;
           double invD7 = invD6*invD;
           double invD12 = invD6*invD6;
           double invD13 = invD12*invD;
            
           energy1 = sigma6EpsilonFour * (sigma6*invD12 - invD6);
           dE1dD = sigma6EpsilonFour * ( minusTwelveSigma6*invD13 + 6*invD7);
            if(! (energy1 >0 | energy1 <= 0) )
            System.out.print("eee");
        }
        if (energy1 > 0)  
            throw new RuntimeException("In hydrogenBondsEnergyElement: Energy should not get this values ! "+ energy1);
        
        //quench to zero in rMax
        double rMaxMinusDisSquare = rMaxMinusDis*rMaxMinusDis;
        double rMaxMinusDisSquarePlusAlpha = rMaxMinusDisSquare+ALPHA;
        double rMaxMinusDisTimesAlpha = rMaxMinusDis*ALPHA;
        double rMaxMinusDisSquarePlusAlphaSquare = rMaxMinusDisSquarePlusAlpha*rMaxMinusDisSquarePlusAlpha;

        contact = rMaxMinusDisSquare/rMaxMinusDisSquarePlusAlpha;
        dCdD = -2*rMaxMinusDisTimesAlpha/rMaxMinusDisSquarePlusAlphaSquare; //The minus comes from rMaxMinusDis
        double AnsEnergy = energy1 * contact;
        dEdD = dE1dD*contact + dCdD*energy1;
        
        dEdX = dEdD*distance.dDistanceDx();
        dEdY = dEdD*distance.dDistanceDy();
        dEdZ = dEdD*distance.dDistanceDz();

        deDxOAtom =    dEdX * oFactor;
        deDyOAtom =    dEdY * oFactor;
        deDzOAtom =    dEdZ * oFactor;
        deDxHAtom = -1*dEdX * hFactor;
        deDyHAtom = -1*dEdY * hFactor;
        deDzHAtom = -1*dEdZ * hFactor;

// following fields of the variable pair is need for the unit hydrogenBondsPairs
        hb_Attribute.set(deDxOAtom,
			      deDyOAtom,
			      deDzOAtom,
			      deDxHAtom,
			      deDyHAtom,
			      deDzHAtom,
			      AnsEnergy);

        if(! (AnsEnergy >0 | AnsEnergy <= 0) )
            System.out.print("eee");
         return AnsEnergy;
   }

    private void updateAtoms(){
        updateAtoms(weight);
    }

    private void updateAtoms(double weight){
       if (! oAtom.frozen()) {
           oAtom.addToFx(-1 * deDxOAtom * weight); // force = -derivative
           oAtom.addToFy(-1 * deDyOAtom * weight); // force = -derivative
           oAtom.addToFz(-1 * deDzOAtom * weight); // force = -derivative   
       }
       if (! hAtom.frozen()) {
           hAtom.addToFx(-1 * deDxHAtom * weight);
           hAtom.addToFy(-1 * deDyHAtom * weight);
           hAtom.addToFz(-1 * deDzHAtom * weight);
       }
// following fields of the variable pair is need for the unit hydrogenBondsPairs
        hb_Attribute.set (oAtom,
			      hAtom);
   }
	
   public final double deDxOAtom(){return deDxOAtom;}
   public final double deDyOAtom(){return deDyOAtom;}
   public final double deDzOAtom(){return deDzOAtom;}
   public final double deDxHAtom(){return deDxHAtom;}
   public final double deDyHAtom(){return deDyHAtom;}
   public final double deDzHAtom(){return deDzHAtom;}
    //    public UnderZeroLJEnergyElement ljElement(){return LJElement;} //  
   public final double energy(){return energy;}
   public final Atom oAtom() {return oAtom;}
   public final Atom hAtom() {return hAtom;}  
             
   public String toString() {
       if ((oAtom == null) & (hAtom == null)) return "HydrogenBondEnergyElement: - atoms not set yet";
       if ((oAtom == null) | (hAtom == null)) throw new RuntimeException("HydrogenBondEnergyElement: This is weird\n"+
                                                                         "hAtom = "+hAtom+"\n"+
                                                                         "oAtom = "+oAtom); 
       if (free) return "HydrogenBondsEnergyElement was free (no values)!";
       return "HydrogenBondsEnergyElement (h,o): "+hAtom.residueNumber()+
           " "+oAtom.residueNumber()+" dis: "+distanceValue()+"\n"+
						 hAtom+"\n"+oAtom;
       //return "HydrogenBondsEnergyElement";
   }

    public final double distanceValue() {
        return distance.distance();
    }
    protected void setAtoms(){
	throw new RuntimeException("setAtoms() may not be used by HydrogenBondsEnergyElement for "+
				   "efficiency.");
    }
    
}
