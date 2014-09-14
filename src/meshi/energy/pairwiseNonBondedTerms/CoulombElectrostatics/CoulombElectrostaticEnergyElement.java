package meshi.energy.pairwiseNonBondedTerms.CoulombElectrostatics;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;

/**
 *  This object represents Electrostatic Energy between two Atoms. 
 *  It updates the energy according to the change in their distance and position.
 **/

public  class CoulombElectrostaticEnergyElement extends NonBondedEnergyElement {
	
    public static final double MAX_ENERGY = 100;
    public static final double ALPHA = 0.5;
    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected int atom1Number, atom2Number;
    protected double dielectricConstant, q1, q2 ;//q1 = charge of atom # 1, etc.
    protected boolean frozen;
    protected double dEdD; // Derivation of the electrostatic energy according to the distance between two atoms.
    protected double dEdX; // Derivation of the electrostatic energy according to the X coordinate
    protected double dEdY; // Derivation of the electrostatic energy according to the Y coordinate
    protected double dEdZ;// Derivation of the electrostatic energy according to the Z coordinate
    protected double energy; 
    protected double weight; 
    protected double rMax; // The maximum distance between two atoms
    //that will be considered when calculating the electrostatics energy. 
    protected double contact, dCdD, dCdX, dCdY, dCdZ;
    private final int FIRST = 0,SECOND=1;
    protected double[] charges;

    /**
     * default constructor
     *
     **/
    public  CoulombElectrostaticEnergyElement() {}
	
    /**
     * constructor
     * 
     * @param distanceMatrix
     * @param weight
     * @param dielectricConstant
     * @param charges
     **/
    public  CoulombElectrostaticEnergyElement(DistanceMatrix distanceMatrix, 
                                              double weight, double dielectricConstant, double[] charges) {
        this.charges = charges;
        this.weight = weight;
        this.distanceMatrix = distanceMatrix;
        this.dielectricConstant = dielectricConstant;
        rMax = distanceMatrix.rMax();
    }
	
    /**
     * setAtoms
     **/
    protected void setAtoms(){
        throw new RuntimeException("setAtoms() may not be used by ElectrostaticEnergyElement for efficiency.");
    }
    
    // 
    /**
     * Sets the relevant charge values for each atom in an atom pair.
     * The data is taken from the parameter list.
     * @param obj  an AtomPair
     **/
    public void set(Object obj) {
        Distance nonBonded = (Distance)obj;
        atoms = new AtomList(2);
        atom1 = nonBonded.atom1();
        atom2 = nonBonded.atom2();
	atoms.add(atom1);
	atoms.add(atom2);
        atom1Number = atom1.number();
        atom2Number = atom2.number();
        q1 = charges[atom1.type().ordinal()];
	q2 = charges[atom2.type().ordinal()];
	if ((q1 == Double.MAX_VALUE) ||(q2 == Double.MAX_VALUE))
	    throw new RuntimeException("No parameters for "+atom1+" and "+atom2);
    }
   
    /**
     * evaluate -
     * 1) Invokes updateEnergy() - updates the energy,
     * 2) Invokes updateAtoms() - updates the atoms' position.
     * @return double - energy*weight
     **/	
    public double evaluate() {
        updateEnergy();
        updateAtoms();
        return energy*weight;
    }
    
    /**
     * Updates the energy
     * @return double - energy*weight
     **/	
    public double updateEnergy(){
        double EL;	
        double dELdD;
        double rMaxMinusDis, rMaxMinusDisPlus;
        double energy1, dE1dD;

        // if one of the atoms has zero charge, then the electostatic energy is 0.
        if(q1==0 || q2==0){
            energy = dEdD = dEdX = dEdY = dEdZ = contact = dCdX = dCdY = dCdZ = dELdD = 0;
        }

        else {
            //double EL;	//electrostatics
            //double dELdD;
            double invD = 0;
            double invD2;
	    double invD3;
            double multipleCharges = -99999; 
            double dis = -1;
			
            Distance distance;
			
            distance = distanceMatrix.distance(atom1Number, atom2Number);
            dis = distance.distance();
            rMaxMinusDis = rMax - dis;
            if (rMaxMinusDis <= 0) {
                energy = dEdD = dEdX = dEdY = dEdZ = contact = dCdX = dCdY = dCdZ = dELdD =  0;
                EL =0;
            }
            else {
                invD = distance.invDistance();		
                invD2 = invD*invD;
		invD3 = invD*invD2;
                multipleCharges = q1 * q2;

                EL =  multipleCharges * invD2 / dielectricConstant;
                dELdD = -2 * multipleCharges * invD3 / dielectricConstant ; 
		
            }

            energy1 = EL;
            dE1dD = dELdD;


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
				
        }
        return energy;
    }
   
   
    /**
     * Updates the atoms position
     *
     **/

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

    /**
     * @param weight
     * @return double
     **/
    public double contactAtoms(double weight) {
        double e = contact(weight);
        atom1.addEnergy(e/2);
        atom2.addEnergy(e/2);
        return e;
    }
	
    /**
     * @param index
     * @return double
     **/
    public double dEdXAtom(int index){
        if (index == FIRST)
            return dEdX;
        if (index == SECOND)
            return -1*dEdX;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }
	
    /**
     * 
     * @param index
     * @return double
     **/
    public double dEdYAtom(int index){
        if (index == FIRST)
            return dEdY;
        if (index == SECOND)
            return -1*dEdY;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }
	
    /**
     * 
     * @param index
     * @return double
     **/
    public double dEdZAtom(int index){
        if (index == FIRST)
            return dEdZ;
        if (index == SECOND)
            return -1*dEdZ;
        else
            throw new RuntimeException("index must be 0 or 1 !"); 
    }

    /**
     * 
     * @param weight
     * @return double
     **/
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
    
    /**
     * 
     * @return double
     **/
    public double contact() {
        return contact;
    }
    
    /**
     * toString
     * @return String
     **/
    public String toString() {
        if ((atom1 == null) & (atom2 == null)) return "ElectrostaticEnergyElement - atoms not set yet";
        if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
                                                                          "atom1 = "+atom1+"\n"+
                                                                          "atom2 = "+atom2); 
        Distance distance = distanceMatrix.distance(atom1Number, atom2Number);
        double dis = distance.distance();
        return ("ElectrostaticEnergyElement q1 = "+q1+" q2 = "+ q2 +" dielectricConstant = "+dielectricConstant+" Distance = "+
                dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
    }

}
