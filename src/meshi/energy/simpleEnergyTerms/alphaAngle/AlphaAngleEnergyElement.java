package meshi.energy.simpleEnergyTerms.alphaAngle;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.parameters.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import java.util.*;
//---------------------------------------------------------------------
public class AlphaAngleEnergyElement extends EnergyElement implements MeshiPotential {
    private Atom atom1,atom2,atom3;
    private Angle angle;
    private double loAng,hiAng;
    private double weight;
    private final double BUFF = 0.05;
    private final double INV_BUFF = 1/BUFF;
    private final double infForce = 1.0/7.0; // About ~10 degree around PI and 0 will have affected forces
    private final double infForce3 = infForce * infForce * infForce;	
    public double force, force2;
    private boolean  on;

    /**
     *The constructor needs a CA3 angle. 
     *It checks for the secondary structure the angle is in, and then assign the proper 
     *CA3 angle range.
     **/	
     
    public AlphaAngleEnergyElement(Angle angle,
					     AlphaAngleParameters param, double weight) {
	force = param.weightAA*weight;
	force2 = 2*param.weightAA*weight;
	this.angle = angle;
	if (! angle.getAngleName().equals("CA3"))
	    throw new RuntimeException("\nError: Only CA3 angles can be treated\n");	        
	Residue residue = angle.atom2.residue();
	if (residue.secondaryStructure().equals(SecondaryStructure.HELIX)) {
		loAng = param.startAlphaHELIX;
		hiAng = param.endAlphaHELIX;
	}
	else if (residue.secondaryStructure().equals(SecondaryStructure.SHEET)) {
		loAng = param.startAlphaSHEET;
		hiAng = param.endAlphaSHEET;
	}
	else if (residue.secondaryStructure().equals(SecondaryStructure.COIL) || 
			residue.secondaryStructure().equals(SecondaryStructure.HELIX_OR_COIL) || 
			residue.secondaryStructure().equals(SecondaryStructure.SHEET_OR_COIL)) {
		loAng = param.startAlphaCOIL;
		hiAng = param.endAlphaCOIL;
	}
	else
	    throw new RuntimeException("\nError: residue with unassigned secondary structure\n"+
				       residue+" : "+residue.secondaryStructure());
	if ((loAng<infForce) || (hiAng<infForce) || (loAng + infForce > Math.PI) || (hiAng + infForce > Math.PI)) 
	    throw new RuntimeException("The allowed region " + loAng + "-" + hiAng + " is too close to 0 or PI." +
				       "currently MESHI can not handle those values."); 	    

	setAtoms();
	updateFrozen();
	on = true;	
    }
	

    /**
     * Angle energy calculation and atom forces updating.
     **/
    public double evaluate() {
	double d;
	double deDd;
	double energy = 0;
	double angleValue;
	    
	if (frozen()) return 0;
	if (!on) return 0;
	angleValue = angle.angle();
	if ((angleValue>loAng) && (angleValue<hiAng))
	   return 0;
	if ((angleValue>infForce) && (angleValue<=loAng)) {
	    d = angleValue - loAng;
	    energy = d * d * force;	
	    deDd =  -1 * d * force2; // force = -derivative  
	    if (! atom1.frozen()) {
		atom1.addToFx(deDd*angle.dangleDx1());  
		atom1.addToFy(deDd*angle.dangleDy1());    
		atom1.addToFz(deDd*angle.dangleDz1());
	    }
	    if (! atom2.frozen()) {
		atom2.addToFx(deDd*angle.dangleDx2());   
		atom2.addToFy(deDd*angle.dangleDy2());    
		atom2.addToFz(deDd*angle.dangleDz2());    
	    }    
	    if (! atom3.frozen()) {
		atom3.addToFx(deDd*angle.dangleDx3()); 
		atom3.addToFy(deDd*angle.dangleDy3());    
		atom3.addToFz(deDd*angle.dangleDz3()); 
	    }
	}
	else if ((angleValue<Math.PI-infForce) && (angleValue>=hiAng)) {
	    d = angleValue - hiAng;
	    energy = d * d * force;	
	    deDd =  -1 * d * force2; // force = -derivative  
	    if (! atom1.frozen()) {
		atom1.addToFx(deDd*angle.dangleDx1());  
		atom1.addToFy(deDd*angle.dangleDy1());    
		atom1.addToFz(deDd*angle.dangleDz1());
	    }
	    if (! atom2.frozen()) {
		atom2.addToFx(deDd*angle.dangleDx2());   
		atom2.addToFy(deDd*angle.dangleDy2());    
		atom2.addToFz(deDd*angle.dangleDz2());    
	    }    
	    if (! atom3.frozen()) {
		atom3.addToFx(deDd*angle.dangleDx3()); 
		atom3.addToFy(deDd*angle.dangleDy3());    
		atom3.addToFz(deDd*angle.dangleDz3()); 
	    }
	}
	else if (angleValue>=Math.PI-infForce) {
	    double aux1 = Math.PI - infForce - hiAng;
	    double aux2 = force * aux1 * infForce3;
	    double aux3 = force * aux1 * (aux1 - infForce);
	    double aux4 = (angleValue - Math.PI);
	    energy = aux3 + aux2/(aux4 * aux4);
	    deDd = 2*aux2/(aux4 * aux4 * aux4);
	    if (! atom1.frozen()) {
		atom1.addToFx(deDd*angle.dangleDx1());  
		atom1.addToFy(deDd*angle.dangleDy1());    
		atom1.addToFz(deDd*angle.dangleDz1());
	    }
	    if (! atom2.frozen()) {
		atom2.addToFx(deDd*angle.dangleDx2());   
		atom2.addToFy(deDd*angle.dangleDy2());    
		atom2.addToFz(deDd*angle.dangleDz2());    
	    }    
	    if (! atom3.frozen()) {
		atom3.addToFx(deDd*angle.dangleDx3()); 
		atom3.addToFy(deDd*angle.dangleDy3());    
		atom3.addToFz(deDd*angle.dangleDz3()); 
	    }	    	
	}
	else if (angleValue<=infForce) {
	    double aux1 = infForce - loAng;
	    double aux2 = -1 * force * aux1 * infForce3;
	    double aux3 = force * aux1 * (aux1 + infForce);
	    double aux4 = angleValue;
	    energy = aux3 + aux2/(aux4 * aux4);
	    deDd = 2*aux2/(aux4 * aux4 * aux4);
	    if (! atom1.frozen()) {
		atom1.addToFx(deDd*angle.dangleDx1());  
		atom1.addToFy(deDd*angle.dangleDy1());    
		atom1.addToFz(deDd*angle.dangleDz1());
	    }
	    if (! atom2.frozen()) {
		atom2.addToFx(deDd*angle.dangleDx2());   
		atom2.addToFy(deDd*angle.dangleDy2());    
		atom2.addToFz(deDd*angle.dangleDz2());    
	    }    
	    if (! atom3.frozen()) {
		atom3.addToFx(deDd*angle.dangleDx3()); 
		atom3.addToFy(deDd*angle.dangleDy3());    
		atom3.addToFz(deDd*angle.dangleDz3()); 
	    }	    	
	}
	return energy;
    }	 


    public Angle getAngle() { return angle; }

    public String toString() {
    	String str = "\nOne Angles Deep Element. Made of:\n--------------------\nAngle1:";
    	str = str.concat(angle.toString());
    	str = str.concat("\n");
    	return str;
    }
    
    protected void setAtoms() {
    	atoms = new AtomList();
    	atom1 = angle.atom1;
    	atom2 = angle.atom2;
    	atom3 = angle.atom3;
    	atoms.add(angle.atom1);
    	atoms.add(angle.atom2);
    	atoms.add(angle.atom3);
    }	

    public Atom atom1() {return atom1;}
    public Atom atom2() {return atom2;}
    public Atom atom3() {return atom3;}
    public boolean isOn()    {return on;}
    public void    turnOn()  {on = true;}
    public void    turnOff() {on = false;}

}
