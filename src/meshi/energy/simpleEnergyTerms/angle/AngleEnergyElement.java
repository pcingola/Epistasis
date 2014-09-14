//-----------------------------------------------------------------------------------------------
package meshi.energy.simpleEnergyTerms.angle;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.*;
import meshi.molecularElements.*;   
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import java.util.*;
public class AngleEnergyElement extends EnergyElement {
    protected Atom atom1, atom2, atom3;
    protected int number1, number2, number3;
    protected Angle angle;	
    public double target, force, force2;
    protected final double infForce = 1.0/7.0; // About ~10 degree around PI and 0 will have affected forces
    protected final double infForce3 = infForce * infForce * infForce;	
    double weight;
    public AngleEnergyElement(Angle angle, Parameters parameters, double weight) {
	this.weight = weight;
	atom1 = angle.atom1;
	atom2 = angle.atom2;
	atom3 = angle.atom3;
	this.angle = angle;
	target = ((AngleParameters) parameters).target;
	if ((target<infForce) || (target + infForce > Math.PI)) 
	    throw new RuntimeException("The target angle " + target + " is too close to 0 or PI." +
				       "currently MESHI can not handle those values."); 	    
	force = ((AngleParameters) parameters).force*weight;
	force2 = ((AngleParameters) parameters).force2*weight;
	setAtoms();
	updateFrozen();
    }
	
    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
	atoms.add(atom3);
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
	angleValue = angle.angle();
	if ((angleValue>infForce) && (angleValue<Math.PI-infForce)) {
	    d = angleValue - target;
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
	    double aux1 = Math.PI - infForce - target;
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
	    double aux1 = infForce - target;
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
	
	
        public String toString() {
	return ("AngleEnergyElement target = "+dFormatSrt.f(Angle.rad2deg(target))+" force = "+dFormatSrt.f(force)+" angle = "+
		dFormatSrt.f(Angle.rad2deg(angle.angle()))+" energy = "+dFormatSrt.f(evaluate())+"\n"+
		atom1.verbose(1)+"\n"+atom2.verbose(1)+"\n"+atom3.verbose(1));
    }
}
