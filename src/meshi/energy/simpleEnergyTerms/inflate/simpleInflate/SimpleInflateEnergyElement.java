package meshi.energy.simpleEnergyTerms.inflate.simpleInflate;
import  meshi.energy.simpleEnergyTerms.inflate.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import java.util.*;
//-------------------------------------------------------------------------------------------------------
//                                     InflateElement
//-------------------------------------------------------------------------------------------------------
public  class SimpleInflateEnergyElement extends EnergyElement {
    protected Atom atom1, atom2;
    //    protected AtomPair atomPair;
    protected Atom atom1Copy, atom2Copy;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected FreeDistance distance;
    protected double target;
    public static final double ALPHA = 0.1;
     
    public  SimpleInflateEnergyElement(Distance atoms,  DistanceMatrix distanceMatrix, double weight) {
	atom1 = atoms.atom1();
	atom2 = atoms.atom2();
	setAtoms();
	updateFrozen();
 	distance = new FreeDistance(atom1, atom2);
	this.weight = MeshiProgram.randomNumberGenerator().nextDouble()* weight; 
        target = distanceMatrix.rMax()*3;
    }

    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
    }
	    	
	
    public double evaluate() {
	double dis;
	double d, d2, twoD, invDplus, invDplus2;
	double deDd = 0;
	double deDx;
	double deDy;
	double deDz;
	double energy = 0;

	if (frozen()) return 0;
	distance.update();
	dis = distance.distance();
	d = dis - target;
	energy = weight*d*d;
	deDd   = 2 * weight*d;
	deDx = deDd*distance.dDistanceDx();
	deDy = deDd*distance.dDistanceDy();
	deDz = deDd*distance.dDistanceDz();
	if (! atom1.frozen()) {
	    atom1.addToFx(-1*deDx); // force = -derivative   
	    atom1.addToFy(-1*deDy);
	    atom1.addToFz(-1*deDz);
	}
	if (! atom2.frozen()) {
	    atom2.addToFx(deDx);
	    atom2.addToFy(deDy);
	    atom2.addToFz(deDz);
	}
	return energy;
    }

    public void scaleWeight(double factor) {
	weight *= factor;
    }
}
