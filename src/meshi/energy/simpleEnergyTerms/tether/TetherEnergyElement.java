package meshi.energy.simpleEnergyTerms.tether;
import meshi.energy.simpleEnergyTerms.*;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.energy.*;
import meshi.parameters.*;
import meshi.geometry.*;


public class TetherEnergyElement extends EnergyElement {
    protected Atom atom, evaluatedAtom;   
    protected double force;
    protected double evaluatedX;
    protected double evaluatedY;
    protected double evaluatedZ;
    //protected boolean frozen;
    protected FreeDistance distance;
    
    public TetherEnergyElement() {}
    
    public TetherEnergyElement(Atom inputAtom,double x,double y,double z, double weight) {
        evaluatedX=x;
        evaluatedY=y;
        evaluatedZ=z;
	
	atom = inputAtom;
	atoms = new AtomList();
	atoms.add(atom);
	MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
	new MolecularSystem();
	evaluatedAtom=  new Atom("eval", null, AtomType.XXX, new Coordinates(evaluatedX,evaluatedY,evaluatedZ),0.0);
	MolecularSystem.setCurrentMolecularSystem(saveMS);
	distance = new FreeDistance(atom,evaluatedAtom);
	force = weight;
	frozen = atom.frozen();
    }

    protected void setAtoms(){};

    public double evaluate() {
	double d;
	double deDd;
	double deDx;
	double deDy;
	double deDz;
	double energy;
	double dMinusEpsilon;
	
	if (frozen) return 0;	    
	distance.update();
	d = distance.distance();
	if (d ==0) return 0.0;

	energy = d * d * force;
	deDd =  d * force*2;



	deDx = deDd*distance.dDistanceDx();
	deDy = deDd*distance.dDistanceDy();
	deDz = deDd*distance.dDistanceDz();
	if (! atom.frozen()) {
	    atom.addToFx(-1*deDx); // force = -derivative
	    atom.addToFy(-1*deDy);
	    atom.addToFz(-1*deDz);
	}
	
	return energy;
    }


         /*
    private double distance(Atom atom, double X,double Y, double Z){
    
     return Math.sqrt((atom.x()-X)*(atom.x()-X)+(atom.y()-Y)*(atom.y()-Y)+(atom.z()-Z)*(atom.z()-Z));
    	
    }                                             */

	  public void update() {}

    public String toString() {
	return atom.toString();
    }

    



}
