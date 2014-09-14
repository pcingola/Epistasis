package meshi.energy.simpleEnergyTerms.bond;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.energy.*;

public class BondEnergyElement extends EnergyElement {
    protected Atom atom1, atom2;
    protected int number1, number2;
    protected Distance distance;
    protected double target, force, force2;
    private   boolean  on;
    double weight;
    public BondEnergyElement() {}
    public BondEnergyElement(AtomPair atomPair, Parameters parameters, 
			     DistanceMatrix distanceMatrix, double weight) {
	this.weight = weight;
	atom1 = atomPair.atom1();
	atom2 = atomPair.atom2();
	setAtoms();
	int atom1Number = atomPair.atom1Number();
	int atom2Number = atomPair.atom2Number();
	distance = distanceMatrix.distance(atom1Number, atom2Number);
	try {
	    if (distance == null) distance = new DistanceMirror(distanceMatrix.distance(atom2Number, atom1Number));
	}
	catch (RuntimeException ex) {
	    System.out.println("A problem while building a BondEnergyElement\n"+
			       "atom1 = "+atom1+"\n"+
			       "atom2 = "+atom2+"\n"+
			       "distanceMatrix.distance(atom1Number, atom2Number) = "+distanceMatrix.distance(atom1Number, atom2Number)+"\n"+
			       "distanceMatrix.distance(atom2Number, atom1Number) = "+distanceMatrix.distance(atom2Number, atom1Number)+"\n"+
			       "distanceMatrix.bondedList().size() = "+distanceMatrix.bondedList().size());
	    
	    for (Distance d:distanceMatrix.bondedList()) 
		System.out.println("bonded = "+d);
	    throw ex;
	}
			       
	target = ((BondParameters) parameters).target;
	force = ((BondParameters) parameters).force*weight;
	force2 = ((BondParameters) parameters).force2*weight;	    
	updateFrozen();
	on   = true;
    }
    
    public boolean updateFrozen() {
	super.updateFrozen();
	frozen = frozen || (distance == null);
	return frozen;
    }

    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
    }

    public double evaluate() {
	double d,d2,d3;
	double deDd;
	double deDx;
	double deDy;
	double deDz;
	double energy;
	double ALPHA = 0.0001;
	double d2PlusAlpha;
	double d2PlusAlpha2;
	double dis = distance.distance();
	if (frozen()) return 0.0;
	if (!on)      return 0.0;
	d = dis - target;
	d2 = d*d;
	energy = d2 * force;
    deDd =  d * force2;

	if (dis < 0.2) {
		energy += (dis-0.2)*(dis-0.2)*100000;
        deDd += 200000*(dis-0.2);
	}

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

    public String toString() {
	return ("BondEnergyElement target = "+dFormatSrt.f(target)+" force = "+dFormatSrt.f(force)+" distance = "+
		dFormatSrt.f(distance.distance())+" energy = "+dFormatSrt.f(evaluate())+"\n"+
		atom1.verbose(1)+"\n"+atom2.verbose(1));
    }

    public Atom atom1() {return atom1;}
    public Atom atom2() {return atom2;}
    public boolean isOn()    {return on;}
    public void    turnOn()  {on = true;}
    public void    turnOff() {on = false;}

}
