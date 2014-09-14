package meshi.energy.linearRG;
import java.util.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.*;

public class LinearRgEnergy extends CooperativeEnergyTerm{
    private double targetRG;
    private RgCalculator rgCalculator;
    public LinearRgEnergy() {}
    
    public LinearRgEnergy(AtomList atomList, 
			  double weight,
			  double targetRG, RgCalculator rgCalculator) {
	super(toArray(rgCalculator),atomList, null, null , weight);
	comment = "LinearRG";
	this.targetRG = targetRG;
	this.rgCalculator = rgCalculator;
    }
 
    public void evaluateAtoms() {
    	double e = evaluate();
    	for (int c=0 ; c<atomList.size() ; c++)
    		atomList.atomAt(c).addEnergy(e/atomList.size());
    }
    
    public double evaluate() {
	double ALPHA = 1;
	if (! on) return 0.0;
	double rg = rgCalculator.rg();
	double cmx = rgCalculator.cmx();
	double cmy = rgCalculator.cmy();
	double cmz = rgCalculator.cmz();
	double rgMtarget = rg - targetRG;
	double size = atomList.size();
	if (rgMtarget <= 0) return 0;
	double rgMtargetPalpha = rgMtarget+ALPHA;
	double factor = weight/(rg*size);
	double energy = weight*rgMtarget*rgMtarget/rgMtargetPalpha;
	double deDrg = (2*rgMtarget*rgMtargetPalpha-rgMtarget*rgMtarget)/(rgMtargetPalpha*rgMtargetPalpha);
	
	for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    if (! atom.frozen()) {
		atom.addToFx(-deDrg*factor*(atom.x()-cmx)); // Negating so that it is force
		atom.addToFy(-deDrg*factor*(atom.y()-cmy)); // Negating so that it is force
		atom.addToFz(-deDrg*factor*(atom.z()-cmz)); // Negating so that it is force
	    }
	}
	return energy;
    }
}
