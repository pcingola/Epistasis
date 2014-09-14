package meshi.energy.simpleEnergyTerms.alphaTorsion;
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
public class AlphaTorsionEnergyElement extends EnergyElement implements MeshiPotential {
    private Atom atom1,atom2,atom3,atom4;
    private Torsion torsion;
    private double loTor,hiTor,midTor,cyclicMidTor;
    private double weight;
    private final double BUFF = 0.05; // the size of the smoothing region
    private final double INV_BUFF = 1/BUFF; 
    private boolean  on;

    /**
     *An energy element handling one alpha torsion. 
     *It checks for the secondary structure the torsion is in, and then assign the proper 
     *CA torsion range.
     **/	
     
    public AlphaTorsionEnergyElement(Torsion torsion,
					     AlphaTorsionParameters param, double weight) {
	this.weight = param.weightAA*weight;
	this.torsion = torsion;
	if (! torsion.getTorsionName().equals("ALPHA"))
	    throw new RuntimeException("\nError: Only CA torsions can be treated\n");	        
	Residue residue = torsion.atom2.residue();
	if (residue.secondaryStructure().equals(HELIX)) {
		loTor = param.startAlphaHELIX;
		hiTor = param.endAlphaHELIX;
	}
	else if (residue.secondaryStructure().equals(SHEET)) {
		loTor = param.startAlphaSHEET;
		hiTor = param.endAlphaSHEET;
	}
	else
	    throw new RuntimeException("\nError: CA torsion is valid only to HELIX or SHEET secondary structure\n");
	midTor = (hiTor+loTor)/2;
	cyclicMidTor = (hiTor-2*Math.PI+loTor)/2;
	setAtoms();
	updateFrozen();
	on = true;
    }
    
    
    private int calc(double tor, double loTor, double hiTor, double[] ans) {
	   ans[0] = 0.0;
	   ans[1] = 0.0; 	
	   
	   if (loTor>hiTor) {
	   	  if ((tor<hiTor) | (tor>loTor))
	   	     return 0;
	   	  if ((loTor-hiTor)<4*BUFF)
	   	     return 0;
	   	  if ((tor>=hiTor) && (tor<hiTor+BUFF)) {
	   	  	ans[0] = 0.5*INV_BUFF*(tor-hiTor)*(tor-hiTor);
	   	  	ans[1] = INV_BUFF*(tor-hiTor);
	   	  }
	   	  else if ((tor>=hiTor+BUFF) && (tor<midTor-BUFF)) {
	   	  	ans[0] = 0.5*BUFF + (tor - hiTor - BUFF);
	   	  	ans[1] = 1.0;
	   	  }
	   	  else if ((tor>=midTor-BUFF) && (tor<midTor+BUFF)) {
	   	  	ans[0] = -0.5*INV_BUFF*(tor-midTor)*(tor-midTor) + (midTor - hiTor - BUFF);
	   	  	ans[1] = -INV_BUFF*(tor-midTor);
	   	  }
	   	  else if ((tor>=midTor+BUFF) && (tor<loTor-BUFF)) {
	   	  	ans[0] = 0.5*BUFF - (tor - loTor + BUFF);
	   	  	ans[1] = -1.0;	   	  	
	   	  }
	   	  else if (tor>=loTor-BUFF) {
	   	  	ans[0] = 0.5*INV_BUFF*(tor-loTor)*(tor-loTor);
	   	  	ans[1] = INV_BUFF*(tor-loTor);	   	  	
	   	  }
	   	  return 0;
	   }
	   else {
	   	  if ((tor<hiTor) && (tor>loTor))
	   	     return 0;
          if (tor>hiTor)
	   	     tor -= 2*Math.PI;
	   	  hiTor -= 2*Math.PI;
	   	  if ((loTor-hiTor)<4*BUFF)
	   	     return 0;
	   	  if ((tor>=hiTor) && (tor<hiTor+BUFF)) {
	   	  	ans[0] = 0.5*INV_BUFF*(tor-hiTor)*(tor-hiTor);
	   	  	ans[1] = INV_BUFF*(tor-hiTor);
	   	  }
	   	  else if ((tor>=hiTor+BUFF) && (tor<cyclicMidTor-BUFF)) {
	   	  	ans[0] = 0.5*BUFF + (tor - hiTor - BUFF);
	   	  	ans[1] = 1.0;
	   	  }
	   	  else if ((tor>=cyclicMidTor-BUFF) && (tor<cyclicMidTor+BUFF)) {
	   	  	ans[0] = -0.5*INV_BUFF*(tor-cyclicMidTor)*(tor-cyclicMidTor) + (cyclicMidTor - hiTor - BUFF);
	   	  	ans[1] = -INV_BUFF*(tor-cyclicMidTor);
	   	  }
	   	  else if ((tor>=cyclicMidTor+BUFF) && (tor<loTor-BUFF)) {
	   	  	ans[0] = 0.5*BUFF - (tor - loTor + BUFF);
	   	  	ans[1] = -1.0;	   	  	
	   	  }
	   	  else if (tor>=loTor-BUFF) {
	   	  	ans[0] = 0.5*INV_BUFF*(tor-loTor)*(tor-loTor);
	   	  	ans[1] = INV_BUFF*(tor-loTor);	   	  	
	   	  }	   	
	   	  return 0;
	   }
	}

    public double weight() { return weight;}
    public double evaluate() {
	double tor = torsion.torsion();
	double[] ans = new double[2]; 
	double energy, dEdTorsion1;
    	
	if (frozen()) return 0.0;
	if (!on)      return 0.0;
	
	calc(tor,loTor,hiTor,ans);
	energy = weight*ans[0];
	dEdTorsion1 = -weight*ans[1];  // minus - so it is force
    

	if (atom1.frozen() == false) {
	    atom1.addToFx(dEdTorsion1*torsion.dTorsionDx1());
	    atom1.addToFy(dEdTorsion1*torsion.dTorsionDy1());
	    atom1.addToFz(dEdTorsion1*torsion.dTorsionDz1());
	}
	if (atom2.frozen() == false) {
	    atom2.addToFx(dEdTorsion1*torsion.dTorsionDx2());
	    atom2.addToFy(dEdTorsion1*torsion.dTorsionDy2());
	    atom2.addToFz(dEdTorsion1*torsion.dTorsionDz2());
	}
	if (atom3.frozen() == false) {
	    atom3.addToFx(dEdTorsion1*torsion.dTorsionDx3());
	    atom3.addToFy(dEdTorsion1*torsion.dTorsionDy3());
	    atom3.addToFz(dEdTorsion1*torsion.dTorsionDz3());
	}	        	    
	if (atom4.frozen() == false) {
	    atom4.addToFx(dEdTorsion1*torsion.dTorsionDx4());
	    atom4.addToFy(dEdTorsion1*torsion.dTorsionDy4());
	    atom4.addToFz(dEdTorsion1*torsion.dTorsionDz4());
	}	        	    
    	return energy;
    }


    public Torsion getTorsion() { return torsion; }

    public String toString() {
    	String str = "\nOne Torsions Deep Element. Made of:\n--------------------\nTorsion:";
    	str = str.concat(torsion.toString());
    	str = str.concat("\n");
    	return str;
    }
    
    protected void setAtoms() {
    	atoms = new AtomList();
    	atom1 = torsion.atom1;
    	atom2 = torsion.atom2;
    	atom3 = torsion.atom3;
    	atom4 = torsion.atom4;
    	atoms.add(torsion.atom1);
    	atoms.add(torsion.atom2);
    	atoms.add(torsion.atom3);
    	atoms.add(torsion.atom4);
    }	

    public Atom atom1() {return atom1;}
    public Atom atom2() {return atom2;}
    public Atom atom3() {return atom3;}
    public Atom atom4() {return atom4;}
    public boolean isOn()    {return on;}
    public void    turnOn()  {on = true;}
    public void    turnOff() {on = false;}

}
