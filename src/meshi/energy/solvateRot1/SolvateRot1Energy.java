package meshi.energy.solvateRot1;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.util.mathTools.Sigma;
import meshi.util.mathTools.Spline1D;
import meshi.util.rotamericTools.Rot1Arrays;
import meshi.util.rotamericTools.RotamericTools;


public final class SolvateRot1Energy extends CooperativeEnergyTerm implements Rot1Arrays{

     
    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays.
     **/
     private double[] residueVals;
     private double[] splineVals;     
     private double[] dSplineVals;     
     private double[] AtomSumSigmC;
     private double[] dAtomEnergydEnvior;
     private double[] forceX;
     private double[] forceY;       
     private double[] forceZ;
       
    /** 
     * These fields are for general use in the class
     **/       
     private int[] lut; 	// The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.      
     private int atomListSize;
     private SolvateRot1ParametersList parameters; // The instance of the parameter list object.
     private Spline1D[] splines; // The splines array, that can associate each atom in the protein with the spline that suits its type.
     private double[][] means;
     private double[][] oneOverSTDs;
     private double maximalZ;
     
    public SolvateRot1Energy() {}
    
    public SolvateRot1Energy(AtomList atomList, 
                    DistanceMatrix dm,
				    SolvateRot1ParametersList parameters,
				    double[][] pp,
				    double maximalZ,
				    double weight) {
	super(toArray(dm),atomList, dm, parameters, weight);
	int c;
	int maxAtomNum=-1;
	comment = "Undefined Solvation";
	atomListSize = atomList.size();
	this.parameters = parameters;
	this.maximalZ = maximalZ;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	if (parameters.maxEnd > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than:" + parameters.maxEnd);
	calcMeansAndStd(pp);
	residueVals = new double[pp.length];

    // Creating the auxilary arrays
    splineVals = new double[atomListSize];
    dSplineVals = new double[atomListSize];
    AtomSumSigmC = new double[atomListSize];
    dAtomEnergydEnvior = new double[atomListSize];
    forceX = new double[atomListSize];
    forceY = new double[atomListSize];
    forceZ = new double[atomListSize];
    
    
	// Creating the lookup table for the atom numbers.
	// The table converts the atom internal number (field of Atom) to its index  
	// in the atom list given to the constructor.      
	for (c=0; c<atomListSize ; c++) {
	    if (atomList.atomAt(c).number() > maxAtomNum)
		maxAtomNum = atomList.atomAt(c).number();
	}
	lut = new int[maxAtomNum+1];
	for (c=0; c<maxAtomNum ; c++) {
	    lut[c] = -1;
	}
	for (c=0; c<atomListSize ; c++) {
	    lut[atomList.atomAt(c).number()] = c;
	}
	
	// Setting up the splines array, that can associate each atom in the protein with the spline that suits its type.
	splines = new Spline1D[atomListSize];
	for (c=0; c<atomListSize ; c++) 
	    splines[c] = parameters.atomTypeSplines[atomList.atomAt(c).type().ordinal()];
										
    } // of the constructor
    
    public void setComment(String str) {
    	comment = str;
    }        
 
    public void evaluateAtoms() {
		evaluate(true);
    }
    
    public double evaluate() {
    	return evaluate(false);
    }

    public final double evaluate(boolean updateAtoms) {
	if (! on) return 0.0;
	double energy = 0;
	double atomEnergy=0; 
	int cc;        
	DistanceList dislist = dm.nonBondedList();
	Iterator iter;
	Distance dis;
	SolvateRot1DistanceAttribute sigmaValues; 
	Atom atom;
	int ind1,ind2;


	//Reseting the auxilary arrays and variables
	for (cc=0 ; cc<atomListSize ; cc++) {
	    AtomSumSigmC[cc] = 0; 
	    forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
	} 
	for (cc=0 ; cc<residueVals.length ; cc++) {
	    residueVals[cc] = 0.0; 
	} 
	
	// First pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
    if (dis.getAttribute(SolvateRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE) == null) {
      sigmaValues = new SolvateRot1DistanceAttribute();
      dis.addAttribute(sigmaValues);
      if (dis.mode().frozen)
      	updateSigmVals(dis);
    }
    if (!dis.mode().frozen) {
    	updateSigmVals(dis);
    	sigmaValues = (SolvateRot1DistanceAttribute) dis.getAttribute(SolvateRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
    	ind1 = lut[dis.atom1().number()];
    	ind2 = lut[dis.atom2().number()];
    	AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
    	AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
    }
	}

	//Calculating the energy values. Looping on all the atoms in the protein... TWICE!
	for (cc=0 ; cc<atomListSize ; cc++) {
		if (!atomList.atomAt(cc).frozen()) {
			splines[cc].calc(AtomSumSigmC[cc]);
			splineVals[cc] = splines[cc].s;
			dSplineVals[cc] = splines[cc].s_tag;
			residueVals[atomList.atomAt(cc).residueNumber()] += splines[cc].s;
		}
	}
	// A second time	
	for (cc=0 ; cc<atomListSize ; cc++) {
		atom = atomList.atomAt(cc);
		  if (!atom.frozen()) { // This line was deleted on 19.6.2006. It looks like a bug...
		if (means[atom.residueNumber()]!=null) {
			if (atom.name.equals("CA")) {
				energy += weight*
Math.min((residueVals[atom.residueNumber()] -
means[atom.residueNumber()][0])*
oneOverSTDs[atom.residueNumber()][0] , maximalZ);
				if (updateAtoms)
					atom.addEnergy(weight*
Math.min((residueVals[atom.residueNumber()] -
means[atom.residueNumber()][0])*
oneOverSTDs[atom.residueNumber()][0] , maximalZ));
			}
			dAtomEnergydEnvior[cc] = weight*oneOverSTDs[atom.residueNumber()][0]*dSplineVals[cc];
		}
		} // This line was deleted on 19.6.2006. It looks like a bug... 
	}
	   	  
	// Second pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
    if (!dis.mode().frozen) {
      sigmaValues = (SolvateRot1DistanceAttribute) dis.getAttribute(SolvateRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()]; 
	          		   
	   // Doing the self derivatives
	   forceX[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dx1;
	   forceY[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dy1;
	   forceZ[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dz1;
	   
	   forceX[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dy2;
	   forceZ[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dz2;

	   // Doing the cross derivatives
	   forceX[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dy2;
	   forceZ[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dz2;
	                  
	   forceX[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dx1;
	   forceY[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dy1;
	   forceZ[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dz1;
	}
	}  // second pass on the non bonded list                	
    
	// Finally, the appropriate forces are assigned for every atom. 
	for (cc=0 ; cc<atomListSize ; cc++) {
		atom = atomList.atomAt(cc);
		if (!atom.frozen()) {
		    atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
		    atom.addToFy(-forceY[cc]); // Negating so that it is realy force
		    atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
		}
	}
    return energy;
    }
 
 
    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and 
     * atom2 in the Distance - dis. The results are updated in the fields of the 
     * SolvateDistanceAttribute of dis - sigmaValues.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	SolvateRot1DistanceAttribute sigmaValues = 
    		(SolvateRot1DistanceAttribute) dis.getAttribute(SolvateRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE); 
    	int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
       	
    	sigmaValues.resetAllSigmVals();

    	// Hydrogens are not treated currently
    	// -------------------------------
    	if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen()) {     	
    		return;
    	}
    	    	

    	// Calculating the carbon sigmoid of atom1. 
    	// ---------------------------------------
    	if (dis.atom2().type().isCarbon()) {
   			Sigma.sigma(dis.distance(),
   						parameters.Cend[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.Cp1[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.Cp2[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.CvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.CvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
    		sigmaValues.sigmCa1 = Sigma.s;
    		sigmaValues.dsigmCa1dx1 = Sigma.s_tag*dis.dDistanceDx();
  		  	sigmaValues.dsigmCa1dy1 = Sigma.s_tag*dis.dDistanceDy();
  		  	sigmaValues.dsigmCa1dz1 = Sigma.s_tag*dis.dDistanceDz();
 		   	sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
 		   	sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
   		 	sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;
   		}

    	// Calculating the carbon sigmoid of atom2. 
    	// ---------------------------------------
    	if (dis.atom1().type().isCarbon()) {	
   			Sigma.sigma(dis.distance(),
   						parameters.Cend[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.Cp1[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.Cp2[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.CvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.CvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
    		sigmaValues.sigmCa2 = Sigma.s;
    		sigmaValues.dsigmCa2dx1 = Sigma.s_tag*dis.dDistanceDx();
  		  	sigmaValues.dsigmCa2dy1 = Sigma.s_tag*dis.dDistanceDy();
  		  	sigmaValues.dsigmCa2dz1 = Sigma.s_tag*dis.dDistanceDz();
 		   	sigmaValues.dsigmCa2dx2 = -sigmaValues.dsigmCa2dx1;
 		   	sigmaValues.dsigmCa2dy2 = -sigmaValues.dsigmCa2dy1;
   		 	sigmaValues.dsigmCa2dz2 = -sigmaValues.dsigmCa2dz1;
   		}
    	
   	} // of updateSigmVals


	public void calcMeansAndStd(double[][] pp) {
		means = new double[pp.length][];
		oneOverSTDs = new double[pp.length][];
		for (int c=0 ; c<pp.length ; c++) {
			if ((pp[c]!=null) && (((int) pp[c][2]) != 0) && (((int) pp[c][2]) != 5)) {
		    	double[] tmp = RotamericTools.getMean(((int) pp[c][2]), pp[c][3] ,0);
				means[c] = new double[1];
				oneOverSTDs[c] = new double[1];
				means[c][0] = tmp[0];
				oneOverSTDs[c][0] = 1/(tmp[1]+0.00000000001); 
			}
		}
	}
    
}
	
