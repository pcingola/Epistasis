package meshi.energy.solvate;
import java.util.Iterator;
import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBond;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBondList;
import meshi.molecularElements.hydrogenBonds.HydrogenBondDahiyatList;
import meshi.parameters.*;
import meshi.util.mathTools.Sigma;
import meshi.util.mathTools.Spline1D;
import meshi.util.mathTools.Spline2D;

/**
 * The implementation of the cooperative solvation term for proteins as described in Kalisman & Keasar (2008). 
 * Since the cooperative solvation is described in the above paper, we bring here only the implementaion
 * details. Especially regarding the calculation of the derivatives, which was too lengthy for the paper. 
 * 
 * The class allows to give a different weight to the solvation energies of certain atom types in the final
 * summation described in Eq. 5. These atom types are side-chain polars, side-chain carbons, and backbone 
 * polars. The class also include a regular hydrogen bond term. The functional form is therefore: 
 *
 * Esolv = weightSCPolarSolvate*Eside_chain_polars + 
 * 		weightSCCarbonSolvate*Eside_chain_carbons + 
 *      weightBBPolarSolvate*Ebackbone_polars + 
 *      weightHB*Ehb
 *      
 * Where Ehb is the negative of the HBC summation over all atoms. The weights are defined in the "Creator" class, 
 * and passed to the constructor as parameters. 
 *
 * General remarks:
 * ----------------
 * 1) This term is derivable twice (because the splines are derivable twice).
 * 2) This term is using hydrogen bond description that is dependent on two angles in the bond. This decription 
 *    follows that of McDonald and Thornton (1994) and Dahiyat et al. (1996). The only place where the hydrogen bond
 *    list is declared explicitly is in line 204. This means that any hydrogen bond implementation that extends the 
 *    "AbstractHydrogenBondList" template can be used, by correcting line 204. 
 * 3) We calculate the regular hydrogen bond energy term (Ehb) together with the solvation terms themselves, since 
 *    the hydrogen bonds calculation is a by-product of the first step in the solvation evaluation, and is thus for free.
 * 4) Disulfide bonds are treated as "hydrogen bonds" between the SG's of two cystines.
 * 5) The SD sulfor of methionine is treated as a hydrophobic carbon.
 * 6) See the remarks in the "Creators" classes for a quick start on how to create a working instance of the term.
 *  
 *
 * The energy evaluation:
 * ----------------------
 * The energy value and derivatives is calculated in 3 steps:
 * 1) A first pass over the non-bonded list. Each Distance instance in the non-bonded-list, is used to update the
 * CNC's and HBC's of its atom pairs (Eqs. 1 and 2, respectively). The partial derivatives of the CNC's and HBC's
 * with respect to the distance atoms are alsoclaculated.  Since some of this values will be needed also in step 3, 
 * we save them in an instance of "SolvateDistanceAttribute" that is attached as an "Attribute" to the Distance instance.
 * 2) A pass on the atom list. Once we have the CNC and HBC of every atom in the protein, we can proceed to calculate the solvation energy 
 * associated with every atom. In this implementation we combined the EI(CNC,HBC) evaluation (Eq. 3) of every atom and 
 * the -log(spline(EI)) evaluation (Eq. 4) into a single step by using a 2D spline, i.e. spline2D(CNC,HBC). The 2D spline 
 * is, of course, atom type specific. The derivatives of each atom solvate energy value with respect to the HBC and CNC  
 * are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect 
 * to the CNC's and HBC's) from step 2, we can now calculate the energy derivative with respect to the atomic 
 * coordinates. In this step we simply make sure that every term that arises from the derivative chain rule is accounted for. 
 *
 **/
public final class SolvateEnergy extends CooperativeEnergyTerm {
    
    // Relative strength of SALT BRIDGES compared with regular HYDROGEN BONDS for desolvation purposes, i.e. in
    // regard to the effect on observed CNC medians. Following Table 1 in the paper. 
    public static final double SALT_BRIDGE_STRENGTH_ASP_OD = 0.45; 	 
    public static final double SALT_BRIDGE_STRENGTH_GLU_OE = 0.45; 	 
    public static final double SALT_BRIDGE_STRENGTH_LYS_NZ = 1.65; 	 
    public static final double SALT_BRIDGE_STRENGTH_ARG_NH = 1.5; 	 
    public static final double SALT_BRIDGE_STRENGTH_TRO = 0.45; 	 
    public static final double SALT_BRIDGE_STRENGTH_TRN = 1.8; 	 
    /** The following parameter allow for a different weighting of SALT BRIDGES compared with regular HYDROGEN BONDS for the 
	Ehb energy, that is also claculated. 
    **/
    public static final double SALT_BRIDGE_STRENGTH_GENERAL = 1.0; 	 

    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays. The lengths of these arrays is the length of the atom list. 
     * The indexing to these arrays is by the index of the atom in the atom list.
     **/
    private double[] CNC;  // Eq. 1
    private double[] HBC;  // Eq. 2
    private double[] HBCforHBenergy; // May be different because we allow different weighting of the salt bridges in the regular HB energy.
    private double[] dSplineDCNC;     
    private double[] dSplineDHBC;     
    private double[] forceX;
    private double[] forceY;       
    private double[] forceZ;
       
    /** 
     * These following fields are for general use in the class
     **/       
    /** Size of the atom list. **/
    /**  The instance of the parameter list object. **/
    private SolvateParametersList parameters; 
    /** The look-up table (lut) converts the atom internal number (field of Atom), which is the index of the array, to its 
     * index in the atom list given to the constructor. **/
    private int[] lut; 	  
    /** Setting the general type for each atom in the atom list: (0) Carbon (1) Backbone polar, (2) Sidechain polar (3) Hydrogens 
     * This is done to save time on type checking. The index to the array is the index of an atom in the atom list. **/
    private int[] superType;  	
    /** The 2D spline array (i.e. spline(CNC,HBC)) for the polar side-chain atoms. The indexing in the array is the 
     * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) **/
    private Spline2D[] splinesSCPolar; 
    /** The 2D spline array (i.e. spline(CNC,HBC)) for the polar backbone atoms. The indexing in the array is the 
     * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) **/
    private Spline2D[] splinesBB; 
    /** The 1D spline array (i.e. spline(CNC,HBC)) for the polar carbon atoms. The indexing in the array is the index of 
     * the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded). The splines
     * are 1D because HBC is 0, for carbons. **/
    private Spline1D[] splinesSCCarbon; 
    /** The only hydrogen bond list class in the term. **/ 
    private AbstractHydrogenBondList solvateHB;

    // Weights
    private double weightSCPolarSolvate;
    private double weightBBPolarSolvate;
    private double weightSCCarbonSolvate;
    private double weightHB;
    public final MolecularSystem molecularSystem;
    public final int molecularSystemSize;
    private static Atom atom;
    private static int atomTypeNumber;
    private static AtomType atomType;
    public SolvateEnergy() {molecularSystem = null; molecularSystemSize = -1;}
    
    /** 
     * See the comment at the top of the class for descriptions on the weights.
     **/
    public SolvateEnergy(AtomList atomList, 
			 DistanceMatrix dm,
			 SolvateParametersList parameters,
			 double weightSCPolarSolvate,
			 double weightBBPolarSolvate,
			 double weightSCCarbonSolvate,
			 double weightHB) {
	this(atomList, dm, parameters, weightSCPolarSolvate, weightBBPolarSolvate,
	     weightSCCarbonSolvate, weightHB,"Solvation");
    }
    public SolvateEnergy(AtomList atomList, 
			 DistanceMatrix dm,
			 SolvateParametersList parameters,
			 double weightSCPolarSolvate,
			 double weightBBPolarSolvate,
			 double weightSCCarbonSolvate,
			 double weightHB,String comment) {
	super(toArray(), atomList, dm, parameters, weightSCPolarSolvate);

	molecularSystem = atomList.molecularSystem();

	this.weightSCPolarSolvate = weightSCPolarSolvate;
	this.weightBBPolarSolvate = weightBBPolarSolvate;
	this.weightSCCarbonSolvate = weightSCCarbonSolvate;
	this.weightHB = weightHB;
	this.comment  = comment;
	int c;
	
	molecularSystemSize = molecularSystem.size();
	this.parameters = parameters;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	if (parameters.maxEnd > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than:" + parameters.maxEnd);
				       


	// Creating the auxilary arrays
	CNC = new double[molecularSystemSize];
	HBC = new double[molecularSystemSize];
	HBCforHBenergy = new double[molecularSystemSize];
	dSplineDCNC = new double[molecularSystemSize];
	dSplineDHBC = new double[molecularSystemSize];
	forceX = new double[molecularSystemSize];
	forceY = new double[molecularSystemSize];
	forceZ = new double[molecularSystemSize];
	superType = new int[molecularSystemSize];
    
	// Setting up the splines array, that can associate each residue in the protein with the spline that 
	// suits its type.
	splinesSCPolar = new Spline2D[molecularSystemSize];
	splinesSCCarbon = new Spline1D[molecularSystemSize];
	splinesBB = new Spline2D[molecularSystemSize];
	for (c=0 ; c<molecularSystemSize ; c++)  {
	    atom = molecularSystem.get(c).atom;
	    if (! atom.nowhere()) {
		atomTypeNumber = atom.type().ordinal();
		splinesSCPolar[c] = parameters.scPolarSplines[atomTypeNumber];
		splinesSCCarbon[c] = parameters.scCarbonSplines[atomTypeNumber];
		splinesBB[c] = parameters.bbSplines[atomTypeNumber];
	    }
	}
		
	// Setting up the HB list
	solvateHB = new HydrogenBondDahiyatList(dm, atomList, parameters);
	updateableResources.add(solvateHB);
	
	// Determining the superType for each atom: (0) Hydrophobic, (1) Backbone polar, (2) Sidechain polar (3) Hydrogens
	for (c=0 ; c<molecularSystemSize ; c++)  {
	    atom = molecularSystem.get(c).atom;
	    atomType = atom.type();
	    if (atomType.isCarbon() || (atomType == AtomType.MSD))
		superType[c] = 0;
	    else if ((atomType.isOxygen() || atomType.isNitrogen()) && (atom.name().length()==1))
		superType[c] = 1;
	    else if ((atomType.isOxygen() || atomType.isNitrogen() || (atomType == AtomType.CSG)) 
		     && (atom.name().length()>1))
		superType[c] = 2;
	    else if (atomType.isHydrogen())
		superType[c] = 3;
	    else
		throw new RuntimeException("The following if's should have covered all cases. But due to my carelessness this atom was left out:\n" +
					   atom + "\n");
	}
    }  // Of constructor
    
    public void evaluateAtoms() {
	evaluate(true,weightSCPolarSolvate,weightSCCarbonSolvate,weightBBPolarSolvate,weightHB);
    }
    
    /**
     * Calculates Esolv with the weights given in the constructor.
     **/
    public double evaluate() {
    	return evaluate(false,weightSCPolarSolvate,weightSCCarbonSolvate,weightBBPolarSolvate,weightHB);
    }

    /**
     * Calculates Esolv with the weights you give as parameters!
     **/
    public final double evaluate(boolean updateAtoms, double W_SCPolarSolvate, double W_SCCarbonSolvate, double W_BBPolarSolvate,
				 double W_HB) {
	if (! on) return 0.0;
	double energy = 0;
	int cc;        
	Iterator iter;
	int ind1,ind2;

	resetAuxilaryArrays();
	firstPassOverTheNonBondedList();
	energy = calculatingTheSolvationEnergiesOfEachAtom(updateAtoms, W_SCPolarSolvate, W_SCCarbonSolvate, W_BBPolarSolvate,
							   W_HB);
	secondPassOverTheNonBondedList(W_HB);
	assignForcesToAtoms();
    	return energy;
    }
 
    private void resetAuxilaryArrays() {
	int cc; 
	//Reseting the auxilary arrays
	for (cc=0 ; cc<molecularSystem.size() ; cc++) {
	    CNC[cc] = 0;
	    HBC[cc] = 0; 
	    HBCforHBenergy[cc] = 0; 
	    dSplineDCNC[cc] = 0;
	    dSplineDHBC[cc] = 0;
	    forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
	} 
    }

    public void firstPassOverTheNonBondedList(){
	int ind1,ind2;
	DistanceList dislist = dm.nonBondedList();
	SolvateDistanceAttribute sigmaValues=null;
	SolvateDistanceAttributeBetweenPolars sigmaValuesBP = null;
	SolvateDistanceAttributeWithNonPolar  sigmaValuesWNP = null;
	// ***********************************
	// First pass over the non-bonded list
	// ***********************************
	for (Distance dis:dislist) {
	    sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
	    // In case there is not a "SolvateDistanceAttribute" in the distance, we create one.
	    if (sigmaValues == null) {
		sigmaValues = SolvateDistanceAttributeCreator.create(dis, parameters,solvateHB);
		dis.addAttribute(sigmaValues);
	    }
	    if (!sigmaValues.isDisInvolvesHydrogens) {
		ind1 = dis.atom1.number;
		ind2 = dis.atom2.number;
		if (sigmaValues.isDisBetween2Polars) {
		    sigmaValuesBP = (SolvateDistanceAttributeBetweenPolars) sigmaValues;
		    sigmaValuesBP.update();
		    HBC[ind1] += sigmaValuesBP.saltBridgeFactorA1 * sigmaValuesBP.sigmHB;
		    HBCforHBenergy[ind1] += sigmaValuesBP.saltBridgeFactorForHBenergyA1 * sigmaValuesBP.sigmHB;
		    HBC[ind2] += sigmaValuesBP.saltBridgeFactorA2 * sigmaValuesBP.sigmHB;
		    HBCforHBenergy[ind2] += sigmaValuesBP.saltBridgeFactorForHBenergyA2 * sigmaValuesBP.sigmHB;
		}
		else {
		    sigmaValuesWNP = (SolvateDistanceAttributeWithNonPolar) sigmaValues;
		    sigmaValuesWNP.update();
		    CNC[ind1] += sigmaValuesWNP.sigmCa1;
		    CNC[ind2] += sigmaValuesWNP.sigmCa2;
		}
	    }
	}
    }

    public double calculatingTheSolvationEnergiesOfEachAtom(boolean updateAtoms, double W_SCPolarSolvate, double W_SCCarbonSolvate, double W_BBPolarSolvate,
							    double W_HB){
	int cc;
	double energy = 0;
	double energySCPolar=0; 
	double energySCCarbon=0; 
	double energyBBON=0; 
	// ***************************************************************************
	// A pass over the atom list. Calculating the solvation energies of each atom.
	// ***************************************************************************
	for (cc=0 ; cc<molecularSystem.size() ; cc++) {
	    if (! molecularSystem.get(cc).atom.nowhere()) {
		energySCPolar=energySCCarbon=energyBBON=0.0; 
		if (superType[cc]==0) {		// non-polar atom 
		    try {
			splinesSCCarbon[cc].calc(CNC[cc]);
		    }
		    catch (Exception e) {
			throw new RuntimeException("SCC: " + CNC[cc] + " " + HBC[cc] + "\n" +
						   molecularSystem.get(cc).atom + "\n");
		    }
		    energySCCarbon = W_SCCarbonSolvate*splinesSCCarbon[cc].s;
		    dSplineDCNC[cc] = W_SCCarbonSolvate*splinesSCCarbon[cc].s_tag;
		    dSplineDHBC[cc] = 0.0;
		}
		else if (superType[cc]==1) { 		// O and N in the backbone
		    try {
			splinesBB[cc].calc(CNC[cc],HBC[cc]);
		    }
		    catch (Exception e) {
			throw new RuntimeException("BBON: " + CNC[cc] + " " + HBC[cc] + "\n" +
						   molecularSystem.get(cc).atom+ "\n");
		    }
		    energyBBON = W_BBPolarSolvate*(splinesBB[cc].s);
		    dSplineDCNC[cc] = W_BBPolarSolvate*splinesBB[cc].s_tag_x;
		    dSplineDHBC[cc] = W_BBPolarSolvate*splinesBB[cc].s_tag_y;
		}
		else if (superType[cc]==2) {	// polar groups in side-chains
		    try {
			splinesSCPolar[cc].calc(CNC[cc],HBC[cc]);
		    }
		    catch (Exception e) {
			throw new RuntimeException("SCP: " + CNC[cc] + " " + HBC[cc] + "\n" +
						   molecularSystem.get(cc).atom+ "\n");
		    }
		    energySCPolar = W_SCPolarSolvate*splinesSCPolar[cc].s;
		    dSplineDCNC[cc] = W_SCPolarSolvate*splinesSCPolar[cc].s_tag_x;
		    dSplineDHBC[cc] = W_SCPolarSolvate*splinesSCPolar[cc].s_tag_y;
		}
		
		
		energy += (energySCPolar + energySCCarbon + energyBBON - 0.5*W_HB*HBCforHBenergy[cc]); // We want to count each HB once.
		if (updateAtoms) 
		    molecularSystem.get(cc).atom.addEnergy(energySCPolar + energySCCarbon + energyBBON - W_HB*HBCforHBenergy[cc]);				
	    }
	}
	return energy;
    }
     
    public void secondPassOverTheNonBondedList(double W_HB){
	int ind1,ind2;
	DistanceList dislist = dm.nonBondedList();
	SolvateDistanceAttribute sigmaValues=null; 
	SolvateDistanceAttributeBetweenPolars sigmaValuesBP = null;
	SolvateDistanceAttributeWithNonPolar  sigmaValuesWNP = null;
	// ************************************
	// Second pass over the non-bonded list
	// ************************************
	for (Distance dis:dislist){
	    sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
	    if (!sigmaValues.isDisInvolvesHydrogens) {
		ind1 = dis.atom1.number;
		ind2 = dis.atom2.number;
		if (sigmaValues.isDisBetween2Polars) {
		    sigmaValuesBP = (SolvateDistanceAttributeBetweenPolars) sigmaValues;
		    if (sigmaValuesBP.sigmHB>0.0) { // The HB related derivatives 
			solvateHB.findBondByPolars(dis.atom1(),dis.atom2()).applyForcesToAtoms(
											       dSplineDHBC[ind1]* sigmaValuesBP.saltBridgeFactorA1 + 
											       dSplineDHBC[ind2]* sigmaValuesBP.saltBridgeFactorA2 +
											       - 0.5 * W_HB * (sigmaValuesBP.saltBridgeFactorForHBenergyA1 + sigmaValuesBP.saltBridgeFactorForHBenergyA2));
		    }
		}
		else {  // Carbon related derivatives 
		    sigmaValuesWNP = (SolvateDistanceAttributeWithNonPolar) sigmaValues;
		    // Doing the self derivatives
		    forceX[ind1] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dx1;
		    forceY[ind1] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dy1;
		    forceZ[ind1] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dz1;
	   			
		    forceX[ind2] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dx2;
		    forceY[ind2] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dy2;
		    forceZ[ind2] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dz2;
		
		    // Doing the cross derivatives
		    forceX[ind2] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dx2;
		    forceY[ind2] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dy2;
		    forceZ[ind2] += dSplineDCNC[ind1] * sigmaValuesWNP.dsigmCa1dz2;
	   			
		    forceX[ind1] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dx1;
		    forceY[ind1] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dy1;
		    forceZ[ind1] += dSplineDCNC[ind2] * sigmaValuesWNP.dsigmCa2dz1;   
		}
	    }  // No hydrogens are involved	
	}   // Of second iteration                 	
    }
    public void assignForcesToAtoms(){
	int cc;
	Atom atom;
	// Finally, the appropriate forces are assigned for every atom. 
	for (cc=0 ; cc<molecularSystem.size() ; cc++) {
	    if (! molecularSystem.get(cc).atom.nowhere()) {	
		atom = molecularSystem.get(cc).atom;
		if (! atom.frozen()) {
		    atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
		    atom.addToFy(-forceY[cc]); // Negating so that it is realy force
		    atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
		}
	    }
	}
    }
}
	

















/*
  This is what I put in the spline parameters:

  cp scModeling/newForClust/SolvateCarbonSplinesCorrection1.txt ../meshi/parameters/meshiPotential/SolvateCarbonSideChainSplines.dat

  cp scModeling/newForClust/SolvatePolarBackboneSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines.dat

  cp scModeling/newForClust/SolvatePolarBackboneSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines_alt2.dat

  cp scModeling/newForClust/SolvatePolarSideChainSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines.dat

  cp scModeling/newForClust/SolvatePolarSideChainSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines_alt2.dat

*/
