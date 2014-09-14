package meshi.util;
import meshi.util.rotamericTools.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.geometry.rotamers.*;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.*;
import meshi.energy.sideChainModelingSolvate.*;


/**
 eThis class perform conccurent sidechain modeling as describe in Kalisman and Keasar (2006). There are two modes
 *of operations. One as a stand alone application that prints the new pdb coordinates. Run as follow:
 *java -Xmx300m SCMOD <commands file name> <pdb file name> <number of iterations>
 *
 *The second mode is to give an instance of Protein to the static method "scmod". This mode is usful if ones wish 
 *to do sidechain modeling as part of a larger application. The coordinates of the Protein instance are changed.
 *
 *Important points:
 *-----------------
 *1) The number of iterations can be 0-Inf. If it is 0 all rotamers are put into their most probable back bone 
 *dependent rotamer. If it is 1, the sidecahin are modeled only by their excluded volume clashes and the BBDEP 
 *rotamer probability. If it is 2 or more, also solvation and HB energy terms are included. 
 *2) In the second mode, residues in the protein that have their CA atoms frozen are not modeled. The stand alone 
 *application models all the residues.
 *3) In the second mode, the protein instance must contain ALL the atom including hydrogens.
 *4) This class uses the package "meshi.energy.sideChainModelingSolvate" which is a truncated form of the package
 *"meshi.energy.solvate", that is suitable for this application ONLY.
 **/


public final class Scmod extends MeshiProgram implements KeyWords, MeshiPotential { /**
										     * The implemented
										     * interfaces defines the 
										     * names of atom and residue 
										     * types. 
										     **/
    private static DunbrackLib lib = null; 
    private static double excludedVolumeWeight;
    private static double safeRadius = 10; 
    //                            W  Y  F L M  P  V  I H E K R  D Q  N  C T  S  
    private static int[] seder = {18,19,4,9,10,12,17,7,6,3,8,14,2,13,11,1,16,15}; /*The order (seder in Hebrew) of the residue
										    types*/

    private static double[] w1;
    private static double[] w2;
    private static double[] w3;
    private static double[] w4;
    private static double[] w5;
    private static double[] w6;

    // Weights for iterations 2 till infinity
    // 20 weights for EV1
    private static double[] i2w1 = {0.00,	100.14,	26.27,	54.54,	355.74,	0.00,	262.64,	159.10,	160.97,	400.19,	42.37,	310.64,	228.00,	42.01,	104.85,	627.60,	92.48,	154.73,	1538.80, 5320.36};
    // 20 weights for EV2
    private static double[] i2w2 = {0.00,	0.00,	8.99,	36.93,	8.10,	0.00,	25.98,	114.81,	58.02,	30.58,	18.22,	10.55,	137.75,	53.84,	17.03,	0.00,	0.00,	0.00,	17.42,	26.12};	
    // 20 weights for EV3
    private static double[] i2w3 = {0.00,	0.00,	0.00,	1.40,	2.11,	0.00,	2.99,	0.00,	-0.0,	0.00,	2.64,	0.00,	0.00,	0.45,	0.05,	0.00,	0.00,	0.00,	0.55,	2.72};	
    // 20 weights for BBDEP rotamer prob
    private static double[] i2w4 = {0.00,	3.17,	1.26,	0.77,	1.25,	0.00,	1.29,	2.44,	0.77,	2.58,	0.89,	1.10,	4.10,	0.80,	0.66,	2.68,	5.40,	6.29,	1.32,	1.24};	
    // 20 weights for Solvate
    private static double[] i2w5 = {0.00,	0.78,	0.15,	-0.0,	0.50,	0.00,	0.09,	3.81,	0.05,	4.10,	1.10,	-0.0,	0.35,	-0.0,	-0.0,	0.82,	1.23,	6.52,	0.60,	0.65};	
    // 20 weights for HB
    private static double[] i2w6 = {0.00,	0.78,	0.66,	0.50,	0.00,	0.00,	0.83,	0.00,	1.14,	0.00,	0.00,	0.63,	0.00,	0.44,	0.32,	1.27,	1.68,	0.00,	1.18,	1.08};	

    // Weights for the first iteration
    // 20 weights for EV1
    private static double[] i1w1 = {0.00,	6.76,	2.86,	4.68,	-0.0,	0.00,	1.04,	4.52,	3.53,	6.20,	3.03,	-0.0,	77.30,	18.71,	15.62,	22.50,	8.02,	6.12,	-0.0,	-0.0};	
    // 20 weights for EV2
    private static double[] i1w2 = {0.00,	0.00,	1.98,	4.89,	7.19,	0.00,	0.79,	3.94,	0.69,	1.72,	1.39,	7.31,	30.13,	34.60,	10.19,	0.00,	0.00,	0.00,	0.20,	1.05};	
    // 20 weights for EV3
    private static double[] i1w3 = {0.00,	0.00,	0.00,	0.97,	0.16,	0.00,	0.07,	0.00,	0.04,	0.00,	0.26,	0.00,	0.00,	0.39,	0.17,	0.00,	0.00,	0.00,	0.04,	0.03};	
    // 20 weights for BBDEP rotamer prob
    private static double[] i1w4 = {0.00,	2.93,	1.19,	0.73,	1.25,	0.00,	1.23,	2.25,	0.75,	2.40,	0.81,	1.03,	4.04,	0.77,	0.65,	2.54,	5.18,	6.02,	1.30,	1.24};	
    // 20 weights for solvate
    private static double[] i1w5 = {0.00,	0.0,	0.0,	0.0,	0.0,	0.00,	0.0,	0.0,	0.0,	0.00,	0.00,	0.0,	0.0,	0.0,	0.0,	0.00,	0.00,	0.00,	0.0,	0.0};	
    // 20 weights for HB
    private static double[] i1w6 = {0.00,	0.0,	0.0,	0.0,	0.0,	0.00,	0.0,	0.0,	0.0,	0.00,	0.00,	0.0,	0.0,	0.0,	0.0,	0.00,	0.00,	0.00,	0.0,	0.0};




    public static void scmod(CommandList commands, Protein protein , int maxIter)  {
	scmod(commands, protein , maxIter, 0.5); // the original function
    }
    public static void scmod(CommandList commands, Protein protein , int maxIter, double excludedVolumeWeight)  {
	scmod(commands,  protein, maxIter,excludedVolumeWeight, null);
    }
    public static void scmod(CommandList commands, Protein protein , int maxIter, double excludedVolumeWeight, 
			     TotalEnergy totalEnergy)  {
	int ind;
	int pred;
	int numOfRotamers;
	double[] rotamer;
	double probrot;
	double score,minscore;
	double[] badPoints;
	DistanceMatrix distanceMatrix;
	TotalEnergy energy;
	SideChainSolvateEnergy energyTerm;

	Scmod.excludedVolumeWeight = excludedVolumeWeight;

	// freezing completely residues that have their CA's frozen.
	for (int residueIndex=0; residueIndex<protein.residues().size() ; residueIndex++) {
	    if (! protein.residues().get(residueIndex).dummy())
		if (protein.residues().get(residueIndex).ca().frozen())
		    protein.residues().get(residueIndex).atoms().freeze();
		
	}
	 	
	ExcludedVolParametersList parametersList =
	    new ExcludedVolParametersList(commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/" +
					      EXCLUDED_VOL_PARAMETERS);

	if (lib == null) lib = new DunbrackLib(commands,1.0,30);
		
	if (totalEnergy == null) 
	    distanceMatrix = new DistanceMatrix(protein.atoms().molecularSystem(), 5.5, 2.0 , DistanceMatrix.DEFAULT_EDGE(5.5, 2.0),4);
	else    distanceMatrix = totalEnergy.distanceMatrix();

	double[][] pp = RotamericTools.putIntoRot1(protein,distanceMatrix,lib);
	distanceMatrix = new DistanceMatrix(protein.atoms().molecularSystem(), 5.5, 2.0 , DistanceMatrix.DEFAULT_EDGE(5.5, 2.0), 4); 

   
	if (totalEnergy == null) {
	    // The energy creators that will be used from now on!
	     EnergyCreator[] energyCreators1 = { new SideChainSolvateCreator() };
	     energy = new TotalEnergy(protein, distanceMatrix, energyCreators1, commands);
	     energyTerm = (SideChainSolvateEnergy) energy.getEnergyTerm(new SideChainSolvateEnergy());
	}
	else {
	    energy = totalEnergy;
	    energy.off();
	    energyTerm = (SideChainSolvateEnergy) energy.getEnergyTerm(new SideChainSolvateEnergy());    
	    energyTerm.on();
	}

	for (int stage=0 ; stage<maxIter ; stage++) {
	    System.out.println("\nDoing iteration: " + (stage+1));
	    if (stage==0) {
		w1 = i1w1;
		w2 = i1w2;    
		w3 = i1w3;    
		w4 = i1w4;    
		w5 = i1w5;    
		w6 = i1w6;
	    }
	    else {
		w1 = i2w1;
		w2 = i2w2;    
		w3 = i2w3;    
		w4 = i2w4;    
		w5 = i2w5;    
		w6 = i2w6;
	    }    	
	    for (int residueTypeIndex=0 ; residueTypeIndex<seder.length ; residueTypeIndex++) {
		int residueType = seder[residueTypeIndex];
		for (int residueIndex=0 ; residueIndex<pp.length ; residueIndex++) 
		    if ((pp[residueIndex]!=null) && (pp[residueIndex][2]==residueType) && 
			(!protein.residue(residueIndex).ca().frozen()) && (!protein.residue(residueIndex).ca().nowhere())) {
			minscore = 1e15;
			pred = -1;
			if ((stage>0) & (energyTerm != null)){
			    energyTerm.inactivateFarFromAtom(protein.residue(residueIndex).ca() , safeRadius);
			}
			//*****************************************************************************    
			      // Going over the rotamers
			      numOfRotamers = lib.getRotamerNum((int) pp[residueIndex][2] , pp[residueIndex][0] , pp[residueIndex][1]);
			for (ind=0 ; ind<numOfRotamers ; ind++) {
			    rotamer = lib.getRotamer((int) pp[residueIndex][2] , pp[residueIndex][0] , pp[residueIndex][1] , ind);
			    probrot = lib.getRotamerProb((int) pp[residueIndex][2] , pp[residueIndex][0] , pp[residueIndex][1] , ind);

			    ResidueBuilder.build(protein.residue(residueIndex), 
						 protein.residue(residueIndex).type, 
						 rotamer);
    
                           
			    badPoints = calcVariousEnergies(residueIndex,probrot,stage,protein,energy,distanceMatrix,energyTerm,parametersList);     
			    score = w1[(int) pp[residueIndex][2]]*badPoints[0] +
				w2[(int) pp[residueIndex][2]]*badPoints[1] + 
				w3[(int) pp[residueIndex][2]]*badPoints[2] - 
				w4[(int) pp[residueIndex][2]]*Math.log(0.002+badPoints[3]) +
				w5[(int) pp[residueIndex][2]]*badPoints[4] + 
				w6[(int) pp[residueIndex][2]]*badPoints[5]; 
			    if (score<minscore) {
				pred = ind;
				minscore = score;
			    }
			}  // The loop on the rotamers 
			// END - Going over the rotamers
			//*****************************************************************************    
    
			      if (pred<0)
			      throw new RuntimeException("Serious problems " + pred);
           
			rotamer = lib.getRotamer((int) pp[residueIndex][2] , pp[residueIndex][0] , pp[residueIndex][1] ,  pred );
			ResidueBuilder.build(protein.residue(residueIndex),protein.residue(residueIndex).type,rotamer);
	
    
		    } // Going for the next residue.
	    } // Going for the next residue type.
	} // Going for the next iteration.
        
    }  // scmod


    protected static double[] calcVariousEnergies(int resnum , double prob, int stage , Protein protein , 
						  TotalEnergy energy , DistanceMatrix distanceMatrix, SideChainSolvateEnergy energyTerm , ExcludedVolParametersList parametersList) {
	double tmpenergy;
	ExcludedVolParameters parameters;
	double[] result = new double[6];
	for (int c=0 ; c<6 ; c++)
	    result[c] = 0;

	result[3] = prob;
	if (stage==0) {
	    try { energy.update(); }
	    catch (Exception e) { throw new RuntimeException(e); }
	}
	else {
	    result[4] = energy.evaluate();
	    if (energyTerm != null) result[5] = energyTerm.getNonWeightedHBenergy();
	    else result[5] = 0;
	}

	// reseting the relevent atom energies
	AtomList al = protein.residue(resnum).atoms();
	Atom atom1,atom2;
	for(Atom atom:al)
	    atom.resetEnergy();
	// Calculating the EV on the relevent atoms  
	for(Distance distance:distanceMatrix.nonBondedList()){
	    if (!distance.mode().frozen) {
	    	atom1 = distance.atom1();
	    	atom2 = distance.atom2();
		if ((atom1.residueNumber() == resnum) ||
		    (atom2.residueNumber() == resnum)) {
		    parameters = (ExcludedVolParameters) parametersList.parameters(distance);
		    if (distance.distance()<parameters.sigma) {
			if ((atom1.name().length()==1) || (atom2.name().length()==1) || atom1.name().equals("CA") ||
			    atom2.name().equals("CA") || atom1.name().equals("CB") || atom2.name().equals("CB")) {
			    tmpenergy = excludedVolumeWeight*3*(parameters.sigma-distance.distance());
                	}
                	else {
			    tmpenergy = excludedVolumeWeight*(parameters.sigma-distance.distance());
               		}
               		atom1.addEnergy(tmpenergy);
               		atom2.addEnergy(tmpenergy);
		    }
		}
	    }
	}
	// Getting the energy values from the sidechain
	for (Atom atom:al)
	    if ((atom.name().length()>1) && 
		(atom.name().charAt(1)=='G')) {
		result[0] += atom.energy()*atom.energy()*atom.energy();
	    }
	    else if ((atom.name().length()>1) && 
		     (atom.name().charAt(1)=='D')) { 
          	result[1] += atom.energy()*atom.energy()*atom.energy();
	    }
	    else if ((atom.name().length()>1) && 
		     (!atom.name().equals("CA")) && 
		     (!atom.name().equals("CB")) &&
		     (atom.name().charAt(1)!='X')) {
          	result[2] += atom.energy()*atom.energy()*atom.energy();
	    }
 

	return result;
    }
    
} // of class

