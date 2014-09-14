package meshi.energy.pairwiseNonBondedTerms.CoulombElectrostatics;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.util.file.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.parameters.*;
import java.util.*;

/** 
* 
*
* Creates a CoulombElectrostaticTerm for the protein with a given weight.
* The creator: 
*  1. Reads data from the file, and creates the charge parameters list.
*  2. Initializes the DielectricConstsnt from the command file (given by user)
*  3. Creates a CoulombElectrostatics object that will evaluate the total electrostatics energy.
*	
* Constructors
*	public CoulombElectrostaticCreator(double weight)	-  A CoulombElectrostaticCreator that takes a weight
*                                                          parameter for energy term.                                                       
*	public CoulombElectrostaticCreator()    		    -  The default CoulombElectrostaticCreator.
*
* Object methods 
*	AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands)
*   public double getDielectricConstsnt(CommandList commands)
*
*  
**/

public class CoulombElectrostaticCreator extends EnergyCreator  implements KeyWords {
    private static double[] charges = null;
	/** 
	 * A CoulombElectrostaticCreator constructor that takes the weight
	 * given to this energy term as a parameter. 
	 * @param weight the weight given to this energy term.
	 **/ 
    public CoulombElectrostaticCreator(double weight) {
  		super(weight);
    }
	
    /** 
     * default constructor
     **/
    public CoulombElectrostaticCreator() {
  		super(ELECTROSTATICS);
    }
	
	/** 
	 * the creator: 
	 * 1. Reads data from the file "CHARGE_PARAMETERS" and creates the charge parameters list.
	 * 2. Initialize the DielectricConstsnt from the command file (given by user).
	 * 3. Creates a CoulombElectrostatics object that will evaluate the total electrostatics energy.
	 * @param protein
	 * @param distanceMatrix
	 * @param commands
	 * @return AbstractEnergy
	 **/	    
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					                        CommandList commands) {
	if (charges == null) {
	    charges = new double[AtomType.values().length];
	    for (int i = 0 ; i < AtomType.values().length;i++)
		charges[i] = Double.MAX_VALUE;
	    String parametersFileName = parametersDirectory(commands)+"/"+ELECTROSTATICS_PARAMETERS;
	    MeshiLineReader lines = new MeshiLineReader(parametersFileName);
	    String line;
	    while ((line = lines.readLine("#")) != null) {
		StringTokenizer st = new StringTokenizer(line);
		int typeOrdinal = AtomType.type(st.nextToken()).ordinal();
		double charge = Double.valueOf(st.nextToken().trim()).doubleValue();
		charges[typeOrdinal] = charge;
	    }
	    try {
		lines.close();
	    } catch (Exception ex) {throw new RuntimeException("Cannot close "+lines+"\n"+ex);}
	}
	double dielectricConstant = getDielectricConstant(commands);
	return new CoulombElectrostatics(distanceMatrix, weight(), dielectricConstant, charges);
    }
    
	/**
	 * This takes the value of the dielectricConstant from the commands file (the user may choose to change it)
	 * @param commands
	 * @return double 
	 **/ 
	public double getDielectricConstant(CommandList commands) {
		return commands.firstWordFilter(ELECTROSTATICS).secondWord(DIELECTRIC_CONSTANT).thirdWordDouble();
    }
	
}//end
