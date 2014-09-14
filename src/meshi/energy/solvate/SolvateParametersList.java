package meshi.energy.solvate;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.parameters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.mathTools.Spline1D;
import meshi.util.mathTools.Spline2D;

/**
 * The parameter class required for the SolvateEnergy class. This class reads the parameters from 14 different files.
 * The 14 file names are given as a String array to the constructor. 
 * 
 * The first 10 files, are files that read parameters for the different sigmoid functions used in the calculations 
 * of the carbon and HB indices in SolvateEnergy. See the documentation of the Sigma class to see what are the 
 * parameters required to define a sigmoid. The formats of all these files are similar. Each file contains a
 * 14x14 matrix of doubles (14 is the number of types in Tsai 99'. They are: 1-PRO-N ; 2-backbone-N ; 3-ASN-ND ; 4-LYS-NZ ;
 * 5-backbone-O ; 6-ASP-OD ; 7-PHE-CG ; 8-PHE-CD ; 9-VAL-CB ; 10-backbone-CA ; 11-ALA-CB ; 12-MET-SD ; 13-CYS-SG ; 14-hydrogens).
 * Value i,j in the matrix, is the relevent sigmoid property of atom type j on atom type i. For example, the sixth value in the 
 * second row, corresponds to a sigmoid of hydroxyl oxygen (type 6) on a backbone nitrogen (type 2).     
 *
 * The 13 file types:
 * ------------------
 * 1) SolvateExtResCend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 2) SolvateExtResCp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 3) SolvateExtResCp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 4) SolvateExtResCvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 5) SolvateExtResCvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 6) SolvateExtResHBend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the HB index. 
 *
 * 7) SolvateExtResHBp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the HB index. 
 *
 * 8) SolvateExtResHBp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the HB index. 
 *
 * 9) SolvateExtResHBvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the HB index. 
 *
 * 10) SolvateExtResHBvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the HB index. 
 *
 * 11) SolvateMESHI2Tsai.dat - In MESHI there are 190 defined atom types. In the solvate energy, we sometime use 
 * a reduce representation of 14 atom types as defined in Tsai et. al. 99'. The format of this file is:
 * {MESHI atom type (string)} {Tsai type number (int)}
 * {MESHI atom type (string)} {Tsai type number (int)}  
 * {MESHI atom type (string)} {Tsai type number (int)}
 *  ...
 * {MESHI atom type (string)} {Tsai type number (int)}  
 *
 * 12) SolvateSCpolarSplines.dat - 2D spline parameters for the solvation of side-chain polar atoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are 
 * not side-chain polar atoms. The format of this file is:
 * TRN {2D spline initialization string for TRN (see Spline2D.java)}
 * TRC {zero spline}
 * TRO {2D spline initialization string for TRO (see Spline2D.java)}
 * AH {zero spline)
 * AN {zero spline)
 *  ...
 * CSG {2D spline initialization string for CSG (see Spline2D.java)} 
 *  ...
 * YOH {2D spline initialization string for YOH (see Spline2D.java)}
 *
 * 13) SolvateSCcarbonSplines.dat - 1D spline parameters for the solvation of side-chain NON-polar atoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are 
 * not side-chain NON-polar atoms. The format of this file is:
 * TRN {zero spline}
 * TRC {zero spline}
 * TRO {zero spline}
 * AH {zero spline}
 *  ...
 * DCG {1D spline initialization string for DCG (see Spline1D.java)} 
 *  ...
 * YCZ {1D spline initialization string for YCZ (see Spline1D.java)}
 * YOH {zero spline}
 *
 * 14) SolvateBBpolarSplines.dat - 2D spline parameters for the solvation of backbone polar atoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are 
 * not backbone polar atoms. The format of this file is:
 * TRN {zero spline}
 * TRC {zero spline}
 * TRO {zero spline}
 * AH {zero spline)
 * AN {2D spline initialization string for AN (see Spline2D.java)}
 *  ...
 * AO {2D spline initialization string for AO (see Spline2D.java)} 
 *  ...
 * YOH {zero spline}
 **/

public class SolvateParametersList implements Parameters {
	
	public final int NTsai = 14; // The number of atom types used in Tsai 99'. Any Hydrogen is type 14.
    public final int[] atomicTypeConverter; // From MESHI atom types to Tsai 99'.
    public final double maxEnd; // The maximal distance where any sigmoid is not zero. This value must be less than the Rmax in the distance matrix. 
    public final Spline2D[] scPolarSplines;
    public final Spline2D[] bbSplines;
    public final Spline1D[] scCarbonSplines;
    public final double[][] Cend;     
    public final double[][] Cp1;
    public final double[][] Cp2;
    public final double[][] CvalAtp1;
    public final double[][] CvalAtp2;
    public final double[][] HBend;     
    public final double[][] HBp1;
    public final double[][] HBp2;
    public final double[][] HBvalAtp1;
    public final double[][] HBvalAtp2;

/**
 * The parameter to the constructor is an array of 14 Strings, giving the 14 files required as parameters.
 * See the MeshiPotential class for a detailed list.
 **/    
    public SolvateParametersList(String[] parameterFiles) {
	super();
	int maxAtomType=-1;
	int tmp;
	BufferedReader br;
    StringTokenizer stok;
    String line = "";
		
	Cend = new double[NTsai][NTsai];
	readSigmoidValueFile(Cend,parameterFiles[0]);

	Cp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(Cp1,parameterFiles[1]);

	Cp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(Cp2,parameterFiles[2]);

	CvalAtp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(CvalAtp1,parameterFiles[3]);

	CvalAtp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(CvalAtp2,parameterFiles[4]);

	HBend = new double[NTsai][NTsai];
	readSigmoidValueFile(HBend,parameterFiles[5]);

	HBp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBp1,parameterFiles[6]);

	HBp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBp2,parameterFiles[7]);

	HBvalAtp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBvalAtp1,parameterFiles[8]);

	HBvalAtp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBvalAtp2,parameterFiles[9]);
		
	// Finding the largest ending distance
	double maxEndTmp = -1.0;
	for (int c1=0 ; c1<NTsai ; c1++) 
	for (int c2=0 ; c2<NTsai ; c2++) {
	   if (Cend[c1][c2]>maxEndTmp)
	      maxEndTmp = Cend[c1][c2];
	   if (HBend[c1][c2]>maxEndTmp)
	      maxEndTmp = HBend[c1][c2];
	}
	maxEnd = maxEndTmp;

	// Converting the 190 MESHI atom types to the 14 mentioned in Tsai 99'
	System.out.println("Reading solvation parameter file: " + parameterFiles[10]);
	try{
		// first pass on the file - to find the maximal atom type
		br = new BufferedReader(new FileReader(parameterFiles[10]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = AtomType.type(stok.nextToken().trim()).ordinal();
	    	if (tmp>maxAtomType)
	    	   maxAtomType = tmp;
	    	line = br.readLine();
	    }
	    br.close();
	    atomicTypeConverter = new int[maxAtomType+1];
	    for(int c=0 ; c<atomicTypeConverter.length; c++)
	    	atomicTypeConverter[c] = -1; 
		// second pass on the file - reading the new types
		br = new BufferedReader(new FileReader(parameterFiles[10]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = AtomType.type(stok.nextToken().trim()).ordinal();
	    	atomicTypeConverter[tmp] = Integer.valueOf(stok.nextToken().trim()).intValue()-1;
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}
	
	// Reading the spline parameters for the solvation of side-chain polar groups (There should 
	// be one spline for each residue type).
	scPolarSplines = new Spline2D[maxAtomType+1];
	System.out.println("Reading solvation parameter file: " + parameterFiles[11]);
	try{
		br = new BufferedReader(new FileReader(parameterFiles[11]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = AtomType.type(stok.nextToken().trim()).ordinal();
	    	scPolarSplines[tmp] = new Spline2D(stok);
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}

	// Reading the spline parameters for the solvation of the side-chain non-polar groups (There should 
	// be one spline for each residue type).
	scCarbonSplines = new Spline1D[maxAtomType+1];
	System.out.println("Reading solvation parameter file: " + parameterFiles[12]);
	try{
		br = new BufferedReader(new FileReader(parameterFiles[12]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = AtomType.type(stok.nextToken().trim()).ordinal();
	    	scCarbonSplines[tmp] = new Spline1D(stok);
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}

	// Reading the spline parameters for the solvation of the backbone polar O and N. There should 
	// be two splines (for N and O) for each residue type. The exception is proline that has only O.
	bbSplines = new Spline2D[maxAtomType+1];
	System.out.println("Reading solvation parameter file: " + parameterFiles[13]);
	try{
		br = new BufferedReader(new FileReader(parameterFiles[13]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = AtomType.type(stok.nextToken().trim()).ordinal();
	    	bbSplines[tmp] = new Spline2D(stok);
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}


    }
    
    public Parameters createParameters(String s) {
    	throw new RuntimeException("This method should not be called");
    }
    
    public Parameters parameters(Object obj) {
    	    	throw new RuntimeException("This method should not be called");
    }
    
    
    // Reading a 14x14 file into the relevent parameter.
	private void readSigmoidValueFile(double[][] ar, String filename) {
		BufferedReader br;
		StringTokenizer stok;
		String line = "";
		System.out.println("Reading solvation parameter file: " + filename);
		try{
			br = new BufferedReader(new FileReader(filename));
			for (int ind1=0 ; ind1<NTsai ; ind1++) {
				line = br.readLine();
				if (line==null) 
					throw new RuntimeException("In " + filename + " there should be a " + NTsai + "x" +
					NTsai + " array of doubles");
				stok = new StringTokenizer(line);
				for (int ind2=0 ; ind2<NTsai ; ind2++) 
					ar[ind1][ind2] = Double.valueOf(stok.nextToken().trim()).doubleValue();
			}
			br.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}    
}
