package meshi.symmetryComplex.energy.edmEnergy;

import meshi.util.file.MeshiLineReader;

import java.io.*;

/* A splined polynomial for a protein's electron density map. */

public class EDMSplinedPolynomial {
	
	private static final int NUM_DIMENSIONS = 3;
		/* all EDM polynomials have three dimensions */
	
	private double breaks[][];		/* polynomial breaks for each dimension */
	private double coefs[];			/* polynomial coefficients */

	/** Create a new EDMSplinedPolynomial from parameters.
	 * 
	 * @param breaks each break has different coefficients
	 * @param coefs coefficienst by breaks
	 */
	public EDMSplinedPolynomial(
			double breaks[][], double coefs[] ) {	
		this.breaks = breaks;
		this.coefs = coefs;
	}
	
	/** Create a new EDMSplinedPolynomial from a file. */
	public EDMSplinedPolynomial( String edmDataFileNum ) throws IOException {
		MeshiLineReader mlr = new MeshiLineReader( edmDataFileNum );

		/* read all breaks arrays */
		breaks = new double[NUM_DIMENSIONS][];
		int breakSize;
		int numCoefs = 1;
		for( int i=0; i<NUM_DIMENSIONS; i++ ) {
			/* read breaks size */
			breakSize = new Integer( mlr.readLine("#") ).intValue();
			
			/* set up breaks array read contents of break */
			breaks[i] = new double[breakSize];
			for( int j=0; j<breakSize; j++ )
				breaks[i][j] = new Double( mlr.readLine("#") ).doubleValue();
			
			/* update number of coefficients */
			numCoefs = numCoefs*(breakSize-1)*4;
		}
		
		/* set up coefficients array and read all coefficients */
		coefs = new double[numCoefs];
		for( int i=0; i<numCoefs; i++ )
			coefs[i] = new Double( mlr.readLine("#") ).doubleValue();
	}
	
	/** Evaluate value of polynomial or any of its first derivations at given point.
	 * 
	 * @param derivVar variable to be derived (zero for none)
	 * @param args list of variable assignments
	 * @return evaluation of polynomial
	 */
	public double value( int derivVar, double ... args ) {
		double result = 0;		/* result of calculation */
		int i,j,k;
		
		/* verify number of arguments is correct */
		if( args.length != NUM_DIMENSIONS )
			throw new RuntimeException( "polynomial number of variables mismatch" );
		/* verify we're not asking to derive a non existent variable */
		if( derivVar > NUM_DIMENSIONS )
			throw new RuntimeException( "derived variable identifier doesn't exist in this polynomial" );
		
		/* calculate 3D polynomial */
		double val1 = args[0];
		double val2 = args[1];
		double val3 = args[2];
		int bin1 = findBin( breaks[0], val1 );
		int bin2 = findBin( breaks[1], val2 );
		int bin3 = findBin( breaks[2], val3 );
			
		/* verify a bin was found */
		if( bin1 == -1 || bin2 == -1 || bin3 == -1 )
			throw new RuntimeException(
					"unable to locate correct bin for: (" + val1 +
					"," + val2 + "," + val3 + ")" );

		/* convert values to distance from beginning of bin */
		val1 = val1-breaks[0][bin1];
		val2 = val2-breaks[1][bin2];
		val3 = val3-breaks[2][bin3];
			
		/* calculation acceleration */
		int bin3x4 = bin3*4;
		int co3lenx4 = (breaks[2].length-1)*4;
		int bin2xco3lenx16 = bin2*co3lenx4*4;
		int co2lenxco3lenx16 = (breaks[1].length-1)*co3lenx4*4;
		int bin1xco2lenxco3lenx64 = bin1*co2lenxco3lenx16*4;
			
		/* calculation for tri-variate splined polynomial */
		if( derivVar == 0 )
			/* no derivation required */
			for( i=0; i<4; i++ )
				for( j=0; j<4; j++ )
					for( k=0; k<4; k++ )
						result = result +
									coefs[ bin1xco2lenxco3lenx64 + co2lenxco3lenx16*i + bin2xco3lenx16 + co3lenx4*j + bin3x4 + k] * 
									quickPower(val1, 4-(i+1)) *
									quickPower(val2, 4-(j+1)) *
									quickPower(val3, 4-(k+1));
			else if( derivVar == 1 )
				/* derive according to first variable */
				for( i=0; i<4; i++ )
					for( j=0; j<4; j++ )
						for( k=0; k<4; k++ )
							result = result +
										coefs[ bin1xco2lenxco3lenx64 + co2lenxco3lenx16*i + bin2xco3lenx16 + co3lenx4*j + bin3x4 + k] * 
										quickPower(val1, 3-(i+1)) * (4-(i+1)) *
										quickPower(val2, 4-(j+1)) *
										quickPower(val3, 4-(k+1));
			else if( derivVar == 2 )
				/* derive according to second variable */
				for( i=0; i<4; i++ )
					for( j=0; j<4; j++ )
						for( k=0; k<4; k++ )
							result = result +
										coefs[ bin1xco2lenxco3lenx64 + co2lenxco3lenx16*i + bin2xco3lenx16 + co3lenx4*j + bin3x4 + k] * 
										quickPower(val1, 4-(i+1))*
										quickPower(val2, 3-(j+1)) * (4-(j+1))  *
										quickPower(val3, 4-(k+1));
			else if( derivVar == 3 )
				/* derive according to third variable */
				for( i=0; i<4; i++ )
					for( j=0; j<4; j++ )
						for( k=0; k<4; k++ )
							result = result +
										coefs[ bin1xco2lenxco3lenx64 + co2lenxco3lenx16*i + bin2xco3lenx16 + co3lenx4*j + bin3x4 + k] * 
										quickPower(val1, 4-(i+1)) *
										quickPower(val2, 4-(j+1)) *
										quickPower(val3, 3-(k+1)) * (4-(k+1));
		
		return result;
	}

	/** Find break index (bin) for value.
	 * 
	 * @return left bound of break, -1 if no break found.
	 */
	public int findBin( double breaks[], double value ) {
		/* update value to be in bin range */
		if( value < breaks[0] )
			value = value + 2*Math.PI;
		else if( value > breaks[breaks.length-1] )
			value = value - 2*Math.PI;
		
		/* search for value in all the bins */
		for( int binCount=0; binCount<breaks.length; binCount++ )
			if( breaks[binCount] <= value && breaks[binCount+1] > value )
				return binCount;
		return -1;
	}

	/** Raises a double by the power, with 0 <= power <= 3. For the special
	 * exception of power=-1, returns 0 (as it eases polynomial derivation
	 * calculations).
	 */
	private static double quickPower( double torsion, int power ) {
		switch( power ) {
		case -1: return 0;
		case 0: return 1;
		case 1: return torsion;
		case 2: return torsion*torsion;
		case 3: return torsion*torsion*torsion;
		}
		
		return -1;
	}

}
