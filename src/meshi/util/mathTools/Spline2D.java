package meshi.util.mathTools;
import java.util.*;

/**  
 * The class creates a 2D cubic spline object from given parameters, and allows for the 
 * calculation of the spline value and derivative at any given point inside the spline 
 * intervals.
 *
 * After the method calc(x,y) is run, the public fields s s_tag_x and s_tag_y are updated:
 * s - the spline value at x,y.
 * s_tag_x - the spline first derivative at x,y with respect to x.
 * s_tag_y - the spline first derivative at x,y with respect to y.
 *
 * Note:
 * 1) calc will not work for values outside the range [break 1 , break n].
 * 2) The class operates on any set of break points. However, if the break points are 
 * evenly spaced (in both axes) then the calculation of calc(x,y) will be faster.
 * 3) See the constructor documention for more details on how to set up the class properly.  
 *
 **/

public class Spline2D {
	private double[] breaksX;
	private double[] breaksY;
	private double[][][] coefs;
	private boolean evenBreaks = false;
	private double breakIntervalX;
	private double breakIntervalY;
	
	public double s = 0.0;
	public double s_tag_x = 0.0;
	public double s_tag_y = 0.0;

	/**
	 * The constructor can create a specific spline object from a string with the format (Note that following 
	 * 6 lines should appear in the SAME line, but I broke it for clarity). Just to make it perfectly clear,
	 * this object is created from ONE very long test-line:
	 * {Number of breaks in X axis} {Number of breaks in Y axis} {break 1X} ... {break nX} {break 1Y} ... {break nY} 
	 * {square(1,1) constant coef} {square(1,1) X coef} {square(1,1) X^2 coef} {square(1,1) X^3 coefficient}  
	 * {square(1,1) Y coef} {square(1,1) Y*X coef} {square(1,1) Y*X^2 coef} {square(1,1) Y*X^3 coefficient} 
	 * {square(1,1) Y^2 coef} {square(1,1) Y^2*X coef} {square(1,1) Y^2*X^2 coef} {square(1,1) Y^2*X^3 coefficient} 
	 * {square(1,1) Y^3 coef} {square(1,1) Y^3*X coef} {square(1,1) Y^3*X^2 coef} {square(1,1) Y^3*X^3 coefficient} 
	 * {square(1,2) constant coef} {square(1,2) X coef} ... 
	 *
	 * Note:
	 * 1) There is an alternative constructor that requires a tokenizer of a string with 
	 * the above format.
	 * 2) The breaks must increase monotonicly.
	 * 3) The spline coefficients are not verified. Therefore, derivability is obtained only if the 
	 * coefficients are of a derivable spline. 
	 * 4) In order to have numerical stability, the relative precision of the values in the 
	 * dataLine must be at least 1e-10.   
	 **/	
	public Spline2D(String dataLine) {
		this(new StringTokenizer(dataLine));
	}
	
	public Spline2D(StringTokenizer st) {
		int i,j,k;
		int numXbreak = Integer.valueOf(st.nextToken().trim()).intValue();
		int numYbreak = Integer.valueOf(st.nextToken().trim()).intValue();
		breaksX = new double[numXbreak];
		breaksY = new double[numYbreak];
		coefs = new double[numXbreak][numYbreak][16];
		for (i=0 ; i<numXbreak ; i++) 
			breaksX[i] = Double.valueOf(st.nextToken().trim()).doubleValue();
		for (i=0 ; i<numYbreak ; i++) 
			breaksY[i] = Double.valueOf(st.nextToken().trim()).doubleValue();
		for (i=0 ; i<(numXbreak-1) ; i++) 
			for (j=0 ; j<(numYbreak-1) ; j++) 
				for (k=0 ; k<16 ; k++)
					coefs[i][j][k] = Double.valueOf(st.nextToken().trim()).doubleValue();
		// testing for even breaks
		evenBreaks = testForEvenBreaks(breaksX) && testForEvenBreaks(breaksY);
		if (evenBreaks) {
		   breakIntervalX = breaksX[1]-breaksX[0];
		   breakIntervalY = breaksY[1]-breaksY[0];
		}
	}
	
	public final void calc(double x, double y) {
		if ((x<breaksX[0]) || (x>breaksX[breaksX.length-1]) || (y<breaksY[0]) || (y>breaksY[breaksY.length-1])) {
			if ((x>(breaksX[0]-0.001*(breaksX[1]-breaksX[0]))) && (y>(breaksY[0]-0.001*(breaksY[1]-breaksY[0])))) 
				System.out.println("Warning: some of the spline inputs are slightly outside the spline range. Running can continue.");
			else
     			throw new RuntimeException("X or Y are outside the spline range.");		
		}
		if (evenBreaks) {
		   calcEvenBreaks(x,y);
		   return;
		}
		int indexX;
		int indexY;
		for (indexX=0 ; (x>breaksX[indexX+1]) ; indexX++) {}
		for (indexY=0 ; (y>breaksY[indexY+1]) ; indexY++) {}
		double offsetX = x - breaksX[indexX];
		double offsetX2 = offsetX*offsetX;
		double offsetX3 = offsetX2*offsetX;
		double offsetY = y - breaksY[indexY];
		double offsetY2 = offsetY*offsetY;
		double offsetY3 = offsetY2*offsetY;
		s = coefs[indexX][indexY][0] + 
			coefs[indexX][indexY][1]*offsetX + 
			coefs[indexX][indexY][2]*offsetX2 +
			coefs[indexX][indexY][3]*offsetX3 + 
		    coefs[indexX][indexY][4]*offsetY + 
			coefs[indexX][indexY][5]*offsetY*offsetX + 
			coefs[indexX][indexY][6]*offsetY*offsetX2 +
			coefs[indexX][indexY][7]*offsetY*offsetX3 + 
		    coefs[indexX][indexY][8]*offsetY2 + 
			coefs[indexX][indexY][9]*offsetY2*offsetX + 
			coefs[indexX][indexY][10]*offsetY2*offsetX2 +
			coefs[indexX][indexY][11]*offsetY2*offsetX3 + 
		    coefs[indexX][indexY][12]*offsetY3 + 
			coefs[indexX][indexY][13]*offsetY3*offsetX + 
			coefs[indexX][indexY][14]*offsetY3*offsetX2 +
			coefs[indexX][indexY][15]*offsetY3*offsetX3;
			 
		s_tag_x = coefs[indexX][indexY][1] + 
			coefs[indexX][indexY][2]*2*offsetX +
			coefs[indexX][indexY][3]*3*offsetX2 + 
			coefs[indexX][indexY][5]*offsetY + 
			coefs[indexX][indexY][6]*offsetY*2*offsetX +
			coefs[indexX][indexY][7]*offsetY*3*offsetX2 + 
			coefs[indexX][indexY][9]*offsetY2 + 
			coefs[indexX][indexY][10]*offsetY2*2*offsetX +
			coefs[indexX][indexY][11]*offsetY2*3*offsetX2 + 
			coefs[indexX][indexY][13]*offsetY3 + 
			coefs[indexX][indexY][14]*offsetY3*2*offsetX +
			coefs[indexX][indexY][15]*offsetY3*3*offsetX2;

		s_tag_y = coefs[indexX][indexY][4] + 
			coefs[indexX][indexY][5]*offsetX + 
			coefs[indexX][indexY][6]*offsetX2 +
			coefs[indexX][indexY][7]*offsetX3 + 
		    coefs[indexX][indexY][8]*2*offsetY + 
			coefs[indexX][indexY][9]*2*offsetY*offsetX + 
			coefs[indexX][indexY][10]*2*offsetY*offsetX2 +
			coefs[indexX][indexY][11]*2*offsetY*offsetX3 + 
		    coefs[indexX][indexY][12]*3*offsetY2 + 
			coefs[indexX][indexY][13]*3*offsetY2*offsetX + 
			coefs[indexX][indexY][14]*3*offsetY2*offsetX2 +
			coefs[indexX][indexY][15]*3*offsetY2*offsetX3;
	}
	
	private final void calcEvenBreaks(double x, double y) {
		int indexX = (int) ((x-breaksX[0])/breakIntervalX);
		int indexY = (int) ((y-breaksY[0])/breakIntervalY);
		double offsetX =  x - breaksX[0] - indexX*breakIntervalX;
		double offsetX2 = offsetX*offsetX;
		double offsetX3 = offsetX2*offsetX;
		double offsetY =  y - breaksY[0] - indexY*breakIntervalY;
		double offsetY2 = offsetY*offsetY;
		double offsetY3 = offsetY2*offsetY;
		s = coefs[indexX][indexY][0] + 
			coefs[indexX][indexY][1]*offsetX + 
			coefs[indexX][indexY][2]*offsetX2 +
			coefs[indexX][indexY][3]*offsetX3 + 
		    coefs[indexX][indexY][4]*offsetY + 
			coefs[indexX][indexY][5]*offsetY*offsetX + 
			coefs[indexX][indexY][6]*offsetY*offsetX2 +
			coefs[indexX][indexY][7]*offsetY*offsetX3 + 
		    coefs[indexX][indexY][8]*offsetY2 + 
			coefs[indexX][indexY][9]*offsetY2*offsetX + 
			coefs[indexX][indexY][10]*offsetY2*offsetX2 +
			coefs[indexX][indexY][11]*offsetY2*offsetX3 + 
		    coefs[indexX][indexY][12]*offsetY3 + 
			coefs[indexX][indexY][13]*offsetY3*offsetX + 
			coefs[indexX][indexY][14]*offsetY3*offsetX2 +
			coefs[indexX][indexY][15]*offsetY3*offsetX3;
			 
		s_tag_x = coefs[indexX][indexY][1] + 
			coefs[indexX][indexY][2]*2*offsetX +
			coefs[indexX][indexY][3]*3*offsetX2 + 
			coefs[indexX][indexY][5]*offsetY + 
			coefs[indexX][indexY][6]*offsetY*2*offsetX +
			coefs[indexX][indexY][7]*offsetY*3*offsetX2 + 
			coefs[indexX][indexY][9]*offsetY2 + 
			coefs[indexX][indexY][10]*offsetY2*2*offsetX +
			coefs[indexX][indexY][11]*offsetY2*3*offsetX2 + 
			coefs[indexX][indexY][13]*offsetY3 + 
			coefs[indexX][indexY][14]*offsetY3*2*offsetX +
			coefs[indexX][indexY][15]*offsetY3*3*offsetX2;

		s_tag_y = coefs[indexX][indexY][4] + 
			coefs[indexX][indexY][5]*offsetX + 
			coefs[indexX][indexY][6]*offsetX2 +
			coefs[indexX][indexY][7]*offsetX3 + 
		    coefs[indexX][indexY][8]*2*offsetY + 
			coefs[indexX][indexY][9]*2*offsetY*offsetX + 
			coefs[indexX][indexY][10]*2*offsetY*offsetX2 +
			coefs[indexX][indexY][11]*2*offsetY*offsetX3 + 
		    coefs[indexX][indexY][12]*3*offsetY2 + 
			coefs[indexX][indexY][13]*3*offsetY2*offsetX + 
			coefs[indexX][indexY][14]*3*offsetY2*offsetX2 +
			coefs[indexX][indexY][15]*3*offsetY2*offsetX3;
	}

	private boolean testForEvenBreaks(double[] vec) {
		double interval = vec[1]-vec[0];
		for (int c=1 ; c<vec.length ; c++)
			if (Math.abs((vec[c]-vec[0]) - interval)/(interval+(vec[c]-vec[0])) > 1e-13)
				return false;
		return true;
	}
}	
