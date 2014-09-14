package meshi.util.mathTools;

/**
 * This class host the static method "ramp" which gives a derivable ramp function
 * of the following form:
 *
 * x<0   - not defined
 * 0<x<S - The function: (S-x)^2/(S-x+alpha) 
 * S<x - 0.0
 * 
 * Where: S is the ramp begining, and alpha is length of the derivable transition region at the ramp begining.
 *
 * After the method sigma(...) is run, the public fields r and r_tag are updated:
 * r - the ramp value at x.
 * r_tag - the ramp first derivative at x.
 *
 * Note: The ramp function is derivable for all x=<0.
 **/   

public class Ramp {	
	public static double r = 0.0; // The value of ramp
	public static double r_tag = 0.0; // The derivative of ramp  
	
	
	public static final void ramp(double x, double S, double alpha) {

		double a1,a2,v1,v2,a,b,c,d; //Auxilary variables 
							
    	if (x>=S) {
    		r = 0;
    		r_tag = 0;
    		return;
    	}
    	if (x>=0.0) {
    		a1 = 1/(S-x+alpha);
    		r  = (S-x)*(S-x)*a1;
    		r_tag = -1 + alpha*alpha*a1*a1;
    		return;
    	}
    	throw new RuntimeException("A negative x-value is not possible.");  
	}

}