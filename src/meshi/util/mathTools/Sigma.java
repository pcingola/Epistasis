package meshi.util.mathTools;

/**
 * This class host the static method "sigma" which gives a derivable sigmoid function
 * of the following form:
 *
 * x<0   - not defined
 * 0<x<p1 - a linear descent from 1.0 to the value of valAtp1
 * p1<x<p2 - a cubic spline descent from (p1,valAtp1) to (p2,valAtp2)
 * p2<x<end - a quadratic descent from (p2,valAtp2) to (end,0.0).  
 * end<x - 0.0
 * 
 * After the method sigma(...) is run, the public fields s and s_tag are updated:
 * s - the sigmoid value at x.
 * s_tag - the sigmoid first derivative at x.
 *
 * Note: The sigmoid function is derivable for all x=<0.
 **/   

public class Sigma {	
	public static double s = 0.0; // The value of sigma
	public static double s_tag = 0.0; // The derivative of sigma  
	
	
	public static final void sigma(double x, double end, double p1, double p2,
								double valAtp1, double valAtp2) {

		double a1,a2,v1,v2,a,b,c,d; //Auxilary variables 
							
    	if (x>=end) {
    		s = 0;
    		s_tag = 0;
    		return;
    	}
    	if (x>=p2) {
    		s = (x-end)*(x-end)*valAtp2/((p2-end)*(p2-end));
    		s_tag = 2*(x-end)*valAtp2/((p2-end)*(p2-end));
    		return;
    	}
    	if (x>=p1) {
    		v1 = valAtp1;
    		v2 = valAtp2;
    		a1 = (valAtp1-1)/p1;
            a2 = 2*valAtp2/(p2-end); 
            a = (p1*a1-2*v1-p2*a1-p2*a2+p1*a2+2*v2)/(-p2*p2*p2-3*p2*p1*p1+p1*p1*p1+3*p1*p2*p2);
            b = -(2*p1*p1*a2+p1*p1*a1+3*p1*v2-p1*p2*a2-3*p1*v1+p1*p2*a1-2*p2*p2*a1+3*p2*v2-p2*p2*a2-3*p2*v1)/(-p2*p2*p2-3*p2*p1*p1+p1*p1*p1+3*p1*p2*p2);
            c = (2*p1*p1*p2*a1+p1*p1*p2*a2+p1*p1*p1*a2-p2*p2*p1*a1+6*p1*p2*v2-2*p2*p2*p1*a2-6*p1*p2*v1-p2*p2*p2*a1)/(-p2*p2*p2-3*p2*p1*p1+p1*p1*p1+3*p1*p2*p2);
            d = -(p1*p1*p1*p2*a2-p1*p1*p1*v2+p1*p1*p2*p2*a1+3*p1*p1*p2*v2-p1*p1*p2*p2*a2-p1*a1*p2*p2*p2+v1*p2*p2*p2-3*v1*p1*p2*p2)/(-p2*p2*p2-3*p2*p1*p1+p1*p1*p1+3*p1*p2*p2);
       		s = a*x*x*x + b*x*x + c*x + d; 
    		s_tag = 3*a*x*x + 2*b*x + c;
    		return;
    	}
    	if (x>=0.0) {
    		s = 1-(1-valAtp1)/p1*x;
    		s_tag = -(1-valAtp1)/p1;
    		return;
    	}
    	throw new RuntimeException("A negative x-value is not possible. Yet x = "+x);  
	}

}