package meshi.energy.simpleEnergyTerms.alphaTorsion;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import java.io.*;
import meshi.util.filters.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import java.util.*;

/**
 *Parsing a line in the parameter file of the alpha torsion energy.
 *Each line in the parameter file of the alpha torsion energy is with the following format:
 *{AA type}
 *{Weight = the parabola height}
 *{Starting torsion value value of the HELIX secondary structure}
 *{Finish torsion value value of the HELIX secondary structure}
 *{Starting torsion value value of the SHEET secondary structure}
 *{Finish torsion value value of the SHEET secondary structure}
 **/
 
 public class AlphaTorsionParameters implements Parameters {

    public String aaLetter;
    public double weightAA;
    public double startAlphaHELIX;
    public double endAlphaHELIX;
    public double startAlphaSHEET;
    public double endAlphaSHEET;


    public AlphaTorsionParameters(String line) {
    	StringTokenizer stok;
       	stok = new StringTokenizer(line);
        aaLetter = stok.nextToken().trim();
        weightAA = Double.valueOf(stok.nextToken()).doubleValue();
        if (weightAA < 0.0)
            throw new RuntimeException("Weight must be non-negative\n");      	        
        startAlphaHELIX = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaHELIX < -Math.PI) || (startAlphaHELIX > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaHELIX = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaHELIX < -Math.PI) || (endAlphaHELIX > Math.PI) || (endAlphaHELIX < startAlphaHELIX))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        startAlphaSHEET = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaSHEET < -Math.PI) || (startAlphaSHEET > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaSHEET = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaSHEET < -Math.PI) || (endAlphaSHEET > Math.PI) || (endAlphaSHEET > startAlphaSHEET))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
    }


	
    public String toString() {
	return "AlphaTorsionParameters\n"+
	    "\t AA Letter   = "+aaLetter+"\n";
    }

}
