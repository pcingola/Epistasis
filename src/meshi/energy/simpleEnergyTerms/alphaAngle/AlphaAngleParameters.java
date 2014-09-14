package meshi.energy.simpleEnergyTerms.alphaAngle;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import java.io.*;
import meshi.util.filters.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import java.util.*;

/**
 *Parsing a line in the parameter file of the alpha angle energy.
 *Each line in the parameter file of the alpha angle energy is with the following format:
 *{AA type}
 *{Weight = the parabola height}
 *{Starting angle value of the ALL secondary structure}
 *{Finish angle value of the ALL secondary structure}
 *{Starting angle value of the HELIX secondary structure}
 *{Finish angle value of the HELIX secondary structure}
 *{Starting angle value of the SHEET secondary structure}
 *{Finish angle value of the SHEET secondary structure}
 *{Starting angle value of the COIL secondary structure}
 *{Finish angle value of the COIL secondary structure}
 **/
 
 public class AlphaAngleParameters implements Parameters {

    public String aaLetter;
    public double weightAA;   
    public double startAlphaALL;
    public double endAlphaALL;
    public double startAlphaHELIX;
    public double endAlphaHELIX;
    public double startAlphaSHEET;
    public double endAlphaSHEET;
    public double startAlphaCOIL;
    public double endAlphaCOIL;


    public AlphaAngleParameters(String line) {
    	StringTokenizer stok;
       	stok = new StringTokenizer(line);
        aaLetter = stok.nextToken().trim();
        weightAA = Double.valueOf(stok.nextToken()).doubleValue();
        if (weightAA < 0.0)
            throw new RuntimeException("Weight must be non-negative\n");      	
        startAlphaALL = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaALL < 0.0) || (startAlphaALL > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaALL = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaALL < 0.0) || (endAlphaALL > Math.PI) || (endAlphaALL < startAlphaALL))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        startAlphaHELIX = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaHELIX < 0.0) || (startAlphaHELIX > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaHELIX = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaHELIX < 0.0) || (endAlphaHELIX > Math.PI) || (endAlphaHELIX < startAlphaHELIX))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        startAlphaSHEET = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaSHEET < 0.0) || (startAlphaSHEET > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaSHEET = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaSHEET < 0.0) || (endAlphaSHEET > Math.PI) || (endAlphaSHEET < startAlphaSHEET))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        startAlphaCOIL = Double.valueOf(stok.nextToken()).doubleValue();
        if ((startAlphaCOIL < 0.0) || (startAlphaCOIL > Math.PI))
            throw new RuntimeException("Wrong values in the parameters file\n");      	
        endAlphaCOIL = Double.valueOf(stok.nextToken()).doubleValue();
        if ((endAlphaCOIL < 0.0) || (endAlphaCOIL > Math.PI) || (endAlphaCOIL < startAlphaCOIL))
            throw new RuntimeException("Wrong values in the parameters file\n");      	    	
    }


	
    public String toString() {
	return "AlphaAngleParameters\n"+
	    "\t AA Letter   = "+aaLetter+"\n";
    }

}
