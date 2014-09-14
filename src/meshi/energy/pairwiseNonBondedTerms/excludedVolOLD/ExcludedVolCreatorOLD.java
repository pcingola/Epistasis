package meshi.energy.pairwiseNonBondedTerms.excludedVolOLD;
import meshi.energy.*;
import meshi.util.*;
import meshi.util.file.*;
import meshi.molecularElements.*;
import meshi.geometry.*;
import meshi.parameters.*;
import meshi.util.filters.*;
import java.util.*;

public class ExcludedVolCreatorOLD extends EnergyCreator  implements KeyWords {
    private static double[][] parameters = null;

    private double Rfac = 1.0;
    private Filter filter = null;
 
    public ExcludedVolCreatorOLD(double weight,double Rfac) {
  	super(weight);
  	this.Rfac = Rfac;
    }
    public ExcludedVolCreatorOLD(double Rfac) {
  	super(EXCLUDED_VOL);
  	this.Rfac = Rfac;
    }

    public ExcludedVolCreatorOLD() {
  	super(EXCLUDED_VOL);
  	}
       
    public ExcludedVolCreatorOLD(Filter filter,double Rfac) {
	super(EXCLUDED_VOL);
	this.Rfac = Rfac;
	this.filter = filter;
    }
       

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parameters == null) {
	    parameters = new double[AtomType.values().length][AtomType.values().length];
	    for (int i = 0 ; i < AtomType.values().length;i++)
		for (int j = 0 ; j < AtomType.values().length;j++) {
		    parameters[i][j] = Double.MAX_VALUE;
		}
	    String parametersFileName = parametersDirectory(commands)+"/"+EXCLUDED_VOL_PARAMETERS;
	    MeshiLineReader lines = new MeshiLineReader(parametersFileName);
	    String line;
	    while ((line = lines.readLine("#")) != null) {
		    StringTokenizer st = new StringTokenizer(line);
		    int typeOrdinal1 = AtomType.type(st.nextToken()).ordinal();
		    int typeOrdinal2 = AtomType.type(st.nextToken()).ordinal();
		    double sigma = Double.valueOf(st.nextToken().trim());
		    parameters[typeOrdinal1][typeOrdinal2] = parameters[typeOrdinal2][typeOrdinal1] = sigma;
	    }
	    try {
		    lines.close();
	    }
        catch (Exception ex) {throw new RuntimeException("Cannot close "+lines+"\n"+ex);}
	}
	return new ExcludedVolOLD(distanceMatrix, Rfac, weight(),filter, parameters);
    }
}
