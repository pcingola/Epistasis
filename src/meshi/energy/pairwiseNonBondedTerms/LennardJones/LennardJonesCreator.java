package meshi.energy.pairwiseNonBondedTerms.LennardJones;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.util.string.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.*;
import java.util.*;
import meshi.parameters.*;

public class LennardJonesCreator extends EnergyCreator  implements KeyWords {
    private static double[][][] parameters = null;
    protected String parametersFileName;
    private static int nInstances = 0;
    public LennardJonesCreator(double weight) {
  	super(weight);
	parametersFileName = LENNARD_JONES_PARAMETERS;
	if (nInstances >= 1) throw new RuntimeException("No more than one LennardJonesCreator instance (1)");
	nInstances++;
    }
    public LennardJonesCreator() {
  	super(LENNARD_JONES);
	parametersFileName = LENNARD_JONES_PARAMETERS;
	if (nInstances >= 1) throw new RuntimeException("No more than one LennardJonesCreator instance (2)");
	nInstances++;
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parameters == null) {
	    parameters = new double[AtomType.values().length][AtomType.values().length][2];
	    for (int i = 0 ; i < AtomType.values().length;i++)
		for (int j = 0 ; j < AtomType.values().length;j++) {
		    parameters[i][j][0] = Double.MAX_VALUE;
		    parameters[i][j][1] = Double.MAX_VALUE;
		}
	    String FullParametersFileName = parametersDirectory(commands)+"/"+parametersFileName;
	    MeshiLineReader lines = new MeshiLineReader(FullParametersFileName);
	    String line;
	    while ((line = lines.readLine("#")) != null) {
		StringTokenizer st = new StringTokenizer(line);
		int typeOrdinal1 = AtomType.type(st.nextToken()).ordinal();
		int typeOrdinal2 = AtomType.type(st.nextToken()).ordinal();
		Double epsilon = Double.valueOf(st.nextToken().trim()).doubleValue();
		Double sigma = Double.valueOf(st.nextToken().trim()).doubleValue();
		parameters[typeOrdinal1][typeOrdinal2][0] = parameters[typeOrdinal2][typeOrdinal1][0] = epsilon;
		parameters[typeOrdinal1][typeOrdinal2][1] = parameters[typeOrdinal2][typeOrdinal1][1] = sigma; 
	    }
	    try {
		lines.close();
	    } catch (Exception ex) {throw new RuntimeException("Cannot close "+lines+"\n"+ex);}
	}
	return new LennardJones(distanceMatrix, weight(), parameters);
    }
}
