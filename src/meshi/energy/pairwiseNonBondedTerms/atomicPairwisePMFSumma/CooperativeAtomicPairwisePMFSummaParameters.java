package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.parameters.*;
import meshi.util.string.*;
import meshi.util.*;
import meshi.util.file.*;
import java.util.*;

public class CooperativeAtomicPairwisePMFSummaParameters implements KeyWords, MeshiPotential {
    protected final double[] minAvgSigma = new double[AtomType.values().length];
    
    public CooperativeAtomicPairwisePMFSummaParameters(CommandList commands ) {
	//String          parametersFileName = commands.firstWord( PARAMETERS_DIRECTORY ).secondWord() + "/" + COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_PARAMETERS;
    String          parametersFileName = commands.firstWord( PARAMETERS_DIRECTORY ).secondWord() + "/meshiPotential/" + commands.firstWord( COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_FILENAME).secondWord();
        
    System.out.println("CooperativeAtomicPairwisePMFSummaParametersFile: "+parametersFileName);
    MeshiLineReader mlrParameters = new MeshiLineReader( parametersFileName );

    StringList lines = new StringList(mlrParameters);
    AtomType type;
    for (String str :lines) {
        StringTokenizer line = new StringTokenizer(str);
        line.nextToken(); // ignore the first word in the line.
        String name = line.nextToken();
        type = AtomType.type(name);
        minAvgSigma[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
    }
    }

}