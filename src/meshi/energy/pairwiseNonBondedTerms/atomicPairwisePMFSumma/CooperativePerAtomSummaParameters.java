package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.parameters.*;
import meshi.util.string.*;
import meshi.util.*;
import meshi.util.file.*;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 01/06/2009
 * Time: 12:47:23
 * To change this template use File | Settings | File Templates.
 */

public class CooperativePerAtomSummaParameters implements KeyWords, MeshiPotential {
    protected final double[] mean = new double[AtomType.values().length];
    protected final double[] std = new double[AtomType.values().length];

    public CooperativePerAtomSummaParameters(CommandList commands ) {
    String          parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY ).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_PERATOM_SUMMA_FILENAME).secondWord();

    System.out.println("CooperativePerAtomSummaParametersFile: "+parametersFileName);
    MeshiLineReader mlrParameters = new MeshiLineReader( parametersFileName );

    StringList lines = new StringList(mlrParameters);
    AtomType type;
    for (String str :lines) {
        StringTokenizer line = new StringTokenizer(str);
        //line.nextToken(); // ignore the first word in the line.
        String name = line.nextToken();
        type = AtomType.type(name);
        mean[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
        std[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
    }
    }

}
