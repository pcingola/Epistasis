package meshi.symmetryComplex.proteinsGJ;

import java.io.IOException;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.symmetryComplex.utils.GJFilters;
import meshi.symmetryComplex.topologyMap.BoundariesMap;

public class RunBackboneCompleterAndMinimizer {
    public static final String NAME = "RunBackboneCompleterAndMinimizer";

    public static void main(String[] args) throws IOException {
        CommandList  commands;

        commands = Utils.init(args,5, Integer.parseInt(args[4]),  // args[0] instead args[2]
                "Usage: java -XmxNNNm "+NAME+"  <commands file name><pdb file for TM helices> <sequence><output(path/fileNameMask><seed>\n\n"+
                "For '-XmxNNNm' NNN is the size of the expected memory requirement in MegaBytes.");
        int seed = Integer.parseInt(args[4]);
        System.out.print(NAME+"\nSeed =\t"+seed+"\n");

        String tmCaLocation = args[1],
            sequence = args[2],
            outputFileMask = args[3];

        BoundariesMap boundsMap = new BoundariesMap(commands); //set static data from file

        BackboneCompleterAndMinimizer backboneCompleterAndMinimizer = new BackboneCompleterAndMinimizer (
                    tmCaLocation,
                    sequence,
                    commands,
                    outputFileMask+seed,
                    true
            );
    }
}
