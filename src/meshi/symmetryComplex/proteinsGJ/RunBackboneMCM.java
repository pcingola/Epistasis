package meshi.symmetryComplex.proteinsGJ;

import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.symmetryComplex.topologyMap.BoundariesMap;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 20/02/2008
 * Time: 14:54:41
 * To change this template use File | Settings | File Templates.
 */
public class RunBackboneMCM {
    public static final String NAME = "RunGJBackboneMCM";

    public static void main(String[] args) throws IOException {
        CommandList  commands;

        commands = Utils.init(args,5, Integer.parseInt(args[4]),  // args[0] instead args[2]
                "Usage: java -XmxNNNm "+NAME+"  <commands file name><pdb file for TM helices> <sequence><output(path/fileNameMask><seed>\n\n"+
                "NNN is the size of the expected memory requirement in MegaBytes.");
        int seed = Integer.parseInt(args[4]);
        System.out.print(NAME+"\nSeed =\t"+seed+"\n");

        String tmCaLocation = args[1],
            sequence = args[2],
            outputFileMask = args[3];


        BoundariesMap boundsMap = new BoundariesMap(commands); //set static data from file

        BackboneMCM backboneMCM = new BackboneMCM (
                    tmCaLocation,
                    sequence,
                    commands,
                    outputFileMask
            );
    }

}
