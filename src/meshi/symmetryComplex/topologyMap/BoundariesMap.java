package meshi.symmetryComplex.topologyMap;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.CommandsException;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:37:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class BoundariesMap implements KeyWords {
  public static int numberOfChains;
  public static int halfNumberOfChains;
  private static int [][] frozenMap = null,
          loopMap = null,
          absentMap = null;

  public BoundariesMap(CommandList commands) {
      this(commands.firstWord(TOPOLOGY_MAP).secondWord());
  }

    public BoundariesMap(String mapFile) {
        CommandList mapCommands = new CommandList(mapFile,new CommandsException("Error in BoundariesMap"));
        numberOfChains = mapCommands.firstWord(NUMBER_OF_CHAINS).secondWordInt();
        halfNumberOfChains = numberOfChains/2;

        frozenMap = map(mapCommands.firstWordFilter("frozen"), "F");
        absentMap = map(mapCommands.firstWordFilter("absent"), "A");
        loopMap = map(mapCommands.firstWordFilter("loop"), "L");
    }

      public static int [][] map(CommandList commands, String symbolOfPart){
      int size = commands.secondWord(NUMBER).thirdWordPositiveInt();
      int [][] map = new int [size][2];
          for (int i = 0; i< size;i++) {
          String wordToSearch = symbolOfPart+String.valueOf(i+1);
          if (wordToSearch.equals(""))
              throw new RuntimeException("Wrong size of some part of "+symbolOfPart+" is missed in the topology file");
          map[i][0] = commands.secondWord(wordToSearch).thirdWordPositiveInt();
          map[i][1] = commands.secondWord(wordToSearch).fourthWordPositiveInt();
          }
          return map;
      }

    public static int [][] frozenMap(){return frozenMap;}
    public static int [][] absentMap(){return absentMap;}
    public static int [][] loopMap(){return loopMap;}

    public static int [] frozenMap(int i){return frozenMap[i-1];}
    public static int [] absentMap(int i){return absentMap[i-1];}
    public static int [] loopMap(int i){return loopMap[i-1];}

}


