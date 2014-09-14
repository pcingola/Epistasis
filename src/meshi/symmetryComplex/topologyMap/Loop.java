package meshi.symmetryComplex.topologyMap;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.Command;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:43:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class Loop implements KeyWords {
    private int number;
    private int size;
    private int [] bounds;
    private double angleXfrom, angleXto, angleZfrom, angleZto, variabilityZ;
    private double minWidthOfHairpin, maxWidthOfHairpin, widthOfHairpin;


    /*public static final int[]
N_TER = { 1, 18 },
M1 = { 19, 40 },
E1 = { 41, 75 },
M2 = { 76, 96 },
IN = { 97, 130 },
M3 = { 131, 152 },
E2 = { 153, 188 },
M4 = { 189, 209 },
C_TER = { 210, 999 };*/
//    public static int[] N_TER, C_TER, IN, Ecur;
//    public static int[][] predictedParts, loops, absentParts;
//    public static int numPredictedParts = 4, numLoops = 2, numAbsentParts = 3;
          //

     public Loop(CommandList loop, int number) {
        this.number = number;
        this.bounds = BoundariesMap.loopMap(number);
        size = bounds[bounds.length-1] - bounds[0] + 1;

        Command angle_X =  loop.secondWord(ANGLE_X);
        angleXfrom = angle_X.thirdWordDouble()*Math.PI/180;
        angleXto = angle_X.fourthWordDouble()*Math.PI/180;
        Command angle_Z =  loop.secondWord(ANGLE_Z);

         angleZfrom = (90-angle_Z.thirdWordDouble())*Math.PI/180;
         angleZto = (90-angle_Z.fourthWordDouble())*Math.PI/180;
     //    variabilityZ = angle_Z.fifthWordDouble()*Math.PI/180;
         if (angleZfrom > angleZto) {
             double tmpZ;
             tmpZ = angleZfrom;
             angleZfrom = angleZto;
             angleZto = tmpZ;
         }

        minWidthOfHairpin = loop.secondWord(MIN_WIDTH_OF_HAIRPIN).thirdWordDouble();
        maxWidthOfHairpin = loop.secondWord(MAX_WIDTH_OF_HAIRPIN).thirdWordDouble();
        widthOfHairpin = loop.secondWord(WIDTH_OF_HAIRPIN).thirdWordDouble();
     }
    public int number() {return number;}
    public int size() {return size;}
    public int[] bounds() {return bounds;}
    public double angleXfrom() {return angleXfrom;}
    public double angleXto() {return angleXto;}
    public double angleZfrom() {return angleZfrom;}
    public double angleZto() {return angleZto;}
  //  public double variabilityZ() {return variabilityZ;}
    public double minWidthOfHairpin() {return minWidthOfHairpin;}
    public double maxWidthOfHairpin() {return maxWidthOfHairpin;}
    public double widthOfHairpin() {return widthOfHairpin;}
}
