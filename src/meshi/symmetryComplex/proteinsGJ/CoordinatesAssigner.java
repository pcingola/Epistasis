package meshi.symmetryComplex.proteinsGJ;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.Command;
import meshi.util.MeshiProgram;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ChainList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.ca.CaResidue;
import meshi.symmetryComplex.topologyMap.Topology;
import meshi.symmetryComplex.topologyMap.Loop;
import meshi.symmetryComplex.utils.GJUtils;
import meshi.symmetryComplex.energy.cylinder.CylinderEnergyElement;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.symmetryComplex.molecularImageElements.ImageAtom;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.geometry.Coordinates;

import java.util.Random;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 15/01/2009
 * Time: 11:24:28
 * To change this template use File | Settings | File Templates.
 */
public class CoordinatesAssigner implements KeyWords {
    private int maxNumberOfClashes;
    public static final double CA_CA_DISTANCE = 3.8;
    public static final double CA_CA_DISTANCE2 = 3.8*3.8;
    public static double DELTA_ANGLE = 10.0/180.0*Math.PI;
    public static double DISTANCE_BetweenLoops;
    public static boolean isDistanceBetweenLoops;

    public final int nLoops = 2;
    private ResidueList residueList;
    private static int currentResidueIndex;
    private int size, sizeMinusOne, sizeMinusTwo;
    private int nTries;
    private double clashDistance;
    private Topology topology;
    private Loop[] loops;
    private final int[] LOOP2_IN_BARREL_RESIDUES = {157, 168, 173, 183};       //these values are actual for GJ Cx32 only!
    private final int[] LOOP1_IN_BARREL_RESIDUES = null;

    private SymmetricComplex gj;
    private double outerBarrelOuterR, outerBarrelInnerR, outerBarrelHeight,
   innerBarrelOuterR, innerBarrelInnerR, innerBarrelHeight;
   static double angleXfrom, angleXto, angleZfrom, angleZto;
   static boolean isSteps;
   static boolean toConstrict;
   static boolean isStrictClashes;
   private List<Transformation> transformations;

public CoordinatesAssigner(SymmetricComplex sc, CommandList commands, ResidueList residueList,
        double clashDistance, int maxNumberOfClashes, int nTries, Topology topology) {
    this(residueList, clashDistance, maxNumberOfClashes, nTries, topology);
    gj = sc;
    transformations = gj.transformations();
    Command innerBarrelCommand = commands.firstWordFilter(CYLINDER_ENERGY).secondWord("innerBarrel");
    innerBarrelOuterR = innerBarrelCommand.thirdWordDouble();
    innerBarrelInnerR = innerBarrelCommand.fourthWordDouble();
    innerBarrelHeight = innerBarrelCommand.fifthWordDouble();
    Command outerBarrelCommand = commands.firstWordFilter(CYLINDER_ENERGY).secondWord("outerBarrel");
    outerBarrelOuterR = outerBarrelCommand.thirdWordDouble();
    outerBarrelInnerR = outerBarrelCommand.fourthWordDouble();
    outerBarrelHeight = outerBarrelCommand.fifthWordDouble();
    isSteps = commands.firstWord(STEPS).secondWord().equals("on");
    if (isSteps) DELTA_ANGLE = commands.firstWord(STEPS).thirdWordDouble()/180.0*Math.PI;
    toConstrict = commands.firstWord(CONSTRICT).secondWord().equals("on");
    isStrictClashes = commands.firstWord(STRICT_CLASHES).secondWord().equals("on");
    isDistanceBetweenLoops = commands.firstWord(CHECK_INTERLOOP_DISTANCE).secondWord().equals("on");
    if (isDistanceBetweenLoops)
        DISTANCE_BetweenLoops = commands.firstWord(CHECK_INTERLOOP_DISTANCE).thirdWordDouble();
}

    public CoordinatesAssigner(ResidueList residueList, double clashDistance, int maxNumberOfClashes, int nTries, Topology topology) {
        this.residueList = residueList;
        this.maxNumberOfClashes = maxNumberOfClashes;
        size = residueList.size();
        if (size < 5)
            throw new RuntimeException("ResidueList must have at least five members.\n");
        sizeMinusOne = size - 1;
        sizeMinusTwo = size - 2;
        this.nTries = nTries;
        this.clashDistance = clashDistance;
        this.topology = topology;
        this.loops = topology.loops();
    }

    public void step() {
        CaResidue currentResidue;
        if (currentResidueIndex >= size)
                currentResidueIndex = 0;
        currentResidue = (CaResidue) residueList.get(currentResidueIndex);
        if (currentResidue.dummy())
                   throw new RuntimeException("Try to assign coordinates for dummy residue "+currentResidue+" It must be nowhere..\n");
            if (currentResidue.ca().nowhere()) {
                tryToGenerateNewAtomLeft(currentResidue);
                if (currentResidue.ca().nowhere()) tryToGenerateNewAtomRight(currentResidue);
                if (!currentResidue.ca().nowhere()) {
                    currentResidueIndex++;
                    if (currentResidueIndex >= size) currentResidueIndex = 0;
                    //gj.buildLocations(); // + 11 transformations for the current residue - to know clashes
                    ChainList ch = gj.chains();
                    ImageAtom ia;
                    for (int i = 1; i< ch.size();i++){
                                 ia = (ImageAtom)ch.get(i).residueAt(currentResidue.number()).ca();
                                 ia.buildLocation();
                    }
                }
            }
        currentResidueIndex++;
    }


    private void tryToGenerateNewAtomLeft( CaResidue currentResidue) {
        Residue from;
        Coordinates coordinates;
        if (currentResidueIndex > 0) {
            from =  residueList.get(currentResidueIndex-1);
            if (!from.ca().nowhere()) {
                coordinates = getAtom((CaResidue) from, currentResidue, 1);
                if (coordinates != null) {
                    currentResidue.ca().setXYZ(coordinates);
                    System.out.println("Coordinates added to "+currentResidue+" based on "+from);
                }
                else {
                    from.ca().setXYZ(new Coordinates());
                    System.out.println("Coordinates removed from  "+from);
    //					if (currentResidueIndex > 1)
    //						((CaResidue) residueList.get(currentResidueIndex-2)).addTrial();
                }
            }
        }
    }

    private void tryToGenerateNewAtomRight( CaResidue currentResidue) {
        Residue from;
        Coordinates coordinates;
        if (currentResidueIndex  < sizeMinusOne) {
            from =  residueList.get(currentResidueIndex+1);
            if (!from.ca().nowhere()) {
                 coordinates = getAtom((CaResidue) from, currentResidue, -1);
                 if (coordinates != null) {
                    currentResidue.ca().setXYZ(coordinates);
                    System.out.println("Coordinates added to "+currentResidue+" based on "+from);
                }
                else {
                    from.ca().setXYZ(new Coordinates());
                    System.out.println("Coordinates removed from  "+from);
    //					if (currentResidueIndex > 1)
    //						((CaResidue) residueList.get(currentResidueIndex-2)).addTrial();
                }
            }
        }
    }


    private Coordinates getAtom(CaResidue from, CaResidue current, int direction){
        Atom atomToFirstNextLoop = null, atomToLastNextLoop = null;
        double dis, dis1 , dis2 , disSum = 0, widthSum = 0;
        int clashes;
        Coordinates coordinates;
        Coordinates[] best = new Coordinates[maxNumberOfClashes+1];
        double[] dbest = new double[maxNumberOfClashes+1];
        double[] d3best = new double[maxNumberOfClashes+1];
        double minWidth, maxWidth, width;

        Residue to = getTo(direction);
        if (to.ca().nowhere())
          throw new RuntimeException("Can't find residueTo to assign coordinates of nowhere atoms.");
        Atom atomTo = ((!to.ca().nowhere())? to.ca():  null);  // atomTo can be null in the common case algorithm but can not be null for GJ extracellular loops

        for (int i = 0; i <= maxNumberOfClashes; i++) {
            dbest[i] = 1000;
            d3best[i] = 1000;
            best[i] = null;
        }

         int loopNum = topology.whichLoop(current.number());
         if (loopNum == 0) {   //toDo to check
             int nAtomsforOneLoop = size/nLoops;
             int partOfchain = from.number()/nAtomsforOneLoop;
             loopNum = partOfchain+1;
         }

          Loop loop = loops[loopNum-1];
          if (loopNum != loop.number())
            throw new RuntimeException("\nOrenWindow.getAtom something weird in the loops assignment\n");
          minWidth = loop.minWidthOfHairpin();
          maxWidth = loop.maxWidthOfHairpin();
          width = loop.widthOfHairpin();
          angleXfrom = loop.angleXfrom();
          angleXto = loop.angleXto();
          angleZfrom = loop.angleZfrom();
          angleZto = loop.angleZto();

        Residue toFirstNextLoopDefined, toLastNextLoopDefined;

        if (isDistanceBetweenLoops){

            int nextLoop = ((loopNum == nLoops) ? 1 : (loopNum+1));
            int nextLoopFirstAtom = topology.loopFirstAtomNumber(nextLoop);
            int nextLoopLastAtom = topology.loopLastAtomNumber(nextLoop);

            toFirstNextLoopDefined = getToNull(nextLoopFirstAtom-2,nextLoopLastAtom+1,1);// == null of the end of the loop fullfilling
            toLastNextLoopDefined = getTo(toFirstNextLoopDefined, nextLoopLastAtom+1, 1);  // == null of the end of the loop fullfilling

            if (toFirstNextLoopDefined != null && toLastNextLoopDefined != null &&
                (!topology.inLoop(toFirstNextLoopDefined.number()+1)
                || !topology.inLoop(toLastNextLoopDefined.number()-1)) )
                throw new RuntimeException("\ndisallignment in the next loop searching algorithm\n");

                     atomToFirstNextLoop =((toFirstNextLoopDefined != null)? toFirstNextLoopDefined.ca():  null);
                     atomToLastNextLoop =((toLastNextLoopDefined != null)? toLastNextLoopDefined.ca():  null);

            widthSum = width+2*DISTANCE_BetweenLoops;
           }

        for (int i = 0; i < nTries; i++) {

             coordinates = getAtom(current, from, loop);
             dis = ((atomTo != null)? GJUtils.quickDistanceFrom(coordinates, new Coordinates(atomTo)):100);
             clashes = getClashes(coordinates, current, gj.residues(), clashDistance);                     //toDo check gj.residues if updated

             boolean atomInCylinder = true;
             if (loopNum == 1)
                 atomInCylinder = CylinderEnergyElement.isPointInCylinder(coordinates, outerBarrelOuterR, outerBarrelInnerR, outerBarrelHeight);
             if (loopNum == 2)
                 atomInCylinder = CylinderEnergyElement.isPointInCylinder(coordinates, innerBarrelOuterR, innerBarrelInnerR, innerBarrelHeight);
             if (! atomInCylinder)
                 clashes += maxNumberOfClashes;

             if ( ((current.number() >= LOOP2_IN_BARREL_RESIDUES[0]) &&(current.number() <= LOOP2_IN_BARREL_RESIDUES[1] ) ||
                      ((current.number() >= LOOP2_IN_BARREL_RESIDUES[2] )&&(current.number() <= LOOP2_IN_BARREL_RESIDUES[3] )) )
                      && (clashes == maxNumberOfClashes))
                 continue;             //these residues are strictly recomended to be in barrel

                 if (clashes <= maxNumberOfClashes) {
                 boolean changeBestDis;
                 if (isDistanceBetweenLoops & (atomToFirstNextLoop != null) & (atomToLastNextLoop != null)){
                   changeBestDis = false;
                   if (Math.abs(dbest[clashes] - width) > Math.abs(dis - width)){
                    dis1 = GJUtils.quickDistanceFrom(coordinates, new Coordinates(atomToFirstNextLoop));
                    dis2 = GJUtils.quickDistanceFrom(coordinates, new Coordinates(atomToLastNextLoop));
                    disSum = dis + dis1 + dis2;
                    changeBestDis = Math.abs(d3best[clashes] - widthSum) > Math.abs(disSum - widthSum);
                   }
                 }
                 else  {
                         if ((dis < minWidth && dis > maxWidth) &&// not fit
                             (dbest[clashes] >= minWidth && dbest[clashes] <= maxWidth))
                             continue;
                         if (toConstrict){
                         //changeBestDis = (Math.abs(dbest[clashes] - width) > Math.abs(dis - width))&& (dbest[clashes] > dis);
                           if (dbest[clashes] < minWidth)
                             changeBestDis = (Math.abs(dbest[clashes] - width) > Math.abs(dis - width));
                           else changeBestDis = (Math.abs(dbest[clashes] - width) > Math.abs(dis - width))
                                   && (dbest[clashes] > dis);
                         }
                         else changeBestDis = Math.abs(dbest[clashes] - width) > Math.abs(dis - width);
                        }

                 if (changeBestDis){
                         dbest[clashes] = dis;
                         d3best[clashes] = disSum;
                         best[clashes] = coordinates;
                         }
                 }
         }

         int num_best_1 = -1;
         double val_min_1;

         do {
             num_best_1++;
             if (num_best_1 > maxNumberOfClashes)
                        return null;
             val_min_1 = dbest[num_best_1];
         }  while (dbest[num_best_1] == 1000 );


         if (isStrictClashes) return best[num_best_1];

         if (val_min_1 >= minWidth && val_min_1 <= maxWidth)
             return best[num_best_1];

         int num_best_2 = 0;
         double val_min_2 = dbest[0];
             for (int i = 1; i <= maxNumberOfClashes; i++) {
                 if (Math.abs(dbest[i] - width) < Math.abs(val_min_2 - width)) {
             //if (dbest[i] < val_min_2) {
                    val_min_2 = dbest[i];
                    num_best_2 = i;
             }
         }

         if (val_min_2 <= minWidth || val_min_2 >= maxWidth)
             return best[num_best_1];

         return best[num_best_2];
     }


    private Residue getTo(int direction) {
        int i = currentResidueIndex + direction;
        while ((i >= 0) & (i <= sizeMinusOne)) {
            if (!residueList.get(i).ca().nowhere())
                return residueList.get(i);
            i+=direction;
        }
        return null;
    }

    private Residue getTo(Residue from, int lastAtom, int direction) {
        if (from ==  null) return null;
        int i = from.number() + direction;
        while ((i >= 0) & (i <= lastAtom)) {
            if (residueList.get(i).ca() != null)
                return residueList.get(i);
            i+=direction;
        }
        return null;
    }

    private Residue getToNull(int from, int lastAtom, int direction) {
        int i = from + direction;
        while ((i >= 0) & (i <= lastAtom)) {
            if (residueList.get(i).ca() == null)
                return residueList.get(i-1);
            i+=direction;
        }
        return null;
    }

    private Coordinates getAtom(CaResidue current, CaResidue from, Loop loop){
        double[] pointCoors;
        double[] fromXYZ = {from.ca().x(), from.ca().y(), from.ca().z()};
        //from2XYZ = {from2.ca().x(), from2.ca().y(), from2.ca().z()};

        if (angleXfrom == 0 && angleXto == 0) {
            int direction;
            if (loop.number()==1) direction = -1;
            else direction = 1;
            pointCoors = getRandomPointAroundCyllinder(from.ca(), CA_CA_DISTANCE, direction, fromXYZ);
        }
        else {
           //if (isTransform)
             // pointCoors = getRandomPointOnSphereBiasDirection(CA_CA_DISTANCE, fromXYZ, from2XYZ);
              //else
                 pointCoors = getRandomPointOnSphere(CA_CA_DISTANCE);
        }


            return new Coordinates (from.ca().x()+pointCoors[0], from.ca().y()+pointCoors[1], from.ca().z()+pointCoors[2]);
    }

    private static double[] getRandomPointOnSphereBiasDirection(double radius, double[] fromXYZ, double[] from2XYZ) {
        Random rnd = MeshiProgram.randomNumberGenerator();
        //double azimuth = (2*rnd.nextDouble() - 1) * Math.PI;
        //double angle = rnd.nextDouble() * Math.PI/2;
        double azimuth = angleXfrom + (angleXto - angleXfrom)*rnd.nextDouble(); // {-Pi,Pi} - in the plane XOY
        // if angleXfrom == angleXto than angle is constant, angle = angleXfrom
        double angle = angleZfrom + (angleZto - angleZfrom)*rnd.nextDouble();   // {0, Pi/2}  - along OZ
        if (isSteps)   //the "accordion-like" pattern for Ca
                if (currentResidueIndex%2 ==0) angle += DELTA_ANGLE;
                else angle -= DELTA_ANGLE;

        double[] xyz = new double[3];
        double cosAngle = Math.cos(angle);
        xyz[2] = Math.sin(angle);
        xyz[0] = Math.cos(azimuth) * cosAngle;
        xyz[1] = Math.sin(azimuth) * cosAngle;
        // for clashes == symmetry loops
                double[] fromsVector = {fromXYZ[0]-from2XYZ[0],fromXYZ[1]-from2XYZ[1], fromXYZ[2]-from2XYZ[2]};
                normalize(fromsVector);
                double[] perpendicularVector = {-fromsVector[1], fromsVector[0], 0.0};
                double angle2 = Math.asin(magnitude(perpendicularVector));//problem with asin
                normalize(perpendicularVector);
                Transformation trans = new Transformation(perpendicularVector[0], perpendicularVector[1], perpendicularVector[2], angle2);
                xyz = trans.transform(xyz);
     return new double[] {xyz[0] * radius, xyz[1] * radius, xyz[2] * radius};
    }

    private static void normalize(double[] fromsVector) {
        double magnitude = magnitude(fromsVector);
        fromsVector[0] /= magnitude;
        fromsVector[1] /= magnitude;
        fromsVector[2] /= magnitude;
    }


    private static double magnitude(double[] vector) {
        return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    }

    /**
     * Returns a point on a sphere (more precisely, on a shell) of a given
     * radius. Azimuth is uniformly chosen from the range of [-PI,+PI]; angle is
     * chosen from the range of [-PI/2,+PI/2], but uniformly sampling the angle
     * causes bias to the poles. Instead, one call to ArcCos.acos is used to
     * uniformly sample a point on the sphere. Transforming the azimuth and
     * angle to cartesian coordinates takes two calls to Math.sin and Math.cos,
     * each.
     * <P>
     * For further details on uniformly sampling points on a
     * shell see, for example, Kuffner, Effective Sampling and Distance Metrics
     * for 3D Rigid Body Path Planning, ICRA 2004.
     */
    private static double[] getRandomPointOnSphere(double radius) {
        Random rnd = MeshiProgram.randomNumberGenerator();
        //double azimuth = (2*rnd.nextDouble() - 1) * Math.PI;
        //double angle = ArcCos.acos(2*rnd.nextDouble() - 1) - Math.PI/2;
        double azimuth = angleXfrom + (angleXto - angleXfrom)*rnd.nextDouble();// {-Pi,Pi} - in the plane XOY
        // if angleXfrom == angleXto than angle is constant, angle = angleXfrom
        //double var = 5*Math.PI/180;
        double angle = angleZfrom + (angleZto - angleZfrom)*rnd.nextDouble();  // {0, Pi/2}  - along OZ
        if (isSteps)   //the "accordion-like" pattern for Ca
                if (currentResidueIndex%2 ==0) angle += DELTA_ANGLE;
                else angle -= DELTA_ANGLE;
        double x, y, z, cosAngle = Math.cos(angle);
        z = Math.sin(angle);
        x = Math.cos(azimuth) * cosAngle;
        y = Math.sin(azimuth) * cosAngle;
        return new double[] {x * radius, y * radius, z * radius};
    }


    private static double[] getRandomPointAroundCyllinder(Atom from, double radius, int direction, double[] fromXYZ) {
        Random rnd = MeshiProgram.randomNumberGenerator();
        double angle = angleZfrom + (angleZto - angleZfrom)*rnd.nextDouble();   // {0, Pi/2}  - along OZ
        if (isSteps)   //the "accordion-like" pattern for Ca
                if (currentResidueIndex%2 ==0) angle += DELTA_ANGLE;
                else angle -= DELTA_ANGLE;

        double x0 = fromXYZ[0];
        double y0 = fromXYZ[1];
        double radius0 = Math.sqrt(x0*x0+y0*y0);
        double radiusCa = CA_CA_DISTANCE*Math.cos(angle);
        double prevAzimuth = Math.acos(x0/radius0);
        double deltaAzimuth = 2*Math.asin(radiusCa/(2*radius0));
        double azimuth = prevAzimuth + direction*deltaAzimuth; //XOY

        double[] xyz = new double[3];
        xyz[2] = Math.sin(angle);
        double SPRED_MAX_XOY = 1;
        xyz[0] = radius0*Math.cos(azimuth) + rnd.nextDouble()*SPRED_MAX_XOY;
        xyz[1] = radius0*Math.sin(azimuth) + rnd.nextDouble()*SPRED_MAX_XOY;

       return new double[] {xyz[0], xyz[1], from.z()+xyz[2] * radius};
    }


    /**
     * Modified by Oren Wolfshtat.
     */
    public int getClashes(Coordinates coordinates, Residue current, ResidueList residueList, double clashDistance) {
        if (coordinates.nowhere())
            throw new RuntimeException("Tries to find clashed for nowere coordinates.");


        double dis;
        int clashes = 0;
        Atom ca;

        for (Residue residue: residueList) {
                            if (residue.dummy()) continue;              //toDo 
            //throw new RuntimeException("Dummy residue in the  ResidueList"+residue);

            ca = residue.ca();
            if (ca == null)
                            continue;                            //toDo if I need this
            if ((!ca.nowhere()) && (ca.x() != coordinates.x())&& (ca.y() != coordinates.y())&& (ca.z() != coordinates.z())) {
                dis = GJUtils.quickDistanceFrom(coordinates, new Coordinates(ca));

                if ((Math.abs(residue.number()-current.number()) == 2) && (dis < 6))
                    clashes++;
                else
                    if (dis < clashDistance) {
                        clashes++;
                    }
            }
        }
        clashes += clashesWithImages(coordinates, clashDistance);
        return clashes;
    }

    // ------------------------------------------------------------------------
    public int clashesWithImages(Coordinates coordinates, double clashDist) {
        int clashes = 0;
        double[][] xyz = new double[4][1];
        Coordinates tCoordinates;
        for (Transformation t : transformations) {
                xyz[0][0] = coordinates.x();
                xyz[1][0] = coordinates.y();
                xyz[2][0] = coordinates.z();
                xyz[3][0] = 1.0;
                xyz = t.transform(xyz);
                tCoordinates = new Coordinates(xyz[0][0], xyz[1][0], xyz[2][0]);
                if (GJUtils.quickDistanceFrom(coordinates, tCoordinates) < clashDist)
                            clashes++;
        }
        return clashes;
   }

}



