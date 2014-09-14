package meshi.symmetryComplex.utils;



import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.ResidueList;
import meshi.geometry.Coordinates;
import meshi.util.filters.Filter;
import meshi.symmetryComplex.topologyMap.BoundariesMap;

import java.util.Arrays;
import java.lang.Object;



/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Mar 3, 2008
 * Time: 10:36:02 AM
 * To change this template use File | Settings | File Templates.
 */
public class LoopUtils  {

    public static ResidueList residuesOfLoops(ResidueList residues) {
        return residues.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap()));
    }

    public static ResidueList residuesOfLoopNumber(ResidueList residues, int i) {
          return residues.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(i)));
    }

      public static AtomList atomsOfLoops(AtomList atoms) {
          return atoms.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap()));
      }

      public static AtomList atomsOfLoopNumber(AtomList atoms, int i) {
            return atoms.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(i)));
      }

      public static AtomList atomsCAOfLoops(AtomList atoms) {
          return (atoms.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap()))).CAFilter();
      }

      public static AtomList atomsCAOfLoopNumber(AtomList atoms, int i) {
            return (atoms.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(i)))).CAFilter();
      }


    private static class Line{
        double A, B, C;
        public Line(){}
        public void setLine(Coordinates c1,Coordinates c2){
            A = c2.y() - c1.y();
            B = c1.x() - c2.x();
            C = -(c1.y()*B+c1.x()*A);
        }

        public void setLine(double c1x, double c1y, double c2x, double c2y){
            A = c2y - c1y;
            B = c1x - c2x;
            C = -(c1y*B+c1x*A);
        }
    }

    private static class Point{
        double x, y;
        public Point(double x, double y){
            this.x = x;
            this.y = y;
        }
        double x(){return x;}
        double y(){return y;}
    }


    public static Point pointOfCrossing(double c1x, double c1y, double c2x, double c2y,
                                  double c3x, double c3y, double c4x, double c4y){
        Line line1 = new Line();
        Line line2 = new Line();
        line1.setLine(c1x, c1y, c2x, c2y);
        line2.setLine(c3x, c3y, c4x, c4y);
//bounds of projections
        double minX, maxX, minY, maxY;
//point of lines' crossing
        double pointX, pointY;

        double tmp = line1.A*line2.B-line2.A*line1.B;
        if (tmp != 0) {
            pointX = (line1.B*line2.C-line2.B*line1.C)/tmp;
            pointY = (line1.C*line2.A-line2.C*line1.A)/tmp;
        }
        else {
            pointX = pointY = 10e6;
            }

        if ( (pointX > c1x && pointX > c2x) || (pointX < c1x && pointX < c2x) ||
             (pointX > c3x && pointX > c4x) || (pointX < c3x && pointX < c4x))
             return null;  // not crossing lins
        if ( (pointY > c1y && pointY > c2y) || (pointY < c1y && pointY < c2y) ||
             (pointY > c3y && pointY > c4y) || (pointY < c3y && pointY < c4y))
             return null;  // not crossing lins


        double [] x = {c1x,c2x,c3x,c4x};
        double [] y = {c1y,c2y,c3y,c4y};

        Arrays.sort(x);
        minX = x[0];
        maxX = x[x.length-1];
        Arrays.sort(y);
        minY = y[0];
        maxY = y[y.length-1];


        Point point = null;
        if (pointX >= minX && pointX <= maxX && pointY >= minY && pointY <= maxY) {
         point = new Point(pointX, pointY);
        }
        return point;
    }


    public static int numberOfTwists (AtomList allAtoms, int loopNumber){
        AtomList atoms = allAtoms.getSourceAtoms();
        Atom[] loopArray = atomsCAOfLoopNumber(atoms,loopNumber).toArrayOfAtoms();
        int numberOfTwists = 0;
        Atom atom1, atom2, atom3, atom4;
        int lastNum = loopArray.length-1;
        double dis1, dis2, minDis;
        int numOfClosestPair, step = 2;

    //    if (ifpointOfCrossing(-4,-4,2,2,2,0,0,2)
      //  &&
        //        ifpointOfCrossing(-3,-3,3,3,3,0,0,3))
        //numberOfTwists++;

        for (int i = 0; i< lastNum/2 - step; i++) {
            atom1 = loopArray[i];
            atom2 = loopArray[i+step];

            minDis = 10000; numOfClosestPair = -1;
            for (int j = 0; j<lastNum-step; j++) {
                if (Math.abs(i-j)<4) continue;
                atom3 = loopArray[j];
                atom4 = loopArray[j+step];

                dis1 = GJUtils.quickDistanceFrom(new Coordinates(atom3), new Coordinates(atom1))+                     //5.29
                       GJUtils.quickDistanceFrom(new Coordinates(atom3), new Coordinates(atom2));
                dis2 = GJUtils.quickDistanceFrom(new Coordinates(atom4), new Coordinates(atom1))+
                       GJUtils.quickDistanceFrom(new Coordinates(atom4), new Coordinates(atom2));
                if (dis1+dis2 < minDis) {
                    numOfClosestPair = j;
                    minDis = dis1+dis2;
                }
             }

            atom3 = loopArray[numOfClosestPair];
            atom4 = loopArray[numOfClosestPair+step];

            if (ifCrossing(atom1,atom2,atom3,atom4,false))
                numberOfTwists++;
        }

       return numberOfTwists;
    }

    private static boolean ifCrossing(Atom atom1, Atom atom2, Atom atom3, Atom atom4,
                                      boolean full_checking){
    Point p1,p2,p3;

    p1 = pointOfCrossing(atom1.x(),atom1.z(), atom2.x(),atom2.z(),
            atom3.x(),atom3.z(), atom4.x(),atom4.z()); //5.29
    if (p1 == null) return false;

    p2 = pointOfCrossing(atom1.y(),atom1.z(),atom2.y(),atom2.z(),
            atom3.y(),atom3.z(),atom4.y(),atom4.z());
    if (p2 == null) return false;


    if (full_checking) {
            p3 = pointOfCrossing(atom1.x(),atom1.y(),atom2.x(),atom2.y(),
                               atom3.x(),atom3.y(),atom4.x(),atom4.y());
            if (p3 == null) return false;
    }
/*
int cross = 0;
if (Math.abs(p1.x()-p2.x()) < 1) cross++;
if (Math.abs(p1.y()-p3.y()) < 1) cross++;
if (Math.abs(p2.y()-p3.x()) < 1) cross++;

if (cross > 1) numberOfTwists++;
//            if (cross == 3) numberOfTwists++;
*/
//            if (Math.abs(p1.x() - p2.x()) < 1.4)

    return true;
    }


    public static int numberOfChainCrosses(AtomList allAtoms, String chain1, String chain2){
           AtomList atoms = allAtoms.CAFilter();
           Atom[] chain1Array = atomsCAOfLoops(atoms.getChainAtoms(chain1)).toArrayOfAtoms();
           Atom[] chain2Array = atomsCAOfLoops(atoms.getChainAtoms(chain2)).toArrayOfAtoms();
        int numberOfCrosses = 0;
        Atom atom1, atom2, atom3, atom4;
        //int lastNum = loopArray.length-1;
        int numOfClosest, step = 2;

        for (int i = 0; i<chain1Array.length-step; i++) {
            atom1 = chain1Array[i];
            atom2 = chain1Array[i+step];
            if(atom2.residueNumber()-atom1.residueNumber()>step)
                continue; //two different loops

            numOfClosest = numOfClosestPairOfAtoms(atom1, atom2, chain2Array, step, false);

            atom3 = chain2Array[numOfClosest];
            atom4 = chain2Array[numOfClosest+step];

            if (ifCrossing(atom1,atom2,atom3,atom4,true))    {
                numberOfCrosses++;
                //System.out.print(atom1+"\n"+atom2+"\n"+atom3+"\n"+atom4+"\n");
            }

        }

       return numberOfCrosses;

    }

    private static int numOfClosestPairOfAtoms(Atom atom1, Atom atom2, Atom[] atomArray, int step, boolean sameArray){                   //5.29
        int num = -1;
        double dis1, dis2, minDis = 10000;
        Atom atom3, atom4;

            for (int j = 0; j<atomArray.length-step; j++) {
                atom3 = atomArray[j];
                atom4 = atomArray[j+step];
                if (!sameArray)
                    if(atom4.residueNumber()-atom3.residueNumber()>step) continue; //two different loops

                dis1 = GJUtils.quickDistanceFrom(new Coordinates(atom3), new Coordinates(atom1))+  //5.29
                       GJUtils.quickDistanceFrom(new Coordinates(atom3), new Coordinates(atom2));
                dis2 = GJUtils.quickDistanceFrom(new Coordinates(atom4), new Coordinates(atom1))+
                       GJUtils.quickDistanceFrom(new Coordinates(atom4), new Coordinates(atom2));

                if (dis1+dis2 < minDis) {
                    num = j;
                    minDis = dis1+dis2;
                }
             }

        return num;
    }


    public static double getAverageAngleForLoop(AtomList allAtoms, int loopNumber) {
        AtomList atoms = allAtoms.getSourceAtoms();
        final Atom[] loopArray = atomsCAOfLoopNumber(atoms,loopNumber).toArrayOfAtoms();

        double sumAngles = 0.0;
        Atom atom1, atom2;
        for (int i = 0; i<loopArray.length-1; i++) {
            atom1 = loopArray[i];
            atom2 = loopArray[i+1];
            sumAngles += getAngle(new Coordinates(atom1),new Coordinates(atom2));
        }
        return sumAngles / (loopArray.length-1);
    }

    public static double getAngle(Coordinates coord1, Coordinates coord2) {
        double dis = GJUtils.quickDistanceFrom(coord1, coord2);
        return Math.acos(Math.abs(coord1.z()-coord2.z())/dis)*180/Math.PI;
    }

    public static double [] getAverageAngleForLoopParts(AtomList allAtoms, int loopNumber, int numOfParts) {
        AtomList atoms = allAtoms.getSourceAtoms();
        double [] part = new double [numOfParts];
        final Atom[] loopArray = atomsCAOfLoopNumber(atoms,loopNumber).toArrayOfAtoms();
        int lastNumber = loopArray.length-1;
        int partSize = Math.round(loopArray.length /(numOfParts*2) );

        Atom atom1, atom2;        //5.29

        for (int p = 0; p < numOfParts-1; p++) {
            for (int i = 0; i< partSize; i++) {
                atom1 = loopArray[i+p*partSize];
                atom2 = loopArray[i+p*partSize+1];
                part[p] +=  getAngle(new Coordinates(atom1), new Coordinates(atom2));
                atom1 = loopArray[lastNumber-p*partSize-i];
                atom2 = loopArray[lastNumber-1-p*partSize-i];
                part[p] +=  getAngle(new Coordinates(atom1), new Coordinates(atom2));
            }
        }
        //last part
        int p = numOfParts - 1;
        int sizeOfLastPart = 0;
        for (int i = 0; i+p*partSize < lastNumber/2; i++){
                atom1 = loopArray[i+p*partSize];
                atom2 = loopArray[i+p*partSize+1];
                part[p] +=  getAngle(new Coordinates(atom1), new Coordinates(atom2));
                atom1 = loopArray[lastNumber-p*partSize-i];
                atom2 = loopArray[lastNumber-1-p*partSize-i];
                part[p] +=  getAngle(new Coordinates(atom1), new Coordinates(atom2));
                sizeOfLastPart ++;
            }


        for (int i = 0; i < numOfParts-1; i++)
               part[i] /= (2*partSize);
        part[numOfParts-1] /= (2*sizeOfLastPart);

        return part;
    }

    /**
     * Returns the average distance between C-Beta atoms of cysteines which are
     * thought to be bonded by disulfide links.
     */
    public static double getAverageDistanceBetweenCysCBconnexin(final String atomType, final
 AtomList atoms) {
        final Atom[] cBetas = atoms.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }
        }).toArrayOfAtoms();
        double ans = 0.0;
        ans += cBetas[0].distanceFrom(cBetas[5]);
        ans += cBetas[1].distanceFrom(cBetas[4]);
        ans += cBetas[2].distanceFrom(cBetas[3]);
        return ans / 3;
    }

    /**
     * Returns distances between C-Beta atoms of some special cysteines pattern for connexin
     */
    public static double getDistancesBetweenCysCBofConnexin(final String atomType, final AtomList
 atoms) {
        final Atom[] cBetas = atoms.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }
        }).toArrayOfAtoms();
        double dis1, dis2, dis3, dis4;

        dis1 = cBetas[0].distanceFrom(cBetas[4]);  //53-173
        dis2 = cBetas[2].distanceFrom(cBetas[3]); // 64-168
        dis3 = cBetas[0].distanceFrom(cBetas[2]); //53-64
        dis4 = cBetas[3].distanceFrom(cBetas[4]); //168-173


        System.out.println("CysConnexin\t53-173\t64-168\tSum\t53-64\t168-173\tSum");

 System.out.println("CysConnexin"+atomType+"\t"+dis1+"\t"+dis2+"\t"+(dis1+dis2)+"\t"+dis3+"\t"+dis4+"\t"+(dis3+dis4));
        return dis1+dis2;
    }

    /**
     * Returns distances between C-Beta atoms of some special cysteines pattern
     */
    public static double getDistancesBetweenCysCBofConnexons(final String atomType,
        final AtomList atoms1, final AtomList atoms2, final AtomList atoms3, final AtomList atoms4)
 {

        final Atom[] cysCB1 = atoms1.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }}).toArrayOfAtoms();

        final Atom[] cysCB2 = atoms2.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }}).toArrayOfAtoms();

        final Atom[] cysCB3 = atoms3.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }}).toArrayOfAtoms();

        final Atom[] cysCB4 = atoms4.splitToChains()[0].filter(new Filter() {
            public boolean accept(Object obj) {
                return ((Atom) obj).nameIs(atomType)
                        & ((Atom) obj).residueName().equalsIgnoreCase("CYS");
            }}).toArrayOfAtoms();

        double dis1, dis2, dis3, dis4, dis5, dis6, sum12;

        dis1 = cysCB1[0].distanceFrom(cysCB2[4]);
        dis2 = cysCB1[2].distanceFrom(cysCB2[3]);
        sum12 = dis1 + dis2;
        dis3 = cysCB1[0].distanceFrom(cysCB1[2]); //53-64
        dis4 = cysCB1[3].distanceFrom(cysCB1[4]); //168-173
        dis5 = cysCB1[5].distanceFrom(cysCB3[5]);
        dis6 = cysCB1[5].distanceFrom(cysCB4[5]);



 System.out.println("CysConnexon\t53A-173B\t64A-168B\tSum\t53-64\t168-173\tSum\t179A-179G\t179A-179H\tSumG\tSumH");

 System.out.println("CysConnexon"+atomType+"\t"+dis1+"\t"+dis2+"\t"+sum12+"\t"+dis3+"\t"+dis4+"\t"+(dis3+dis4));
 System.out.println("CysConnexon179"+atomType+"\t"+dis5+"\t"+dis6+"\t"+(sum12+dis5)+"\t"+(sum12+dis6));
        return sum12;
    }

 }

