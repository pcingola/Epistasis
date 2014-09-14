package meshi.geometry;
import  meshi.molecularElements.*;
import  meshi.molecularElements.atoms.*;
import  meshi.parameters.*;
import java.util.*;
import meshi.util.filters.*;
import meshi.util.*;
public class Grid {
    public final double edge;
    private final int minSphereRadiusToEdge;
    private final long  MAX_GRID_SIZE = 1000000;
    private final double SPHERE_COND = 1.3;
    int xSize = 0, ySize = 0, zSize = 0;
    int xSizeNew, ySizeNew, zSizeNew;
    public double minX, minY, minZ, maxX, maxY, maxZ;
    MolecularSystem atoms;
    GridCell[][][] cells;
    GridCell defaultCell;
    double[][] prevCoor;
    ArrayList<AtomCore> movingAtoms;
    int cellCapacity;
    private boolean failToBuild = true;

    public Grid(double edge,double minSphereRadius){
        this.edge = edge;
	//set number of the neighboring cells.
	int temp = (int)(minSphereRadius/edge);
	if (temp == (minSphereRadius/edge))
	    minSphereRadiusToEdge = temp;
	else
	    minSphereRadiusToEdge = temp + 1;
    }

    public Grid(MolecularSystem atoms,double edge,double minSphereRadius)   throws UpdateableException {
	if (atoms.size() == 0) 
	    throw new RuntimeException("No atoms for grid");
	this.atoms = atoms;
        this.edge = edge;
	//set number of the neighboring cells.
	int temp = (int)(minSphereRadius/edge);
	if (temp == (minSphereRadius/edge))
	    minSphereRadiusToEdge = temp;
	else
	    minSphereRadiusToEdge = temp + 1;
        double edge3 = edge*edge*edge;       
	//
	prevCoor = new double[atoms.size()][3];
	for (int i = 0; i < prevCoor.length; i++) {
	    prevCoor[i][0] = prevCoor[i][1] = prevCoor[i][2] = -9999.9999;
	    movingAtoms = new ArrayList<AtomCore>();
	    if (edge <= 2) cellCapacity = 2;
	    else 
		if (edge > 500) cellCapacity = atoms.size();
		else cellCapacity = round(edge3/8);
	}
	build();
	for (int i = 0; i < prevCoor.length; i++) {
	    prevCoor[i][0] = prevCoor[i][1] = prevCoor[i][2] = -9999.9999;
	    movingAtoms = new ArrayList<AtomCore>();
	    if (edge <= 2) cellCapacity = 2;
	    else 
		if (edge > 500) cellCapacity = atoms.size();
		else cellCapacity = round(edge3/8);
	}
    }

    /**
     * Builds the grid. Return false upon failure (say, if the atoms are too spread around 
     * in a space that is too large to be devided into cells without a memory failure.
     **/ 

    public boolean build()  throws UpdateableException {
	double bufferOneThirdSqr = DistanceMatrix.bufferOneThirdSqr;
	minX = minY = minZ = 10000;
	maxX = maxY = maxZ = -10000;
	int size = atoms.size();
			
	movingAtoms.clear();
	for (AtomCore atom:atoms){
	    if (atom.status().activeOrImage())  {
		double x = atom.x(); 
		double y = atom.y(); 
		double z = atom.z();
		if (x < minX) minX = x; 
		if (y < minY) minY = y; 
		if (z < minZ) minZ = z;
		if (x > maxX) maxX = x; 
		if (y > maxY) maxY = y; 
		if (z > maxZ) maxZ = z;
		int atomNumber = atom.number;
		double dx = x - prevCoor[atomNumber][0]; 
		double dy = y - prevCoor[atomNumber][1]; 
		double dz = z - prevCoor[atomNumber][2];
		double d2 = dx*dx + dy*dy + dz*dz;
		if (d2 > bufferOneThirdSqr) { 
		    prevCoor[atomNumber][0] = x; 
		    prevCoor[atomNumber][1] = y; 
		    prevCoor[atomNumber][2] = z;
		    movingAtoms.add(atom);
		}
	    }
	}
	xSizeNew = round((maxX-minX)/edge )+1;
	ySizeNew = round((maxY-minY)/edge )+1;
	zSizeNew = round((maxZ-minZ)/edge )+1;
	if ((xSizeNew > xSize)|
	    (ySizeNew > ySize)|
	    (zSizeNew > zSize)|
	    (xSizeNew < xSize-1)|
	    (ySizeNew < ySize-1)|
	    (zSizeNew < zSize-1)) {
	    xSize = xSizeNew;
	    ySize = ySizeNew;
	    zSize = zSizeNew;
	    if (xSize*ySize*zSize > MAX_GRID_SIZE) {
		System.out.println(" An Error wile creating a new grid\n"+
				   "The requested grid size "+xSize+"*"+ySize+"*"+zSize+"="+
				   (xSize*ySize*zSize)+
				   "is larger than MAX_GRID_SIZE="+MAX_GRID_SIZE+"\n"+
				   minX+" "+minY+" "+minZ+" "+maxX+" "+maxY+" "+maxZ+" ");
		for (int i = 0; i < atoms.size(); i++) System.out.println(atoms.get(i));
		failToBuild = true;
		return false;
	    } 
	    if ((xSize < 0) || (ySize < 0)  || (zSize < 0)) {
		System.out.println("************************************************************\n"+" An Error wile creating a new grid");
		for (AtomCore atom:atoms) System.out.println("#### "+atom);
		throw new RuntimeException("Weird grid dimention "+xSize+" "+ySize+" "+zSize+"\n"+
					    minX+" "+minY+" "+minZ+" "+maxX+" "+maxY+" "+maxZ+" ");
	    }
	    cells = new GridCell[xSize][ySize][zSize];
	    if ((cells.length != xSize )| 
		(cells[0].length != ySize )| 
		(cells[0][0].length != zSize )){
		System.out.println(" An Error wile creating a new grid"+
				   "The program failed in building the requested "
				   +xSize+"*"+ySize+"*"+zSize+" array.\n"+
				   "probably the array is too large. You may solve this problem by "+
				   "reducing MAX_GRID_SIZE which is currently "+MAX_GRID_SIZE);
		throw new UpdateableException();
	    } 
	    for (int X = 0; X < xSize; X++)
		for (int Y = 0; Y < ySize; Y++)
		    for (int Z = 0; Z < zSize; Z++)
		        try {
			    cells[X][Y][Z] = new GridCell(cellCapacity);
			}
                        catch (Exception ex) {
                                System.out.println(ex);
				System.out.println("Very weird exception "+X+" "+Y+" "+Z+" "+cells.length+
						   " "+cells[0].length+" "+cells[0][0].length);
				throw new UpdateableException();
			}

	}
	else 
	    for (int X = 0; X < xSize; X++)
		for (int Y = 0; Y < ySize; Y++)
		    for (int Z = 0; Z < zSize; Z++) {
			try {
			    cells[X][Y][Z].clear();
			}
			catch (RuntimeException ex) {
			    System.out.println("cells.length = "+cells.length+" ; xSize = "+xSize+" ; X = "+X);
			    System.out.println("cells[X].length = "+cells[X].length+" ; ySize = "+ySize+" ; Y = "+Y);
			    System.out.println("cells[X][Y].length = "+cells[X][Y].length+" ; zSize = "+zSize+" ; Z = "+Z);
			    for (Object o:atoms)
				System.out.println(o);
			    throw ex;
			}
		    }

	for (AtomCore atom:movingAtoms){
	    int X = round((atom.x()-minX)/edge );
	    int Y = round((atom.y()-minY)/edge );
	    int Z = round((atom.z()-minZ)/edge );
           
	    int xFrom = ((X < minSphereRadiusToEdge)?0:(X-minSphereRadiusToEdge));
	    int xTo = ((X >= xSize - minSphereRadiusToEdge)?xSize:(X+minSphereRadiusToEdge+1));
	    int yFrom = ((Y < minSphereRadiusToEdge)?0:(Y-minSphereRadiusToEdge));
	    int yTo = ((Y >= ySize - minSphereRadiusToEdge)?ySize:(Y+minSphereRadiusToEdge+1));
	    int zFrom = ((Z < minSphereRadiusToEdge)?0:(Z-minSphereRadiusToEdge));
	    int zTo = ((Z >= zSize - minSphereRadiusToEdge)?zSize:(Z+minSphereRadiusToEdge+1));            
            
	    if    (minSphereRadiusToEdge <= 1)  {
		for (int i = xFrom; i < xTo; i++ )
		    for (int j = yFrom; j < yTo; j++ )
			for (int k = zFrom; k < zTo; k++ ) 
			    cells[i][j][k].add(atom);    		     
	    }
            else  {
                double cellR = (minSphereRadiusToEdge+0.5)*(minSphereRadiusToEdge+0.5);
		for (int i = xFrom; i < xTo; i++ )
                    for (int j = yFrom; j < yTo; j++ ) {                    
			if ( ((i-X)*(i-X)+(j-Y)*(j-Y)) / cellR > SPHERE_COND) continue;
			for (int k = zFrom; k < zTo; k++ ) {
			    if ( ((i-X)*(i-X)+(k-Z)*(k-Z)) / cellR > SPHERE_COND) continue;
			    cells[i][j][k].add(atom);
			}
		    }
            }
	}
	failToBuild = false;
	return true;	
    }
    public GridCell getCell(AtomCore atom){
	int X = round((atom.x()-minX)/edge );
	int Y = round((atom.y()-minY)/edge );
	int Z = round((atom.z()-minZ)/edge );
	return cells[X][Y][Z];
    }

    public static int round(double d){
        return (int)(d+0.5);
    }
    public boolean failToBuild() {return failToBuild;}

}
