package meshi.geometry;
import meshi.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class Coordinates{
   public static double NOWHERE_CONST = -999.999;
    protected double[] x,y,z;
    private boolean random;
    private boolean nowhere;
    public  Coordinates() {
	this(-999.999, -999.999, -999.999);
	nowhere = true;
    }
    public Coordinates(Atom atom) {
	this();
	if ((atom != null) && (!atom.nowhere())) {
	    this.x[0] = atom.x();	
	    this.y[0] = atom.y();	
	    this.z[0] = atom.z();
	    nowhere = atom.nowhere();
	}
    }

    public  Coordinates(double x, double y, double z) {
	this.x = new double[2];
	this.y = new double[2];
	this.z = new double[2];
	this.x[0] = x;	
	this.y[0] = y;	
	this.z[0] = z;
	this.x[1] = 0.0;
	this.y[1] = 0.0;
	this.z[1] = 0.0;
	random = false;
	nowhere = false;
    }
    public  Coordinates(Coordinates coordinates) { 
	this(coordinates.x(),
	     coordinates.y(),
	     coordinates.z());
	nowhere = coordinates.nowhere;
    }
    public  Coordinates(Coordinates coordinates, double radius) { 
	this(coordinates.x(),
	     coordinates.y(),
	     coordinates.z());
	double xx=radius, yy=radius, zz=radius;
	while (xx*xx + yy*yy + zz*zz > radius*radius) {
	    xx = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    yy = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    zz = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	}
	addToX(xx);
	addToY(yy);
	addToZ(zz);
	random = true;
 	nowhere = false;
   }

    public boolean random() {return random;}
    public final boolean nowhere() {return nowhere;}
    public void somewhere() {nowhere = false;}
    public void reset() {
	set(new Coordinates());
    }
    public final double x() {return x[0];}
    public final double y() {return y[0];}
    public final double z() {return z[0];}

    public final double fx() {return x[1];}
    public final double fy() {return y[1];}
    public final double fz() {return z[1];}

//     public void setX(double x) {this.x[0] = x;}
//     public void setY(double y) {this.y[0] = y;}
//     public void setZ(double z) {this.z[0] = z;}
    public void setXYZ(double x, double y, double z) {
	this.x[0] = x;
	this.y[0] = y;
	this.z[0] = z;
	nowhere = false; 
    }
    public void set(Coordinates other) {
	x[0] = other.x[0];
	y[0] = other.y[0];
	z[0] = other.z[0];
 	nowhere = other.nowhere;
   }
   public void set(Coordinates other, double radius) {
	set(new Coordinates(other, radius));
   }

    public void setFx(double fx) {this.x[1] = fx;}
    public void setFy(double fy) {this.y[1] = fy;}
    public void setFz(double fz) {this.z[1] = fz;}

    public void addToX(double addMe) {x[0] += addMe;}
    public void addToY(double addMe) {y[0] += addMe;}
    public void addToZ(double addMe) {z[0] += addMe;}

    public final void addToFx(double addMe) {x[1] += addMe;}
    public final void addToFy(double addMe) {y[1] += addMe;}
    public final void addToFz(double addMe) {z[1] += addMe;}

    public void resetForces() {
	x[1] = 0.0;
	y[1] = 0.0;
	z[1] = 0.0;
    }
    public final double[] X() {return x;}
    public final double[] Y() {return y;}
    public final double[] Z() {return z;}

    public boolean equals(Object obj) {
	Coordinates other = (Coordinates) obj;
	return
            ((this.z() == other.z())&
	     (this.y() == other.y())&
	     (this.x() == other.x()));
    }
    
    public String toString() {
	return "Coordinates ("+x()+","+fx()+")("+y()+","+fy()+")("+z()+","+fz()+")";
    }
    
    public final double distanceFrom(Coordinates other) {
	double dx = x[0]-other.x[0];
	double dy = y[0]-other.y[0];
	double dz = z[0]-other.z[0];	
	return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }
}
