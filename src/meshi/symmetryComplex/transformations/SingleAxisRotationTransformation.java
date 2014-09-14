package meshi.symmetryComplex.transformations;

import meshi.symmetryComplex.transformations.Transformation;

/**
 * Source: <a href="http://en.wikipedia.org/wiki/Rotation_matrix">
 * Rotation matrix</a> (Wikipedia).
 * 
 * @author Oren Wolfshtat
 */
public class SingleAxisRotationTransformation extends Transformation {

    public SingleAxisRotationTransformation(String axis, double degrees) {
        super();

        double radians = Math.toRadians(degrees);

        if ( axis.equalsIgnoreCase("x") )
            setXRotation(radians);

        else if ( axis.equalsIgnoreCase("y") )
            setYRotation(radians);

        else if ( axis.equalsIgnoreCase("z") )
            setZRotation(radians);

        else throw new RuntimeException
            ("Illegal axis.");
    }

    private void setXRotation(double radians) {

        this.setMember(1, 1, Math.cos(radians));
        this.setMember(1, 2, -Math.sin(radians));
        this.setMember(2, 1, Math.sin(radians));
        this.setMember(2, 2, Math.cos(radians));
    }

    private void setYRotation(double radians) {

        this.setMember(0, 0, Math.cos(radians));
        this.setMember(0, 2, Math.sin(radians));
        this.setMember(2, 0, -Math.sin(radians));
        this.setMember(2, 2, Math.cos(radians));
    }

    private void setZRotation(double radians) {

        this.setMember(0, 0, Math.cos(radians));
        this.setMember(0, 1, -Math.sin(radians));
        this.setMember(1, 0, Math.sin(radians));
        this.setMember(1, 1, Math.cos(radians));
    }
}