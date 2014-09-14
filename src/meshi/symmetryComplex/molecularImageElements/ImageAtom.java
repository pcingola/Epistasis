package meshi.symmetryComplex.molecularImageElements;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomStatus;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.geometry.Coordinates;


/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 22/12/2008
 * Time: 12:13:21
 * To change this template use File | Settings | File Templates.
 */
public class ImageAtom extends Atom {
    private Atom source;
    private Transformation transformation;

    public ImageAtom(Atom source, Transformation transformation) {
       super(source.name, source.residue(), source.type(),
               new Coordinates(source), source.temperatureFactor());
       this.source = source;
       this.transformation = transformation;
       buildLocation();
    }

    public void buildLocation() {
        if (source.nowhere())
            //setXYZ(source.x(),source.y(),source.z(),AtomStatus.IMAGE);
              setXYZ(-999.999,-999.999,-999.999,AtomStatus.IMAGE);
        else {
         double[][] xyz = new double[4][1];

         xyz[0][0] = source.x();
         xyz[1][0] = source.y();
         xyz[2][0] = source.z();
         xyz[3][0] = 1.0;

         xyz = transformation.transform(xyz);
         setXYZ(xyz[0][0], xyz[1][0], xyz[2][0], AtomStatus.IMAGE);
        }

   //     if (core.status() != AtomStatus.IMAGE)
     //       throw new MeshiException
       //     ("This method should update image atoms only.");
     }

   protected void updateLocation() {
       if (source.nowhere())
                 throw new RuntimeException("This image atom is still nowhere!"+this);
       else {
        double[][] xyz = new double[4][1];

        xyz[0][0] = source.x();
        xyz[1][0] = source.y();
        xyz[2][0] = source.z();
        xyz[3][0] = 1.0;

        xyz = transformation.transform(xyz);
        setXYZ(xyz[0][0], xyz[1][0], xyz[2][0], AtomStatus.IMAGE);
       }

  //     if (core.status() != AtomStatus.IMAGE)
    //       throw new MeshiException
      //     ("This method should update image atoms only.");
    }

  public Atom source(){return source;}
  //public boolean image() {return (core.status() == AtomStatus.IMAGE);}
  public final Transformation transformation() {return transformation;}


}
