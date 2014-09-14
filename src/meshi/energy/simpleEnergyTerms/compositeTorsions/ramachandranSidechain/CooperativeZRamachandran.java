package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.ResidueTorsionsPropensityAttribute;
import meshi.parameters.ResidueType;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiAttribute;
import meshi.util.UpdateableException;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 13:21:29
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZRamachandran extends AbstractEnergy implements CompositeTorsionsDefinitions {
    private static final double[] mean       = new double[ResidueType.values().length];
    private static final double[] std         = new double[ResidueType.values().length];
    
        private RamachandranSidechainEnergy ramachandranSidechainEnergy;
        private ResidueTorsionsList residueTorsionsList;

        public CooperativeZRamachandran(RamachandranSidechainEnergy ramachandranSidechainEnergy, double weight, CooperativeZRamachandranParameters parameters) {
        super(toArray(), weight);
        this.ramachandranSidechainEnergy = ramachandranSidechainEnergy;
        residueTorsionsList = ramachandranSidechainEnergy.residueTorsionsList;
        comment = "CooperativeZRamachandran";
        for (ResidueType type:ResidueType.values()) {
            mean[type.ordinal()] = parameters.mean[type.ordinal()];
            std[type.ordinal()] = parameters.std[type.ordinal()];
        }

        }

       public double evaluate() {
        return evaluate(false);
       }
        public void evaluateAtoms() {
        evaluate(true);
        }

    public double evaluate(boolean evaluateAtoms) {
	if (! on) return 0.0;
        double energy, e = 0, factor;
        Residue residue;
        ResidueType type;

         for (ResidueTorsions rt :residueTorsionsList){
                    type = rt.getResidueType();
                     if (std[type.ordinal()] != 0)
                             e += (rt.energy() -  mean[type.ordinal()])/std[type.ordinal()];
         }

         energy = weight*e*e;

         if (evaluateAtoms) {
             int nAtoms = 0;
             for (ResidueTorsions residueTorsions:residueTorsionsList) {
                                 nAtoms += residueTorsions.getResidue().atoms().size();
             }
             for (ResidueTorsions residueTorsions:residueTorsionsList) {
                 residue = residueTorsions.getResidue();
                 AtomList atoms = residue.atoms();
                 for (Atom atom:atoms)
                     atom.addEnergy((weight * e * e)/nAtoms);
                 }
             }

        ResidueTorsionsAttribute rta;
        for (ResidueTorsions residueTorsions:residueTorsionsList) {
           rta = (ResidueTorsionsAttribute) residueTorsions.getAttribute(MeshiAttribute.RESIDUE_TORSIONS_ATTRIBUTE);
           if (rta == null) continue; //if all atoms of residue are frozen
           residue = residueTorsions.getResidue();
           type   = residue.type();
            if (std[type.ordinal()] == 0 )
                            throw new RuntimeException("Something is weird in Parameters or residueTorsionsList of"+this);
           factor = 2*e/std[type.ordinal()]*weight;
            residueTorsions.applyForce( PHI, -rta.phi_deriv*factor );
            residueTorsions.applyForce( PSI, -rta.psi_deriv*factor);
            if (rta.chi_1_deriv != 0) residueTorsions.applyForce( CHI_1, -rta.chi_1_deriv*factor );
            if (rta.chi_2_deriv != 0) residueTorsions.applyForce( CHI_2, -rta.chi_2_deriv*factor );
            if (rta.chi_3_deriv != 0) residueTorsions.applyForce( CHI_3, -rta.chi_3_deriv*factor );
            if (rta.chi_4_deriv != 0) residueTorsions.applyForce( CHI_4, -rta.chi_4_deriv*factor );

       }
       return energy;
       }


     public void test(TotalEnergy totalEnergy, Atom atom) {
         if (! on) {System.out.println(""+this +" is off"); return;}
     System.out.println("Testing "+this+" using "+atom);
         if (atom == null)
         throw new RuntimeException("Cannot test "+this);

         double[][] coordinates = new double[3][];
         coordinates[0] = atom.X();
         coordinates[1] = atom.Y();
         coordinates[2] = atom.Z();
         for(int i = 0; i< 3; i++) {
             try{totalEnergy.update();}catch(UpdateableException ue){}
             ramachandranSidechainEnergy.evaluate();
             double x = coordinates[i][0];
             coordinates[i][1] = 0;
             double e1 = evaluate();
             double analiticalForce = coordinates[i][1];
             coordinates[i][0] += EnergyElement.DX;
             // Whatever should be updated ( such as distance matrix torsion list etc. )
             try{totalEnergy.update();}catch(UpdateableException ue){}
             ramachandranSidechainEnergy.evaluate();
             double e2 = evaluate();
             double de = e2-e1;
             double numericalForce = - de/ EnergyElement.DX;
             coordinates[i][0] -=  EnergyElement.DX;
             try{totalEnergy.update();}catch(UpdateableException ue){}

             double diff = Math.abs(analiticalForce - numericalForce);

             if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+EnergyElement.VERY_SMALL)) > EnergyElement.relativeDiffTolerance){
                 System.out.println("Testing "+this);
                 System.out.println("Atom["+atom.number()+"]."+EnergyElement.XYZ.charAt(i)+" = "+x);
                 System.out.println("Analytical force = "+analiticalForce);
                 System.out.println("Numerical force  = "+numericalForce);

                 System.out.println("diff = "+diff+"\n"+
                                    "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+EnergyElement.VERY_SMALL) = "+
                                    2*diff/(Math.abs(analiticalForce) + Math.abs(numericalForce)+EnergyElement.VERY_SMALL));
                 System.out.println();
             }
             if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                 System.out.println("Testing "+this+"\ne1 = "+e1);
             if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                 System.out.println("Testing "+this+"\ne2 = "+e2);
             if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                 System.out.println("Testing "+this+"\nanaliticalForce = "+analiticalForce);
         }


    }
}
