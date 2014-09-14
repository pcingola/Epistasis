package meshi.optimizers;
import meshi.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.AtomList;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.inflate.*;
import meshi.energy.sideChainModelingSolvate.*;
import meshi.energy.simpleEnergyTerms.tether.*;
import meshi.symmetryComplex.proteinsGJ.GJModeller;

import java.util.*;
import java.io.IOException;

public class MCM extends Optimizer {
    private TotalEnergy energy;
    private Minimizer minimizer;
    private TemperatureGenerator temperatureGenerator;
    private Perturbation perturbation;
    private Logger log = null;
    private GJModeller builder = null;         //for SymmetricComplex

    public MCM(TotalEnergy energy, Minimizer minimizer, Perturbation perturbation,
	       TemperatureGenerator temperatureGenerator, int maxSteps) {
	super(energy, maxSteps, 1);
	this.energy = energy;
	this.minimizer = minimizer;
	this.temperatureGenerator = temperatureGenerator;
	this.perturbation = perturbation;
    }

    public MCM(GJModeller builder, TotalEnergy energy, Minimizer minimizer, Perturbation perturbation,
	       TemperatureGenerator temperatureGenerator, int maxSteps) {
	super(energy, maxSteps, 1);
    this.builder = builder;
    this.energy = energy;
	this.minimizer = minimizer;
	this.temperatureGenerator = temperatureGenerator;
	this.perturbation = perturbation;
    }


    public OptimizerStatus run(Logger log) throws OptimizerException{
	this.log = log;
	return run();
    }

    public OptimizerStatus run() throws OptimizerException{
	double oldEnergy = 1000000;
	double[][] oldCoordinates;
	double currentEnergy;
	double dE;
	Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

	SideChainSolvateEnergy sideChainTerm = (SideChainSolvateEnergy) energy.getEnergyTerm(new SideChainSolvateEnergy());
	TetherEnergy           tetherEnergy  = (TetherEnergy)           energy.getEnergyTerm(new TetherEnergy());
	if (sideChainTerm != null) sideChainTerm.off();
	if (tetherEnergy  != null)  tetherEnergy.off();

	oldEnergy = energy.evaluate();
	for (int step = 1; step <= maxSteps; step++) {
	    double temperature = temperatureGenerator.next();
	    System.out.println("\n step # "+step+"; temperature "+temperature+" \ncurrent energy = \n"+energy.report(step));
	    
	    oldCoordinates = getOldCoordinates(energy);	
	    try {
		if (step > 1) {
            AtomList atomListCopy = Utils.duplicateInAnewMolecularSystem( energy.atomList());
            int counter = 0;
            boolean success = false;
            while ((counter < 5) & (!success)) {
                try {
                    perturbation.perturb();
                    success = true;
                } catch (RuntimeException ex) {
                    System.out.println("A problem in perturbation #  "+counter);
                     System.out.println(ex);
                    System.out.println("RMS reached\t"+counter+"\t"+step+"\t"+energy.atomList().getRms(atomListCopy ));
                    setCoordinates(energy,oldCoordinates);
                    counter++;
                }
            }
            if (!success) throw new RuntimeException("SHAVROO ET HAKELIM");
            System.out.println("RMS perturbation\t"+step+"\t"+energy.atomList().getRms(atomListCopy ));
        }

        OptimizerStatus os = minimizer.run();
		System.out.println("MCM step "+os);
		if (log != null) log.mcm(energy,step);
		currentEnergy = energy.evaluate();
		dE = currentEnergy - oldEnergy;
		if (dE > 0) {
		    double rnd = randomNumberGenerator.nextDouble();
		    if (rnd > Math.exp(-dE/temperature)) {//That is Metropolis criterion failed
			setCoordinates(energy,oldCoordinates);
			System.out.println("This step failed: dE = "+dE+
					   "  dE/temperature = "+ dE/temperature+
					   "  Math.exp(-dE/temperature) = "+Math.exp(-dE/temperature)+
					   "  rnd = "+rnd);
			if (log != null) log.mcm(energy,step);
		    }
		    else  {
			System.out.println("This step succeeded: dE = "+dE+
				       "  dE/temperature = "+ dE/temperature+
					   "  Math.exp(-dE/temperature) = "+Math.exp(-dE/temperature)+
					   "  rnd = "+rnd); 
			oldEnergy = currentEnergy;
		    }
		}
		else {
		    System.out.println("This step succeeded: dE = "+dE);
		    oldEnergy = currentEnergy;
		}
	    }
	    catch (OptimizerException ox) {
		System.out.println("This step failed due to OptimizerException");
		setCoordinates(energy,oldCoordinates);
	    }
        //for SymmetricComplex
        try{
            if (builder != null) builder.generateFiles(builder.gj, builder.creator());            
        }
        catch(IOException e) {
            System.out.println("This step failed to make pdb of SymmetricComplex");
            e.printStackTrace();
		}
    }
	return OptimizerStatus.DONE;
    }
    
    private static double[][] getOldCoordinates(TotalEnergy energy) {
	double[][] coordinates = energy.coordinates();
	int length = coordinates.length;
	double[][] out = new double[length][2];

	for (int i = 0; i < length; i++) {
	    out[i][0] = coordinates[i][0];
	    out[i][1] = coordinates[i][1];
	}
	
	return out;
    }

    private static void setCoordinates(TotalEnergy energy, double[][] toSet) {
	double[][] coordinates = energy.coordinates();
	int length = coordinates.length;
	if (length != toSet.length) throw new RuntimeException("Weird parameters to MCM.setCoordinates");

	for (int i = 0; i < length; i++) {
	    coordinates[i][0] = toSet[i][0];
	    coordinates[i][1] = toSet[i][1];
	}
    }
}    
	    
	    
	
	
	
    
/*
	    Inflate inflate = (Inflate) energy.getEnergyTerm(new Inflate());
	    if (inflate == null) throw new RuntimeException("No point in MCM without inflate");
	    if (iteration != 1) {
		inflate.on();
		System.out.println(optimizer.run());
	    }
	    inflate.off();
	    System.out.println(optimizer.run());
*/
