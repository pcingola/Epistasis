package meshi.optimizers;
import meshi.util.*;
import meshi.molecularElements.*;
import meshi.energy.*;
import java.util.*;
import meshi.energy.sideChainModelingSolvate.*;
import meshi.energy.simpleEnergyTerms.tether.*;
                     
public class Perturbation  implements KeyWords{
    private TotalEnergy    energy;
    private AbstractEnergy purtubrationTerm;
    private int            maxSteps;
    private Minimizer      minimizer;
    private CommandList    commands;
    private Protein        protein;
    private Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();
    private double         probabilityOfScmod;


    public Perturbation(TotalEnergy energy, CommandList commands, AbstractEnergy purtubrationTerm, Protein protein, double probabilityOfScmod) {
	this.energy           = energy;
	this.maxSteps         = maxSteps;                    
	this.purtubrationTerm = purtubrationTerm;
	this.commands         = commands;
	this.protein          = protein;
	minimizer             = Utils.getLBFGS(energy, commands, MCM_PERTURBATION);
	Minimizer.terminator.reset();
	this.probabilityOfScmod = probabilityOfScmod;
    }
    
    public void perturb() throws OptimizerException {
	SideChainSolvateEnergy sideChainTerm = (SideChainSolvateEnergy)energy.getEnergyTerm(new SideChainSolvateEnergy());
	TetherEnergy           tetherEnergy  = (TetherEnergy)           energy.getEnergyTerm(new TetherEnergy());
	if (sideChainTerm != null) sideChainTerm.off();
	if (tetherEnergy  != null)  tetherEnergy.off();

	Minimizer.terminator.reset();
	purtubrationTerm.on();
	minimizer.run(false);
	purtubrationTerm.off();
    if (randomNumberGenerator.nextDouble() < probabilityOfScmod) {
	    if (sideChainTerm != null) {
		Scmod.scmod(commands, protein, 2,50,energy);
		energy.on();
		sideChainTerm.off();
		tetherEnergy.on();
		Minimizer.terminator.reset();
		minimizer.run(false);
		tetherEnergy.off();
	    }
	}
	Minimizer.terminator.reset();
	energy.on();
	purtubrationTerm.off();
	if (sideChainTerm != null) sideChainTerm.off();
	if (tetherEnergy  != null)  tetherEnergy.off();
    }
}
    
