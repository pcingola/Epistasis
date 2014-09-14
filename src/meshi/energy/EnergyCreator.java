package meshi.energy;
import meshi.parameters.*;
import meshi.util.string.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import java.util.*;
import meshi.util.*;
import meshi.util.filters.*;

/**
 * Factory classes that generate energy terms. This is the carpet under which we hide all the boring details of
 * what is needed for the creation of a given energy term. Each energy term requires a specific subclass. 
 * See meshi.energy.bond.BondCreator for a simple example.
 **/ 
public abstract class EnergyCreator implements MeshiPotential, KeyWords {
    /**
     * The weight of the energy term within the total energy function.
     **/
    protected double weight;
    protected AbstractEnergy term = null;
    public AbstractEnergy term() {return term;}

    /**
     * The keyword for the energy term. All commands relevant to this energy terms will sart with this word.
     **/ 
    protected Key key;
    protected boolean weightWasSet = false;
    

    private static String parametersDirectory = null;

    
    /**
     * Constructs an energy creator object. The key parameter will serve as a 
     * keyword that identify relevant commands in the commands list.
     **/ 
    public  EnergyCreator(Key key){
	this.key = key;
    }
 
   /**
     * Constructs an energy creator object. The key parameter will serve as a 
     * keyword that identify relevant commands in the commands list.
     **/ 
    public  EnergyCreator(Key key, CommandList commands){
	this.key = key;
	getWeight(commands);
    }

   //------------------
    /**
     * Construct a somewhat degenerate creator object that cannot read commands from the commands list.
     **/
    public  EnergyCreator(double weight){
	this.weight = weight;
	weightWasSet = true;
    }

   //------------------
    /**
     * The weight of the energy term.
     **/
    public double weight() { 
	if (! weightWasSet)  
	    throw new RuntimeException(key+" was not set");  
	return weight; 
    }	

    //------------------
    /**
     * Extract the path of the parameters directory from the commands list.
     **/
    protected static String parametersDirectory(CommandList commands){ 
       if (parametersDirectory == null)  
	   getParametersDirectory(commands); 
       return parametersDirectory; 
   }
  
   protected static void getParametersDirectory(CommandList commands){ 
       Command command = commands.firstWord(PARAMETERS_DIRECTORY); 
	parametersDirectory = command.secondWord(); 
    } 
 
    //------------------
    /**
     * Where each sub class gets the chance to show what it knows. Offers a standard interface to many energy 
     * functions that have different requirements. 
     **/
    public abstract AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
	 					   CommandList commands);

    /**
     * Get the weight of the energy term.
     **/
    public void getWeight(CommandList commands) {
	if (! weightWasSet) {
	 weight = commands.getWeight(key);
	 weightWasSet = true;
	 System.out.println("weight = "+weight);
	} 
    } 

   public boolean weightWasSet() {return weightWasSet;}

   public String toString() { 
	return "Creator of "+key;  
    }   
}
