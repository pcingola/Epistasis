package meshi.energy;

import java.util.ArrayList;
import java.util.HashMap;

import meshi.molecularElements.atoms.Atom;
import meshi.util.Attributable;
import meshi.util.MeshiAttribute;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;

/**
 * A super class for all meshi energy terms.
 * Currently oriented towards drivable energy functions. See meshi.energy.bond.BondEnergy.
 **/
public abstract class AbstractEnergy implements Updateable, Attributable {
	private static class IsEnergy implements Filter {
		@Override
		public boolean accept(Object obj) {
			return (obj instanceof AbstractEnergy);
		}
	}

	/**
	 * A list of updateable elements (implementing the meshi.util.Updateable interface). 
	 **/
	protected class UpdateableList extends ArrayList<Updateable> {
		public UpdateableList() {
			super();
		}

		public UpdateableList(Updateable[] array) {
			super();
			for (Updateable updatable : array)
				add(updatable);
		}

		/**
		 * Iterates over the list updates each of its elements.
		 **/
		public void update(int numberOfUpdates) throws UpdateableException {
			for (Updateable resource : this) {
				resource.update(numberOfUpdates);
			}
		}
	}

	/**
	 * Energy function/term name.
	 **/
	/**
	 * A short name for the energy term.
	 **/
	protected String comment = "Abstract energy";

	/** 
	 * The weight of the energy term within the total energy.
	 **/
	protected double weight;

	/**
	 * The energy term is evaluated only if it is on.
	 **/
	protected boolean on = true;
	public static final Filter filter = new IsEnergy();

	/**
	 * 1/0 good for dubugging.
	 **/
	public static final double INFINITY = 1 / 0.0;

	/**
	 * sqrt(-1) good for dubugging.
	 **/
	public static final double NaN = Math.sqrt(-1.0);

	/*
	 * A list of all the updateable resources needed by the energy term. For a detailed description 
	 * of MESHI's treatment of updatable resources see the documentation of TotalEnergy.
	 **/
	protected UpdateableList updateableResources;

	//----------------------------------- Attributes ---------------------------------------     
	private HashMap attributes = new HashMap();

	//--------------------------------------- auxiliary methods--------------------------------------------------
	/**
	 * Generates an empty array.
	 **/
	protected static Updateable[] toArray() {
		Updateable[] out = {};
		return out;
	}

	/**
	* Generates an aray with one element - the parameter.
	**/
	protected static Updateable[] toArray(Updateable o1) {
		Updateable[] out = { o1 };
		return out;
	}

	protected static Updateable[] toArray(Updateable o1, Updateable o2) {
		Updateable[] out = { o1, o2 };
		return out;
	}

	protected static Updateable[] toArray(Updateable o1, Updateable o2, Updateable o3) {
		Updateable[] out = { o1, o2, o3 };
		return out;
	}

	// ------------------------------------- internal auxiliary classes -----------------------------------------

	/**
	 * Construct a dummy object that cannot be evaluated. Such an object is useful as a key 
	 * during searches.
	 **/
	public AbstractEnergy() {
	}

	/** 
	 * construct An energy term. 
	 **/
	public AbstractEnergy(Updateable[] updateableResources, double weight) {
		this.updateableResources = new UpdateableList(updateableResources);
		this.weight = weight;
	}

	@Override
	public final void addAttribute(MeshiAttribute attribute) {
		attributes.put(attribute.key(), attribute);
	}

	public String comment() {
		return comment;
	}

	/**
	    * Evaluates the energy term and <b> update </b> the derivatives.
	    **/
	public abstract double evaluate();

	/**
	    * Evaluates the energy term and devides the energy between the atoms. The energy field of
	    * each atom is assigned a value - its contribution to the total energy sum.
	    **/
	public abstract void evaluateAtoms();

	@Override
	public final MeshiAttribute getAttribute(int key) {
		return (MeshiAttribute) attributes.get(key);
	}

	//----------------------------------------- housekeeping ---------------------------------------------
	public void handleMissingParameters(Object obj) {
		throw new RuntimeException("Missing parameters for:\n" + obj);
	}

	public boolean isOn() {
		return on;
	}

	/**
	 * Turnes the energyTerm OFF.
	 **/
	public void off() {
		on = false;
	}

	/**
	 * Turnes the energyTerm ON.
	 **/
	public void on() {
		on = true;
	}

	/**
	 * Looking for one "criminal" atom whose derivation is wrong.
	 */
	public abstract void test(TotalEnergy totalEnergy, Atom atom);

	@Override
	public String toString() {
		return comment;
	}

	/**
	 * Updates the updatable resources. For a detailed description 
	 * of MESHI's treatment of updatable resources see the documentation of TotalEnergy.
	 **/
	@Override
	public void update(int numberOfUpdates) throws UpdateableException {
		updateableResources.update(numberOfUpdates);
	}
}
