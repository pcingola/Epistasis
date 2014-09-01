package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;

public class UniformTreeValueCache {

	HashMap<String, double[]> cache = new HashMap<String, double[]>();
	int size;

	public UniformTreeValueCache(int size) {
		this.size = size;
	}

	/**
	 * Get a cached value
	 */
	public Double get(PhylogeneticTree tree, int aaCode) {
		String key = key(tree);

		double vals[] = cache.get(key);
		if (vals == null) return null;

		double v = vals[aaCode];
		return Double.isNaN(v) ? null : v;
	}

	/**
	 * Create a hash 'key'
	 */
	String key(PhylogeneticTree tree) {
		return tree.getId() + "\t" + tree.getUniformCode();
	}

	/**
	 * Set a value in the cache
	 */
	public synchronized void set(PhylogeneticTree tree, int aaCode, double value) {
		String key = key(tree);
		double vals[] = cache.get(key);

		if (vals == null) {
			vals = new double[size];
			Arrays.fill(vals, Double.NaN);
			cache.put(key, vals);
		}

		vals[aaCode] = value;
	}

}
