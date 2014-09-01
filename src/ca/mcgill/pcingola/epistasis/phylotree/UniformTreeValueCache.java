package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;

public class UniformTreeValueCache {

	HashMap<String, double[]> cache = new HashMap<String, double[]>();
	int size;

	public UniformTreeValueCache(int size) {
		this.size = size;
	}

	public synchronized void add(PhylogeneticTree tree, int aaCode, double value) {
		String key = key(tree);
		double vals[] = cache.get(key);

		if (vals == null) {
			vals = new double[size];
			Arrays.fill(vals, Double.NaN);
			cache.put(key, vals);
		}

		vals[aaCode] = value;
	}

	String key(PhylogeneticTree tree) {
		return tree.getId() + "\t" + tree.getUniformCode();
	}

	public Double value(PhylogeneticTree tree, int aaCode) {
		String key = key(tree);
		double vals[] = cache.get(key);
		if (vals == null) return null;
		return vals[aaCode];
	}

}
