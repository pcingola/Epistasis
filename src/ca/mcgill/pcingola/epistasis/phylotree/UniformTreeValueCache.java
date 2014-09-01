package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.HashMap;

public class UniformTreeValueCache {

	HashMap<String, Double> cache = new HashMap<String, Double>();

	public void add(String key, double value) {
		cache.put(key, value);
	}

	public String key(PhylogeneticTree tree, int aaCode) {
		return tree.getId() + "\t" + aaCode + "\t" + tree.getUniformCode();
	}

	public Double value(String key) {
		return cache.get(key);
	}

}
