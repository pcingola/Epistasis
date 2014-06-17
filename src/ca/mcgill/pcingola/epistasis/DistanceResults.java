package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A collection of DistanceResult
 */
public class DistanceResults extends ArrayList<DistanceResult> {

	private static final long serialVersionUID = 1L;

	HashMap<String, DistanceResult> byKey;

	DistanceResults() {
		super();
	}

	/**
	 * Add an element only if the 'key' is unique
	 */
	public void addIfUniq(DistanceResult d, String key) {
		if (byKey == null) byKey = new HashMap<>();
		if (!byKey.containsKey(key)) {
			add(d);
			byKey.put(key, d);
		}
	}

	/**
	 * Add elements, keep 'min' by key
	 */
	public void addMin(DistanceResult d, String key) {
		if (byKey == null) byKey = new HashMap<>();
		if (!byKey.containsKey(key)) {
			add(d);
			byKey.put(key, d);
		} else {
			DistanceResult dold = byKey.get(key);
			if (dold.equalPos(d) && d.compareByPos(dold) < 0) // Same position? Keep smallest one
				byKey.put(key, d);
		}
	}

	/**
	 * Load from file
	 */
	public void load(String fileName) {
		for (String line : Gpr.readFile(fileName).split("\n"))
			add(new DistanceResult(line));
	}

	@Override
	public String toString() {
		return stream().map(d -> d.toString()).collect(Collectors.joining("\n"));
	}

}
