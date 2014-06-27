package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Map IDs
 *
 * @author pcingola
 */
public class IdMapper {

	public static String ids(List<IdMapperEntry> ids, Function<IdMapperEntry, String> ime2id) {
		if (ids == null) return null;

		StringBuilder sb = new StringBuilder();

		// Unique names
		HashSet<String> set = new HashSet<String>();
		for (IdMapperEntry ime : ids)
			set.add(ime2id.apply(ime));

		// Sort
		ArrayList<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		// Concatenate
		for (String s : list)
			sb.append((sb.length() <= 0 ? "" : ",") + s);

		return sb.toString();
	}

	int count;
	AutoHashMap<String, ArrayList<IdMapperEntry>> byGeneId, byTrId, byGeneName, byRefSeqId, byPdbId;
	HashSet<IdMapperEntry> entries;

	public IdMapper() {
		this(null);
	}

	public IdMapper(String fileName) {
		ArrayList<IdMapperEntry> emptyList = new ArrayList<IdMapperEntry>();
		byGeneId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byTrId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byGeneName = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byRefSeqId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byPdbId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);;
		entries = new HashSet<>();

		if (fileName != null) load(fileName);
	}

	public void add(IdMapperEntry ime) {
		if (ime.geneId != null) byGeneId.getOrCreate(ime.geneId).add(ime);
		if (ime.trId != null) byTrId.getOrCreate(ime.trId).add(ime);
		if (ime.geneName != null) byGeneName.getOrCreate(ime.geneName).add(ime);
		if (ime.refSeqId != null) byRefSeqId.getOrCreate(ime.refSeqId).add(ime);
		if (ime.pdbId != null) byPdbId.getOrCreate(ime.pdbId).add(ime);
		entries.add(ime);
	}

	public void addAll(Collection<IdMapperEntry> imes) {
		imes.stream().forEach(ime -> add(ime));
	}

	/**
	 * Pick "best" entries
	 */
	public Collection<IdMapperEntry> best(DistanceResults distanceResults) {
		// Sort best entries by geneId
		Map<String, List<IdMapperEntry>> map = entries.stream() //
				.sorted((im1, im2) -> compareToBetterMap(distanceResults, im1, im2)) //
				.collect(Collectors.groupingBy(im -> im.geneId)) //
				;

		// Show
		ArrayList<IdMapperEntry> best = new ArrayList<>();
		map.keySet().forEach(k -> {
			best.add(map.get(k).get(0)); // Add first entry to 'best'
			System.err.println("map\t" + k);
			map.get(k).forEach(im -> System.err.println("\t\t" + im + "\tAA_contacts: " + contacts(distanceResults, im)));
		} //
				);

		return best;
	}

	/**
	 * Comparator: Get "best" mappings first
	 */
	public int compareToBetterMap(DistanceResults distanceResults, IdMapperEntry im1, IdMapperEntry im2) {
		int cmp = im1.geneId.compareTo(im2.geneId);
		if (cmp != 0) return cmp;

		cmp = im2.pdbAaLen - im1.pdbAaLen; // Longer PDB AA sequence first
		if (cmp != 0) return cmp;

		cmp = contacts(distanceResults, im2) - contacts(distanceResults, im1); // Compare number of AA 'in contact' (more is better)
		if (cmp != 0) return cmp;

		cmp = im2.trAaLen - im1.trAaLen; // Longer transcript AA sequence first
		if (cmp != 0) return cmp;

		cmp = im1.trId.compareTo(im2.trId);
		if (cmp != 0) return cmp;

		cmp = im1.pdbId.compareTo(im2.pdbId);
		if (cmp != 0) return cmp;

		cmp = im1.pdbChainId.compareTo(im2.pdbChainId);
		return cmp;
	}

	/**
	 * Number of AA in contact
	 */
	int contacts(DistanceResults distanceResults, IdMapperEntry ime) {
		return distanceResults.contacts(ime.pdbId, ime.pdbChainId);
	}

	public List<IdMapperEntry> getByGeneId(String id) {
		return byGeneId.get(id);
	}

	public List<IdMapperEntry> getByGeneName(String id) {
		return byGeneName.get(id);
	}

	public List<IdMapperEntry> getByPdbId(String id) {
		return byPdbId.get(id);
	}

	public List<IdMapperEntry> getByPdbId(String id, String chainId) {
		if (byPdbId.get(id) == null) return null;
		return byPdbId.get(id) //
				.stream() //
				.filter(ime -> chainId.equals(ime.pdbChainId)) //
				.collect(Collectors.toList());
	}

	public List<IdMapperEntry> getByRefSeqId(String id) {
		return byRefSeqId.get(id);
	}

	public List<IdMapperEntry> getByTrId(String id) {
		return byTrId.get(id);
	}

	public Collection<IdMapperEntry> getEntries() {
		return entries;
	}

	/**
	 * Is there an entry matching all these conditions?
	 */
	public boolean hasEntry(String refSeqId, String pdbId, String pdbChain) {
		List<IdMapperEntry> res = getByRefSeqId(refSeqId);
		if (res == null) return false;

		return res.stream() //
				.filter(ime -> ime.pdbId.equals(pdbId) && ime.pdbChainId.equalsIgnoreCase(pdbChain)) //
				.findAny() //
				.isPresent() //
				;
	}

	void load(String fileName) {
		Timer.showStdErr("Loading IDs from file: " + fileName);
		count = 0;

		Path path = Paths.get(fileName);
		try (Stream<String> lines = Files.lines(path)) {
			lines.filter(l -> !l.startsWith("#")).forEach(l -> parseLine(l));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		Timer.showStdErr("Done. Total entries added: " + count);
	}

	/**
	 * Parse a line and add it to this map
	 */
	void parseLine(String line) {
		IdMapperEntry ime = new IdMapperEntry(line);
		add(ime);
		count++;
	}

	@Override
	public String toString() {
		return entries.stream().sorted().map(im -> im.toString()).collect(Collectors.joining("\n"));
	}
}
