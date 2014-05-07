package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Map IDs
 * 
 * @author pcingola
 */
public class IdMapper {

	int count;
	AutoHashMap<String, ArrayList<IdMapperEntry>> byGeneId, byTrId, byProtId, byGeneName, byPfamId, byPdbId;

	public static String trIds(List<IdMapperEntry> ids) {
		if (ids == null) return null;

		StringBuilder sb = new StringBuilder();

		// Unique names
		HashSet<String> set = new HashSet<String>();
		for (IdMapperEntry id : ids)
			set.add(id.trId);

		// Sort
		ArrayList<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		// Concatenate
		for (String s : list)
			sb.append((sb.length() <= 0 ? "" : ",") + s);

		return sb.toString();
	}

	public IdMapper(String fileName) {
		ArrayList<IdMapperEntry> emptyList = new ArrayList<IdMapperEntry>();
		byGeneId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byTrId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byProtId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byGeneName = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byPfamId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);
		byPdbId = new AutoHashMap<String, ArrayList<IdMapperEntry>>(emptyList);;

		load(fileName);
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

	public List<IdMapperEntry> getByPfamId(String id) {
		return byPfamId.get(id);
	}

	public List<IdMapperEntry> getByProtId(String id) {
		return byProtId.get(id);
	}

	public List<IdMapperEntry> getByTrId(String id) {
		return byTrId.get(id);
	}

	void load(String fileName) {
		Timer.showStdErr("Loading IDs from file: " + fileName);
		count = 0;

		Path path = Paths.get(fileName);
		try (Stream<String> lines = Files.lines(path)) {
			lines//
			.filter(l -> !l.startsWith("#")) //
					.forEach(l -> parseLine(l)) //
			;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		Timer.showStdErr("Done. Total entries added: " + count);
	}

	/**
	 * Parse a line and add it to this map
	 * @param line
	 */
	void parseLine(String line) {
		IdMapperEntry ime = new IdMapperEntry(line);

		if (ime.geneId != null) byGeneId.getOrCreate(ime.geneId).add(ime);
		if (ime.trId != null) byTrId.getOrCreate(ime.trId).add(ime);
		if (ime.protId != null) byProtId.getOrCreate(ime.protId).add(ime);
		if (ime.geneName != null) byGeneName.getOrCreate(ime.geneName).add(ime);
		if (ime.pfamId != null) byPfamId.getOrCreate(ime.pfamId).add(ime);
		if (ime.pdbId != null) byPdbId.getOrCreate(ime.pdbId).add(ime);

		count++;
	}
}
