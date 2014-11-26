package ca.mcgill.pcingola.epistasis.likelihood;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;

public class TrLikelihoodMatrix {

	boolean debug;
	int len1, len2;
	double llmatrix[][];
	String name;
	String trId1, tr2Id2;
	MultipleSequenceAlignmentSet msas;
	Map<String, Integer> idToIdx1, idToIdx2;
	Map<Integer, String> idxToId1, idxToId2;

	public TrLikelihoodMatrix(MultipleSequenceAlignmentSet msas) {
		this.msas = msas;
	}

	/**
	 * Find best entry
	 */
	public void findBest(int neighbours) {
		double maxScore = Double.NEGATIVE_INFINITY;

		int besti = 0, bestj = 0;
		for (int i = 0; i < len1; i++)
			for (int j = 0; j < len2; j++) {
				double score = score(i, j, neighbours);
				if (score > maxScore) {
					besti = i;
					bestj = j;
					maxScore = score;
				}
			}

		System.out.println(name //
				+ "\t" + neighbours //
				+ "\t" + maxScore //
				+ "\t" + besti //
				+ "\t" + bestj //
				+ "\t" + idxToId1.get(besti) //
				+ "\t" + idxToId2.get(bestj) //
		);

		if (debug) showMatrix(besti, bestj, neighbours);
	}

	/**
	 * Load gene-gene interaction file
	 */
	public boolean load(String geneGeneFile) {
		try {
			//---
			// Read file
			//---
			Timer.showStdErr("Loading gene-gene file:" + geneGeneFile);
			name = geneGeneFile;
			String lines[] = Gpr.readFile(geneGeneFile).split("\n");
			Timer.showStdErr("Done: " + lines.length + " lines.");

			// Find transcript ID
			String ft[] = lines[0].split("\t");
			String f[] = ft[0].split("_");
			String trId1 = f[0] + "_" + f[1];
			f = ft[1].split("_");
			String trId2 = f[0] + "_" + f[1];
			Timer.showStdErr("Transcript IDs: " + trId1 + " \t" + trId2);

			// Find MSAs and lengths and crete indexs
			idToIdx1 = msa2matrixIndex(trId1);
			idxToId1 = reverseIndex(idToIdx1); // Create reverse index

			idToIdx2 = msa2matrixIndex(trId2);
			idxToId2 = reverseIndex(idToIdx2); // Create reverse index

			len1 = msas.getTranscriptLength(trId1);
			len2 = msas.getTranscriptLength(trId2);

			int expLines = len1 * len2; // This doesn't account for blank separating lines
			Timer.showStdErr("Total length: " + len1 + "\t" + len2 + "\tExpected lines: " + expLines);
			if (expLines > lines.length) throw new RuntimeException("Number of lines does not match expected ones!");

			// Parse line and create matrix
			parseMatrix(lines);
		} catch (Throwable e) {
			System.err.println("Error parsing file: " + geneGeneFile);
			e.printStackTrace();
			return false;
		}

		return true;
	}

	/**
	 * Create a map between "msaId[idx]" and matrix index
	 */
	Map<String, Integer> msa2matrixIndex(String trId) {
		List<MultipleSequenceAlignment> msasTr = msas.getMsasByTrIdSorted(trId);
		Map<String, Integer> id2idx = new HashMap<>();

		int idx = 0;
		for (MultipleSequenceAlignment msa : msasTr)
			for (int i = 0; i < msa.getAaSeqLen(); i++) {
				String id = msa.getId() + "[" + i + "]";
				id2idx.put(id, idx++);
			}

		return id2idx;
	}

	/**
	 * Parse lines and create matrix
	 */
	void parseMatrix(String lines[]) {
		//---
		// Create matrix and fill it up
		//---
		llmatrix = new double[len1][len2];

		// Initialize as missing data
		for (int i = 0; i < len1; i++)
			for (int j = 0; j < len2; j++)
				llmatrix[i][j] = Double.NaN;

		// Fill up matrix data
		for (String line : lines) {
			if (line.isEmpty()) continue;

			String ft[] = line.split("\t");
			String id1 = ft[0];
			String id2 = ft[1];
			int idx1 = idToIdx1.get(id1);
			int idx2 = idToIdx2.get(id2);
			llmatrix[idx1][idx2] = Gpr.parseDoubleSafe(ft[2]);
		}

		// Make sure there are no missing values
		int countNa = 0;
		for (int i = 0; i < len1; i++)
			for (int j = 0; j < len2; j++) {
				if (Double.isNaN(llmatrix[i][j])) {
					Gpr.debug("Matrix has NaN entry: [" + i + ", " + j + "]");
					countNa++;
				}
			}

		Timer.showStdErr("Matrix filled. Missing data entries: " + countNa);
		if (countNa > 0) throw new RuntimeException("Missing data. This should never happen!");
	}

	/**
	 * Create a map between "msaId[idx]" and matrix index
	 */
	Map<Integer, String> reverseIndex(Map<String, Integer> idToIdx) {
		Map<Integer, String> idxToId = new HashMap<Integer, String>();

		for (String id : idToIdx.keySet()) {
			if (id == null) throw new RuntimeException("Null ID! This should never happen!");
			int idx = idToIdx.get(id); // Convert to int to check. Missing entry will rise an exception
			idxToId.put(idx, id);
		}

		return idxToId;
	}

	double score(int i, int j, int neighbours) {
		if (neighbours == 0) return llmatrix[i][j];

		double scoreFw = score(i, j, neighbours, 1); // Main diagonal: Forward
		double scoreBack = score(i, j, neighbours, -1); // Reverse diagonal: Backwards

		return scoreFw >= scoreBack ? scoreFw : scoreBack;
	}

	/**
	 * Score average on the forward/reverse diagonal at position [i,j] of llmatrix
	 */
	double score(int i, int j, int neighbours, int stepj) {
		int mini = i - neighbours, maxi = i + neighbours;
		int minj = j - neighbours, maxj = j + neighbours;

		if (mini < 0) mini = 0;

		int startj = j - stepj * neighbours;
		if (startj < 0) startj = 0;

		int count = 0;
		double sum = 0;
		for (int ii = mini; mini <= ii && ii <= maxi && ii < len1; ii++) {
			for (int jj = startj; minj <= jj && jj <= maxj && jj >= 0 && jj < len2; jj += stepj) {
				sum += llmatrix[ii][jj];
				count++;
			}
		}

		return sum / count;
	}

	/**
	 * Show part of the llmatrix (for debugging)
	 */
	void showMatrix(int i, int j, int neigh) {
		// Title
		for (int jj = j - neigh; jj <= j + neigh; jj++)
			System.out.print("\t" + idxToId2.get(jj));
		System.out.println("");

		for (int ii = i - neigh; ii <= i + neigh; ii++) {
			System.out.print(idxToId1.get(ii));
			for (int jj = j - neigh; jj <= j + neigh; jj++)
				System.out.print("\t" + llmatrix[ii][jj]);
			System.out.println("");
		}
	}
}
