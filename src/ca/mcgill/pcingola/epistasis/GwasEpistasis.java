package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.util.Tuple;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Perform GWAS using epistasis data
 *
 * @author pcingola
 */
public class GwasEpistasis {

	public static int SHOW_EVERY_VCF = 1000;
	public static int SHOW_EVERY_GENES_LL = 10000;
	public static int SHOW_LINE_GENES_LL_EVERY = 100 * SHOW_EVERY_GENES_LL;
	public static double SHOW_LINE_LL_MIN = 1.0;
	public static int MINOR_ALLELE_COUNT = 5;
	public static double LL_THRESHOLD = 6.0;

	// Splits information
	boolean analyzeAllPairs = false; // Use for testing and debugging
	boolean debug = false;
	boolean verbose = false;
	int splitI = 2;
	int splitJ = 3;
	int numSplits = 100;
	int countOk, countErr;
	double llThreshold = LL_THRESHOLD;
	String configFile;
	String logLikelihoodFile; // Log likelihood file (epistatic model)
	String vcfFile;
	String genomeVer, pdbDir;
	String phenoCovariatesFile;
	List<MarkerPairLikelihood> llpairs; // Gene log-likelihood entries
	List<byte[]> gtsSplitI, gtsSplitJ;
	List<String> gtIdsSplitI, gtIdsSplitJ; // VCF IDs (corresponding to genotypes)
	Map<String, Transcript> trancriptById; // Transcript by (incomplete) transcript ID (no version number is used)
	Map<String, Marker> llmarkerById = new HashMap<String, Marker>(); // log-likelihood markers by ID
	Map<Long, LikelihoodAnalysis2> llAnByThreadId = new HashMap<>();
	AutoHashMap<String, ArrayList<byte[]>> gtById; // Genotypes by ID
	IntervalForest intForest; // Interval forest (MSAs intervals)
	PdbGenomeMsas pdbGenomeMsas;
	InteractionLikelihood interactionLikelihood;

	public GwasEpistasis(PdbGenomeMsas pdbGenomeMsas, InteractionLikelihood interactionLikelihood, String vcfFile, String phenoCovariatesFile, int numSplits, int splitI, int splitJ) {
		this.pdbGenomeMsas = pdbGenomeMsas;
		this.interactionLikelihood = interactionLikelihood;
		this.vcfFile = vcfFile;
		this.phenoCovariatesFile = phenoCovariatesFile;
		this.numSplits = numSplits;
		this.splitI = splitI;
		this.splitJ = splitJ;
	}

	public GwasEpistasis(String configFile, String genomeVer, String logLikelihoodFile, String vcfFile, String phenoCovariatesFile) {
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		this.logLikelihoodFile = logLikelihoodFile;
		this.vcfFile = vcfFile;
		this.phenoCovariatesFile = phenoCovariatesFile;

		// No split
		numSplits = 1;
		splitI = 1;
		splitJ = 1;
	}

	/**
	 * Build interval forest using LogLik markers
	 */
	protected void buildForest() {
		if (intForest != null) return; // Nothing to do

		Timer.showStdErr("Building interval forest: MSAs");

		// Create a collections of 'markers'
		Markers markers = new Markers();
		Genome genome = pdbGenomeMsas.getConfig().getSnpEffectPredictor().getGenome();
		for (MultipleSequenceAlignment msa : pdbGenomeMsas.getMsas()) {
			// Create a marker for this MSA and add it to 'markers'
			Marker m = new Marker(genome.getChromosome(msa.getChromo()), msa.getStart(), msa.getEnd(), false, msa.getId());
			markers.add(m);
		}

		// Create forest and build
		intForest = new IntervalForest(markers);
		intForest.build();
		Timer.showStdErr("Done. Added " + markers.size() + " markers.");
	}

	public int getCountErr() {
		return countErr;
	}

	public int getCountOk() {
		return countOk;
	}

	/**
	 * Get a likelihood analysis object for each thread
	 */
	LikelihoodAnalysis2 getLikelihoodAnalysis2() {
		long threadId = Thread.currentThread().getId();

		LikelihoodAnalysis2 llan = llAnByThreadId.get(threadId);

		if (llan == null) {
			llan = new LikelihoodAnalysis2(phenoCovariatesFile, vcfFile);
			llan.init();
			llAnByThreadId.put(threadId, llan);
		}
		return llan;
	}

	/**
	 * Perform GWAS analysis
	 */
	public void gwas() {
		initialize(); // Initialize
		readVcf(); // Read VCF file

		//---
		// Test VCF entries
		//---
		Counter count = new Counter();
		Counter countLl = new Counter();
		for (int idxi = 0; idxi < gtIdsSplitI.size(); idxi++) {

			// Split_i info
			final int i = idxi;
			String idi = gtIdsSplitI.get(idxi);
			byte gti[] = gtsSplitI.get(idxi);

			int minJ = 0;
			if (splitI == splitJ) minJ = idxi + 1;

			// Parallel on split_j
			IntStream.range(minJ, gtsSplitJ.size()) //
					.parallel() //
					.forEach(j -> {
						double ll[] = gwas(idi, gti, gtIdsSplitJ.get(j), gtsSplitJ.get(j));
						double llTot = ll[0] + ll[1];
						if (llTot > llThreshold) countLl.inc();
						if (llTot != 0.0) Timer.show(count.inc() + " (" + i + " / " + j + ")\t" + countLl + "\tll: " + llTot + "\tll_LogReg: " + ll[0] + "\tll_MSA: " + ll[1] + "\t" + idi + "\t" + gtIdsSplitJ.get(j));
					});
		}
	}

	/**
	 * Perform analysis on genotypes 'i' and 'j'
	 * @return [ll_LogRes, ll_Msa]
	 */
	double[] gwas(String idi, byte gti[], String idj, byte gtj[]) {
		//---
		// Likelihood based on logistic regression
		//---
		LikelihoodAnalysis2 llan = getLikelihoodAnalysis2();
		double llLogReg = llan.logLikelihood(idi, gti, idj, gtj);

		// Store results
		double res[] = new double[2];
		res[0] = llLogReg;

		// Log likelihood form logistic regression is too low?
		// => Don't bother to calculate next part
		if (llLogReg < llThreshold) return res;

		//---
		// Likelihood based on interaction
		//---

		// Find corresponding MSA ID and index for both genotypes
		Tuple<String, Integer> msaIdxI = id2MsaAa(idi);
		if (msaIdxI == null) return res;

		Tuple<String, Integer> msaIdxJ = id2MsaAa(idj);
		if (msaIdxJ == null) return res;

		// Likelihood based on epistatic interaction
		String msaId1 = msaIdxI.first, msaId2 = msaIdxJ.first;
		int msaIdx1 = msaIdxI.second, msaIdx2 = msaIdxJ.second;
		double llMsa = interactionLikelihood.logLikelihoodRatio(msaId1, msaIdx1, msaId2, msaIdx2, false);
		res[1] = llMsa;

		// Combine
		return res;
	}

	/**
	 * Create a marker from the VCF ID
	 */
	protected Marker id2marker(String vcfId) {
		String idRep = vcfId.replace(':', '_').replace('-', '_');
		String f[] = idRep.split("_");
		String chr = f[0];
		int start = Gpr.parseIntSafe(f[1]) - 1;
		String ref = f[0];
		int end = start + ref.length() - 1;

		// Find chromo and create marker
		Chromosome chromo = pdbGenomeMsas.getConfig().getGenome().getChromosome(chr);
		return new Marker(chromo, start, end, false, vcfId);
	}

	/**
	 * Find MSAid and AaIdx for a genomic position (given as an ID string)
	 */
	Tuple<String, Integer> id2MsaAa(String idj) {
		String resMsaId = null;
		int resAaIdx = -1;

		// Create a marker, find all MSAs that intercept the marker
		Marker mj = id2marker(idj);
		Markers resj = intForest.query(mj);

		// We now need to find the AA index for that MSA
		String seqPrev = null;
		for (Marker m : resj) {
			String msaId = m.getId();
			int aaIdx = pdbGenomeMsas.genomicPos2AaIdx(msaId, mj.getStart());

			if (aaIdx >= 0) {
				// Found something: Store result
				// Note: We could just stop here, all the rest is done just to
				//       ensure that the mapping works OK and that we are not
				//       finding inconsistent sequences
				if (resMsaId == null) {
					resMsaId = msaId;
					resAaIdx = aaIdx;
				}

				// Check sequence length
				MultipleSequenceAlignment msa = pdbGenomeMsas.getMsas().getMsa(msaId);
				if (aaIdx >= msa.getSeqLen()) {
					Gpr.debug("ERROR: Index out of range !"//
							+ "\n\tID_J              : " + idj //
							+ "\n\tMarker            : " + m.toStr() //
							+ "\n\tmsa.Id            : " + msaId //
							+ "\n\tmsa.aaIdx         : " + aaIdx //
					);
					return null;
				}

				// Check sequence
				String colSeq = msa.getColumnString(aaIdx);
				if (seqPrev != null && !colSeq.equals(seqPrev)) Gpr.debug("ERROR: Column sequences differ!"//
						+ "\n\tID_J              : " + idj //
						+ "\n\tMarker            : " + m.toStr() //
						+ "\n\tmsa.Id            : " + msaId //
						+ "\n\tmsa.aaIdx         : " + aaIdx //
						+ "\n\tColumn Seq        : " + colSeq //
						+ "\n\tColumn Seq (prev) : " + seqPrev//
				);

				// Store for next iteration
				seqPrev = colSeq;
			}
		}

		// Resturn result, if any
		if (resMsaId == null) return null;
		return new Tuple<String, Integer>(resMsaId, resAaIdx);
	}

	/**
	 * Initialize
	 */
	public void initialize() {
		if (pdbGenomeMsas == null) {
			pdbGenomeMsas = new PdbGenomeMsas(configFile, genomeVer, pdbDir, null);
			pdbGenomeMsas.initialize();
		}

		// Pre-calculate matrix exponentials
		if (interactionLikelihood != null) interactionLikelihood.precalcExps();
	}

	/**
	 * Convert to minor allele (or filter out)
	 * @return A minor allele genotype, or null if it doesn't satisfy some fitlering requirements)
	 */
	byte[] minorAllele(byte gt[]) {
		// Count alleles
		int ac = 0;
		for (int i = 0; i < gt.length; i++)
			if (gt[i] > 0) ac += gt[i]; // Don't count '-1' (i.e. missing genotypes)

		if (ac < MINOR_ALLELE_COUNT) return null; // Too few alleels: Filter out
		if (ac <= gt.length) return gt; // OK, gt[] is mainor allele

		// Convert to minor allele
		for (int i = 0; i < gt.length; i++)
			if (gt[i] >= 0) gt[i] = (byte) (2 - gt[i]);

		return gt;

	}

	/**
	 * Parse MSA_ID from gene likelihood file
	 *
	 * Format: trId _ chr : start - end [ aaIdx ]
	 *
	 * E.g.  : 	NM_004635_3:50679137-50679201[0]
	 * 			NM_006936_21:46226865-46226954[22]
	 *
	 * @return A marker that contains the codon referenced by the ID
	 *			Note: The marker has ALL bases in the codon.
	 *				  For instance, is the codon is split between two exons, the
	 *				  marker will contain the intron
	 *
	 */
	protected Marker parseMsaId(String id, char aaExpected) {
		// Try to find cached copy
		Marker marker = llmarkerById.get(id);
		if (marker != null) {
			countOk++;
			showCount(true);
			return marker;
		}

		String idRep = id.replace(':', '_').replace('-', '_').replace('[', '_').replace(']', '_');
		String f[] = idRep.split("_");
		String trId = f[0] + "_" + f[1];
		String chr = f[2];
		int start = Gpr.parseIntSafe(f[3]);
		int end = Gpr.parseIntSafe(f[4]);
		int aaIdx = Gpr.parseIntSafe(f[5]);

		// Create a marker that lies onto the referred AA
		marker = pdbGenomeMsas.markerMsa(trId, chr, start, end, aaIdx, aaExpected);

		// Some accounting
		if (marker == null) {
			countErr++;
			showCount(false);
		} else {
			countOk++;
			showCount(true);
			llmarkerById.put(id, marker); // Cache marker
		}

		return null;
	}

	/**
	 * Perform GWAS analysis using epistatic data
	 */
	public void readGenesLogLikelihood() {
		llpairs = new ArrayList<MarkerPairLikelihood>();

		//---
		// Read "genes likelihood" file
		//---
		Timer.showStdErr("Reading genes likelihood file '" + logLikelihoodFile + "'.");
		int count = 0;
		LineFileIterator lfi = new LineFileIterator(logLikelihoodFile);
		for (String line : lfi) {
			if (line.isEmpty()) continue;

			// Parse line
			String f[] = line.split("\t");
			String msaId1 = f[0];
			String msaId2 = f[1];
			double logLikRatio = Gpr.parseDoubleSafe(f[2]);
			String seq1 = f[5];
			String seq2 = f[6];

			// Filter by log likelihood
			count++;

			// Create MarkerPair
			Marker m1 = parseMsaId(msaId1, seq1.charAt(0));
			Marker m2 = parseMsaId(msaId2, seq2.charAt(0));

			if (m1 != null && m2 != null) {
				MarkerPairLikelihood llp = new MarkerPairLikelihood(m1, m2, logLikRatio);
				llpairs.add(llp);
				if (debug) Gpr.debug(llp);
			} else if (debug) {
				if (m1 == null) Gpr.debug("Cannot create marker: " + msaId1 + ", AA sequence '" + seq1.charAt(0) + "'");
				if (m2 == null) Gpr.debug("Cannot create marker: " + msaId2 + ", AA sequence '" + seq2.charAt(0) + "'");
			}

		}

		int tot = countErr + countOk;
		Timer.showStdErr("Genes likelihood file '" + logLikelihoodFile + "'." //
				+ "\n\tEntries loaded: " + count //
				+ "\n\tmapping. Err / OK : " + countErr + " / " + tot + " [ " + (countErr * 100.0 / tot) + "% ]" //
		);
	}

	/**
	 * Read VCF file: Only entries matching markers from GenesLogLik file
	 *
	 * Note:	In order to spread the load on many processes in a cluster, we 'split' the
	 * 			VCF file into 'numSplits' (lineNumber % numSplits = nsplit).
	 * 			In this method, we only store lines that match either of two splits, called
	 * 			'split_i' or 'split_j' (note that actually split_i can be equal
	 * 			to split_j).
	 * 			Then we perform likelihood calculations, only comparing VCF entries
	 * 			in 'split_i' to VCF entries in 'split_j'. This way, we can easily run many
	 * 			processes comparing different splits. E.g. we can set numSplits = 100
	 * 			and launch 100 * (100 / 2 + 1) = 5,100 processes to compare each split.
	 */
	public void readVcf() {
		buildForest();

		// Initialize
		gtsSplitI = new ArrayList<byte[]>(); // Store genotypes for split_i
		gtsSplitJ = new ArrayList<byte[]>(); // Store genotypes for split_j
		gtIdsSplitI = new ArrayList<String>(); // Store VCF 'id' for split_j
		gtIdsSplitJ = new ArrayList<String>(); // Store VCF 'id' for split_j

		// Read VCF file
		int count = 1;
		Timer.showStdErr("Reading vcf file '" + vcfFile + "'");
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		for (VcfEntry ve : vcf) {

			int nsplit = count % numSplits;

			// Do we store this VCF entry?
			if (nsplit == splitI || nsplit == splitJ) {
				// Do we have any MSA in this region?
				if (analyzeAllPairs || !intForest.query(ve).isEmpty()) {

					byte gt[] = ve.getGenotypesScores();
					gt = minorAllele(gt);

					// Did gt pass filtering?
					if (gt != null) {
						String vcfId = ve.getChromosomeName() + ":" + ve.getStart() + "_" + ve.getRef() + "/" + ve.getAltsStr();

						// Store in 'splits'
						if (nsplit == splitI) {
							gtsSplitI.add(gt);
							gtIdsSplitI.add(vcfId);
						}

						if (nsplit == splitJ) {
							gtsSplitJ.add(gt);
							gtIdsSplitJ.add(vcfId);
						}
					}
				}
			}

			// Show something every now and then
			if (count % (SHOW_EVERY_VCF * 100) == 0) System.err.print("\n" + count + " [" + gtsSplitI.size() + " , " + gtsSplitJ.size() + "]\t.");
			else if (count % SHOW_EVERY_VCF == 0) System.err.print('.');

			count++;
		}

		Timer.showStdErr("Done. Total " + count + " VCF entries, added " + (gtsSplitI.size() + gtsSplitJ.size()) + " genotypes [" + gtsSplitI.size() + " , " + gtsSplitJ.size() + "].");
	}

	public void setAnalyzeAllPairs(boolean analyzeAllPairs) {
		this.analyzeAllPairs = analyzeAllPairs;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Show information is being loaded
	 */
	void showCount(boolean ok) {
		if (!debug) {
			int tot = countOk + countErr;

			if (ok) {
				if (tot % SHOW_EVERY_GENES_LL == 0) System.err.print('.'); // Show OK
			} else System.err.print('*'); // Show error

			// Add a newline every now and then
			if (tot % SHOW_LINE_GENES_LL_EVERY == 0) System.err.print("\n" + tot + "\t");
		}
	}
}
