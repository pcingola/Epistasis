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
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Perform GWAS using epistasis data
 *
 * @author pcingola
 */
public class GwasEpistasis extends SnpEff {

	public static int SHOW_EVERY_VCF = 1000;
	public static int SHOW_EVERY_GENES_LL = 10000;
	public static int SHOW_LINE_GENES_LL_EVERY = 100 * SHOW_EVERY_GENES_LL;
	public static double SHOW_LINE_LL_MIN = 1.0;
	public static int MINOR_ALLELE_COUNT = 5;
	public static double LL_THRESHOLD = 6.0;

	// Splits information
	boolean analyzeAllPairs = false; // Use for testing and debugging
	int splitI = 2;
	int splitJ = 3;
	int numSplits = 30;
	int countOk, countErr;
	double llThreshold = LL_THRESHOLD;
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
	IntervalForest llforest; // Interval forest of ll-markers

	public GwasEpistasis(String configFile, String genomeVer, String genesLikeFile, String vcfFile, String phenoCovariatesFile) {
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		logLikelihoodFile = genesLikeFile;
		this.vcfFile = vcfFile;
		this.phenoCovariatesFile = phenoCovariatesFile;
	}

	/**
	 * Build interval forest using LogLik markers
	 */
	protected void buildForest() {
		Timer.showStdErr("Building Log-likelihood marker forest");
		Markers markers = new Markers();
		markers.addAll(llmarkerById.values());

		// Create forest and build
		llforest = new IntervalForest(markers);
		llforest.build();
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
	 * Perform GWAS analysis using epistatic information
	 */
	public void gwas() {
		initialize();
		readGenesLogLikelihood();
		readVcf();

		if (analyzeAllPairs) gwasAll();
		else gwasMatching();
	}

	/**
	 * Analyze all pairs (mostly for debugging)
	 */
	void gwasAll() {
		for (String idi : gtById.keySet()) {
			for (String idj : gtById.keySet()) {
				if (idi.compareTo(idj) < 0) {
					LikelihoodAnalysis2 llan = getLikelihoodAnalysis2();
					for (byte gti[] : gtById.get(idi))
						for (byte gtj[] : gtById.get(idj)) {
							double ll = llan.logLikelihood(idi, gti, idj, gtj);
							System.out.println("ll:" + ll + "\t" + idi + "\t" + idj);
						}
				}
			}
		}
	}

	/**
	 *  Analyze pairs in VCF file that match enriched region (form epistatic analysis)
	 */
	void gwasMatching() {
		llpairs.stream() //
				.parallel() //
				.forEach(llpair -> {

					// Find genotypes in under markers
						String idi = llpair.getMarker1().getId();
						String idj = llpair.getMarker2().getId();

						// No genotypes in any of those regions? Nothing to do
						if (!gtById.containsKey(idi) || !gtById.containsKey(idj)) {
							if (debug) Gpr.debug("Nothing found:\t" + llpair);
						} else {
							LikelihoodAnalysis2 llan = getLikelihoodAnalysis2();

							// Analyze all genotype pairs within those regions
							for (byte gti[] : gtById.get(idi)) {
								for (byte gtj[] : gtById.get(idj)) {
									double ll = llan.logLikelihood(idi, gti, idj, gtj);
									if (debug || ll > SHOW_LINE_LL_MIN) System.out.println("ll:" + ll + "\t" + idi + "\t" + idj);
								}
							}
						}
					} //
				);
	}

	/**
	 * Load all data
	 */
	public void initialize() {
		// Initialize SnpEff
		String argsSnpEff[] = { "eff", "-v", "-c", configFile, genomeVer };
		args = argsSnpEff;
		setGenomeVer(genomeVer);
		parseArgs(argsSnpEff);
		loadConfig();

		// Load SnpEff database
		if (genomeVer != null) loadDb();

		// Initialize trancriptById
		trancriptById = new HashMap<>();
		for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes())
			for (Transcript tr : g) {
				String id = tr.getId();
				if (id.indexOf('.') > 0) {
					// When using RefSeq transcripts, we don't store sub-version number
					id = id.substring(0, id.indexOf('.'));
				}
				trancriptById.put(id, tr);
			}

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
		int idx = Gpr.parseIntSafe(f[5]);

		//---
		// Calculate position within CDS
		//---

		// Find transcript and exon
		Transcript tr = trancriptById.get(trId);
		if (tr == null) return null;
		Exon ex = tr.findExon(start);
		if (ex == null) return null;

		// Calculate start position
		int startPos;
		int fr = 0;
		if (ex.getFrame() != 0) {
			if (ex.isStrandPlus()) {
				idx--;
				if (ex.getFrame() == 2) idx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			} else {
				idx--;
				if (ex.getFrame() == 2) idx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			}
		}

		// Find AA start position
		if (ex.isStrandPlus()) {
			int exStart = Math.max(start, tr.getCdsStart());
			startPos = exStart + (idx * 3 + fr);
		} else {
			int exEnd = Math.min(end, tr.getCdsStart());
			startPos = exEnd - (idx * 3 + fr);
		}

		// Get position within CDS
		int cdsBase = tr.baseNumberCds(startPos, false);
		int cds2pos[] = tr.baseNumberCds2Pos();
		if ((ex.isStrandPlus() && (startPos < ex.getStart())) //
				|| (ex.isStrandMinus() && (startPos > ex.getEnd()))) {
			// If the position is outside the exon, then we must jump to previous exon
			startPos = cds2pos[cdsBase - ex.getFrame()];
			cdsBase = tr.baseNumberCds(startPos, true);
		}

		//---
		// Sanity check: Make sure that AA matches between transcript model and MSA data from 'genes likelihood' file
		//---

		// Extract codon
		String cdsSeq = tr.cds();
		String codonStr = cdsSeq.substring(cdsBase, cdsBase + 3);
		String aa = genome.codonTable().aa(codonStr);

		if (aa.equals("" + aaExpected)) {
			countOk++;
			if (debug) Gpr.debug("OK: " + id + " : " + aa);
			else showCount(true);
		} else {
			countErr++;
			if (debug) Gpr.debug("Entry ID     : " + id //
					+ "\ntr ID        : " + trId + ", chr: " + chr + ", start: " + start + ", end: " + end + ", idx: " + idx + ", fr: " + fr//
					+ "\nTranscript : " + tr //
					+ "\nExon       : " + ex //
					+ "\nStart pos: " + startPos //
					+ "\nCodon    : " + codonStr + ", aa (real): " + aa + ", aa (exp): " + aaExpected //
			);
			else showCount(false);
		}

		//---
		// Create marker
		// Important: The marker has ALL bases in the codon.
		//            For instance, is the codon is split between two exons, the
		//            marker will contain the intron
		//---
		int markerStart, markerEnd;
		if (tr.isStrandPlus()) {
			markerStart = cds2pos[cdsBase];
			markerEnd = cds2pos[cdsBase + 2];
		} else {
			markerStart = cds2pos[cdsBase + 2];
			markerEnd = cds2pos[cdsBase];
		}

		Chromosome chromo = genome.getChromosome(chr);
		marker = new Marker(chromo, markerStart, markerEnd, ex.isStrandMinus(), id);
		llmarkerById.put(id, marker); // Cache marker
		return marker;
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

	/**
	 * Read VCF file: Only entries matching markers from GenesLogLik file
	 * TODO: We could optimize this by using an index and reading only the regions we need
	 */
	public void testVcf() {
		// Read VCF file
		readVcf();

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
						LikelihoodAnalysis2 llan = getLikelihoodAnalysis2();
						String idj = gtIdsSplitJ.get(j);
						double ll = llan.logLikelihood(idi, gti, idj, gtsSplitJ.get(j));

						if (ll > llThreshold) countLl.inc();
						if (ll != 0.0) Timer.show(count.inc() + " (" + i + " / " + j + ")\t" + countLl + "\t" + ll + "\t" + idi + "\t" + idj);
					} //
					);
		}
	}
}
