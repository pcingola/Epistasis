package ca.mcgill.pcingola.epistasis.gwas;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.Genotype;
import ca.mcgill.pcingola.epistasis.coordinates.GenomicCoordinates;
import ca.mcgill.pcingola.epistasis.likelihood.InteractionLikelihood;
import ca.mcgill.pcingola.epistasis.likelihood.LikelihoodAnalysisGtPair;
import ca.mcgill.pcingola.epistasis.likelihood.MarkerPairLikelihood;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;

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
	public static double LL_THRESHOLD_LOGREG = 6.0;
	public static double LL_THRESHOLD_MSA = 0.0;
	public static double LL_THRESHOLD_TOTAL = 5.0;

	// Splits information
	boolean analyzeAllPairs = false; // Use for testing and debugging
	boolean debug = false;
	boolean verbose = false;
	int splitI = 2;
	int splitJ = 3;
	int numSplits = 100;
	int countOk, countErr;
	double llThresholdLogReg = LL_THRESHOLD_LOGREG;
	double llThresholdMsa = LL_THRESHOLD_MSA;
	String configFile;
	String logLikelihoodFile; // Log likelihood file (epistatic model)
	String vcfFile;
	String genomeVer, pdbDir;
	String phenoCovariatesFile;
	List<MarkerPairLikelihood> llpairs; // Gene log-likelihood entries
	List<Genotype> gtsSplitI, gtsSplitJ;
	Map<String, Transcript> trancriptById; // Transcript by (incomplete) transcript ID (no version number is used)
	Map<String, Marker> llmarkerById = new HashMap<String, Marker>(); // log-likelihood markers by ID
	Map<Long, LikelihoodAnalysisGtPair> llAnByThreadId = new HashMap<>();
	AutoHashMap<String, ArrayList<byte[]>> gtById; // Genotypes by ID
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

	public int getCountErr() {
		return countErr;
	}

	public int getCountOk() {
		return countOk;
	}

	/**
	 * Get a likelihood analysis object for each thread
	 */
	LikelihoodAnalysisGtPair getLikelihoodAnalysis2() {
		long threadId = Thread.currentThread().getId();

		LikelihoodAnalysisGtPair llan = llAnByThreadId.get(threadId);

		if (llan == null) {
			llan = new LikelihoodAnalysisGtPair(phenoCovariatesFile, vcfFile);
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
		for (int idxi = 0; idxi < gtsSplitI.size(); idxi++) {

			// Split_i info
			final int i = idxi;
			Genotype gti = gtsSplitI.get(idxi);

			int minJ = 0;
			if (splitI == splitJ) minJ = idxi + 1;

			// Parallel on split_j
			IntStream.range(minJ, gtsSplitJ.size()) //
					.parallel() //
					.forEach(j -> {
						GwasResult gwasRes = gwas(gti, gtsSplitJ.get(j));
						double llTot = gwasRes.logLik();
						if (llTot > llThresholdLogReg) countLl.inc();
						if (llTot != 0.0) Timer.show(count.inc() + " (" + i + " / " + j + ")\t" + countLl + "\t" + gwasRes);
					});
		}
	}

	/**
	 * Perform analysis on genotypes 'i' and 'j'
	 */
	GwasResult gwas(Genotype genoi, Genotype genoj) {
		//---
		// Likelihood based on logistic regression
		//---
		LikelihoodAnalysisGtPair llan = getLikelihoodAnalysis2();
		GwasResult gwasRes = llan.logLikelihood(genoi, genoj);

		// Log likelihood form logistic regression is too low?
		// => Don't bother to calculate next part
		if (gwasRes.logLikelihoodRatioLogReg < llThresholdLogReg) return gwasRes;

		// Calculate p-value form logistic regression likelihood test
		gwasRes.pvalueLogReg();

		//---
		// Likelihood based on interaction
		//---

		// Find corresponding MSA ID and index for both genotypes
		genoi.mapGenomic2Msa(pdbGenomeMsas);
		if (!genoi.hasMsaInfo()) return gwasRes;

		genoj.mapGenomic2Msa(pdbGenomeMsas);
		if (!genoj.hasMsaInfo()) return gwasRes;

		// Likelihood based on epistatic interaction
		String msaId1 = genoi.getMsaId(), msaId2 = genoj.getMsaId();
		int msaIdx1 = genoi.getAaIdx(), msaIdx2 = genoj.getAaIdx();
		interactionLikelihood.logLikelihoodRatio(msaId1, msaIdx1, msaId2, msaIdx2, gwasRes);

		// Epistatic likelihood model too low? => Don't bother to calculate next part
		if (gwasRes.logLik() < LL_THRESHOLD_TOTAL && gwasRes.logLikelihoodRatioMsa < llThresholdMsa) return gwasRes;

		//---
		// Calculate Bayes Factor using Laplace approximation maethod
		//---
		double h1 = 1.0; // P(theta_1 | M_1) : This is the a-priory distribution
		double h0 = 1.0; // P(theta_0 | M_0) : This is the a-priory distribution
		gwasRes.bayesFactor(h1, h0);

		return gwasRes;
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
		Chromosome chromo = pdbGenomeMsas.getConfig().getGenome().getOrCreateChromosome(chr);
		GenomicCoordinates gp = new GenomicCoordinates(chromo, start, end, id);
		if (gp.markerTrAaIdx(pdbGenomeMsas, trId, aaIdx, aaExpected)) {
			countOk++;
			showCount(true);
			llmarkerById.put(id, gp); // Cache marker
			return gp;
		} else {
			countErr++;
			showCount(false);
			return null;
		}
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
		MultipleSequenceAlignmentSet msas = pdbGenomeMsas.getMsas();
		msas.buildForest();

		// Initialize
		gtsSplitI = new ArrayList<Genotype>(); // Store genotypes for split_i
		gtsSplitJ = new ArrayList<Genotype>(); // Store genotypes for split_j

		// Read VCF file
		int count = 1;
		Timer.showStdErr("Reading vcf file '" + vcfFile + "'");
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		for (VcfEntry ve : vcf) {

			int nsplit = count % numSplits;

			// Do we store this VCF entry?
			if (nsplit == splitI || nsplit == splitJ) {
				// Do we have any MSA in this region?
				if (analyzeAllPairs || !msas.query(ve).isEmpty()) {

					Genotype geno = new Genotype(ve);

					// Do we have at least a few alleles? Store in 'splits'
					if (geno.getMinorAlleleCount() >= MINOR_ALLELE_COUNT) {
						if (nsplit == splitI) gtsSplitI.add(geno);
						if (nsplit == splitJ) gtsSplitJ.add(geno);
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
