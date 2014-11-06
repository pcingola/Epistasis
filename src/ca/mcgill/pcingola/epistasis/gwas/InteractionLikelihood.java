package ca.mcgill.pcingola.epistasis.gwas;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.RealVector;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.GenotypePos;
import ca.mcgill.pcingola.epistasis.IdMapper;
import ca.mcgill.pcingola.epistasis.IdMapperEntry;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.DistanceResults;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;
import ca.mcgill.pcingola.epistasis.phylotree.EstimateTransitionMatrix;
import ca.mcgill.pcingola.epistasis.phylotree.EstimateTransitionMatrixPairs;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeAa;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix;
import ca.mcgill.pcingola.epistasis.phylotree.UniformTreeValueCache;

/**
 * Calculate epistatic likelihood based on:
 * 		- Multiple sequence alignments
 * 		- Phylogenetic trees
 * 		- Transition matrices Q and Q2
 *
 * @author pcingola
 */
public class InteractionLikelihood {

	public static int MAX_RAND_ITER = 1000;

	boolean debug = false;
	int cpus = -1; // Limit the number of parallel threads?
	double aaFreqs[], aaFreqsContact[];
	String treeFile;
	Random random = new Random(20140716);
	LikelihoodTreeAa tree;
	DistanceResults aaContacts;
	TransitionMatrix Q, Q2;
	MultipleSequenceAlignmentSet msas;
	IdMapper idMapper;
	PdbGenomeMsas pdbGenomeMsas;
	HashMap<Thread, LikelihoodTreeAa> treeNullByThread = new HashMap<Thread, LikelihoodTreeAa>();
	HashMap<Thread, LikelihoodTreeAa> treeAltByThread = new HashMap<Thread, LikelihoodTreeAa>();
	UniformTreeValueCache lcacheNull = new UniformTreeValueCache(GprSeq.AMINO_ACIDS.length);
	UniformTreeValueCache lcacheAlt = new UniformTreeValueCache(GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length);
	Set<String> done = new HashSet<>();

	public InteractionLikelihood(int cpus, String treeFile, DistanceResults aaContacts, TransitionMatrix Q, TransitionMatrix Q2, double aaFreqs[], double aaFreqsContact[], MultipleSequenceAlignmentSet msas, IdMapper idMapper, PdbGenomeMsas pdbGenomeMsas) {
		this.cpus = cpus;
		this.treeFile = treeFile;
		this.aaContacts = aaContacts;
		this.Q = Q;
		this.Q2 = Q2;
		this.aaFreqs = aaFreqs;
		this.aaFreqsContact = aaFreqsContact;
		this.msas = msas;
		this.idMapper = idMapper;
		this.pdbGenomeMsas = pdbGenomeMsas;

		loadTree(treeFile);
	}

	/**
	 * Estimate Q matrix from MSA and Phylogenetic-Tree
	 */
	public void estimateQhat() {
		// Estimate
		EstimateTransitionMatrix mltm = new EstimateTransitionMatrix(tree, msas);
		mltm.setVerbose(true);
		Q = mltm.estimateTransitionMatrix();

		// Sanity checks
		RealVector z = Q.operate(mltm.calcPi());
		System.out.println("Q matrix:\n" + Gpr.prependEachLine("Q_HAT_MATRIX" + "\t", Q));

		System.out.println("Q's Eigenvalues: ");
		double maxLambda = Q.checkEien(true);

		System.out.println("Q_HAT_SUMMARY" //
				+ "\tNorm( Q * pi ):\t" + z.getNorm() //
				+ "\tmax_lambda:\t" + maxLambda //
				+ "\thas_complex_eigenvalues:\t" + Q.hasComplexEigenvalues() //
				+ "\thas_negative_off_diagonal_entries:\t" + Q.hasNegativeOffDiagonalEntries() //
				+ "\tis_zero:\t" + Q.isZero() //
				+ "\tis_symmetric:\t" + Q.isSymmetric() //
		);
	}

	/**
	 * Calculate Qhat2 (400x400 AA-Pairs transition matrix)
	 */
	public void estimateQhat2() {
		// Estimate
		EstimateTransitionMatrixPairs etm = new EstimateTransitionMatrixPairs(tree, msas, aaContacts);
		etm.setVerbose(true);
		Q2 = etm.estimateTransitionMatrix();

		// Sanity checks
		RealVector z = Q2.operate(etm.calcPi());
		System.out.println("Q2 matrix:\n" + Gpr.prependEachLine("Q_HAT2_MATRIX\t", Q2));

		System.out.println("Q's Eigenvalues: ");
		double maxLambda = Q2.checkEien(true);

		System.out.println("Q_HAT2_SUMMARY" //
				+ "\tNorm( Q2 * pi ):\t" + z.getNorm() //
				+ "\tmax_lambda:\t" + maxLambda //
				+ "\thas_complex_eigenvalues:\t" + Q2.hasComplexEigenvalues() //
				+ "\thas_negative_off_diagonal_entries:\t" + Q2.hasNegativeOffDiagonalEntries() //
				+ "\tis_zero:\t" + Q2.isZero() //
				+ "\tis_symmetric:\t" + Q2.isSymmetric() //
		);

	}

	/**
	 * Find all MSAS for all transcripts in 'gene'
	 */
	Set<MultipleSequenceAlignment> findMsasGene(String gene) {
		List<IdMapperEntry> list = idMapper.getByGeneName(gene);
		if (list == null || list.isEmpty()) return null;

		Set<MultipleSequenceAlignment> set = list.stream() //
				.map(ime -> ime.refSeqId) //
				.filter(tid -> msas.getMsasByTrId(tid) != null) //
				.flatMap(tid -> msas.getMsasByTrId(tid).stream()) //
				.collect(Collectors.toSet()) //
		;
		return set;
	}

	/**
	 * Get a tree for the current thread (alt model)
	 */
	LikelihoodTreeAa getTreeAlt() {
		LikelihoodTreeAa currentTree = treeAltByThread.get(Thread.currentThread());

		if (currentTree == null) {
			// No tree? load one
			Timer.showStdErr("Loading phylogenetic tree from " + treeFile);
			currentTree = new LikelihoodTreeAa();
			currentTree.load(treeFile);
			currentTree.setLcache(lcacheAlt);
			treeAltByThread.put(Thread.currentThread(), currentTree);
		}

		return currentTree;
	}

	/**
	 * Get a tree for the current thread (null model)
	 */
	LikelihoodTreeAa getTreeNull() {
		LikelihoodTreeAa currentTree = treeNullByThread.get(Thread.currentThread());

		if (currentTree == null) {
			// No tree? load one
			Timer.showStdErr("Loading phylogenetic tree from " + treeFile);
			currentTree = new LikelihoodTreeAa();
			currentTree.load(treeFile);
			currentTree.setLcache(lcacheNull);
			treeNullByThread.put(Thread.currentThread(), currentTree);
		}

		return currentTree;
	}

	/**
	 * Likelihhod for all AA 'in contact' pairs
	 */
	public void likelihoodAaContacts() {
		// Pre-calculate matrix exponentials
		Timer.showStdErr("Pre-calculating matrix exponentials");
		precalcExps();

		// Calculate likelihoods
		Timer.showStdErr("Calculating likelihoods");
		aaContacts.parallelStream() //
				.filter(d -> msas.getMsa(d.msa1) != null && msas.getMsa(d.msa2) != null) //
				.map(d -> logLikelihoodRatioStr(msas.getMsa(d.msa1), d.msaIdx1, msas.getMsa(d.msa2), d.msaIdx2, false)) //
				.forEach(System.out::println) //
		;
	}

	/**
	 * Likelihood for all AA in geneName Pairs pairs in 'geneNamePairsFile'
	 * File format: "gene1 \t gene2 \n" (spaces added for legibility)
	 */
	public void likelihoodAllAminAcidsInGenes(String genesFile) {
		// Load gene names
		Set<String> geneLines = new HashSet<String>();
		String genesDir = Gpr.dirName(genesFile);
		for (String line : Gpr.readFile(genesFile).split("\n")) {
			String f[] = line.split("\t");
			String gene1 = f[0].trim();
			String gene2 = f[1].trim();

			// Swap?
			if (gene1.compareTo(gene2) > 0) {
				String t = gene1;
				gene1 = gene2;
				gene2 = t;
			}

			geneLines.add(gene1 + "\t" + gene2);
		}
		Timer.showStdErr("Total number of unique pairs: " + geneLines.size());

		// Pre-calculate matrix exponentials
		Timer.showStdErr("Pre-calculating matrix exponentials");
		precalcExps();

		// Calculate likelihoods
		Timer.showStdErr("Calculating likelihood on all pairs");
		Set<String> genes = idMapper.getEntries().stream().map(im -> im.geneName).collect(Collectors.toSet());

		// Pre calculate all MSAs by gene
		HashMap<String, Set<MultipleSequenceAlignment>> msasByGeneName = new HashMap<>();
		for (String gene : genes) {
			Set<MultipleSequenceAlignment> msas = findMsasGene(gene);
			if (msas != null && !msas.isEmpty()) msasByGeneName.put(gene, msas);
		}

		// Calculate likelihoods
		geneLines.parallelStream() //
				.forEach(str -> {
					String f[] = str.split("\t");
					String g1 = f[0], g2 = f[1];
					logLikelihoodGenes(g1, g2, msasByGeneName.get(g1), msasByGeneName.get(g2), genesDir); //
					} //
				);
	}

	/**
	 * Calculate likelihood for the 'alternative model' (H1, i.e. using Qhat2)
	 */
	double likelihoodAltModel(LikelihoodTreeAa tree, MultipleSequenceAlignment msa1, int idx1, MultipleSequenceAlignment msa2, int idx2) {
		// Set sequence and calculate likelihood
		byte seq1[] = msa1.getColumn(idx1);
		byte seq2[] = msa2.getColumn(idx2);
		tree.setLeafSequenceAaPair(seq1, seq2);
		double lik = tree.likelihood(Q2, aaFreqsContact);
		return lik;
	}

	/**
	 * Likelihood 'null distribution
	 */
	public void likelihoodNullModel(int numSamples) {
		// Pre-calculate matrix exponentials
		Timer.showStdErr("Pre-calculating matrix exponentials");
		precalcExps();

		// Calculate likelihoods
		Timer.showStdErr("Calculating likelihood on AA pairs in contact");
		IntStream.range(0, numSamples).parallel() //
				.mapToObj(i -> likelihoodRatioRand()) // Calculate likelihood
				.forEach(System.out::println);
	}

	/**
	 * Calculate likelihood for the 'null model' (H0, i.e. using Qhat)
	 */
	double likelihoodNullModel(LikelihoodTreeAa tree, MultipleSequenceAlignment msa1, int idx1, MultipleSequenceAlignment msa2, int idx2) {
		// Get sequences
		byte seq1b[] = msa1.getColumn(idx1);
		byte seq2b[] = msa2.getColumn(idx2);

		// Set sequence and calculate likelihood
		tree.setLeafSequenceCode(sequenceGaps(seq1b, seq2b));
		double lik1 = tree.likelihood(Q, aaFreqs);

		// Set sequence and calculate likelihood
		tree.setLeafSequenceCode(sequenceGaps(seq2b, seq1b));
		double lik2 = tree.likelihood(Q, aaFreqs);

		double lik = lik1 * lik2;
		return lik;
	}

	/**
	 * Pick two random columns from MSA and calculate likelihood ratio
	 */
	String likelihoodRatioRand() {
		// Pick different transcripts
		MultipleSequenceAlignment msa1, msa2;
		do {
			msa1 = msas.rand(random);
			msa2 = msas.rand(random);
		} while (msa1.getTranscriptId().equals(msa2.getTranscriptId()));

		int msaIdx1 = -1, msaIdx2 = -1;
		for (int i = 0; (msaIdx1 < 0 || msaIdx2 < 0 || msa1.isSkip(msaIdx1) || msa2.isSkip(msaIdx2)) && i < MAX_RAND_ITER; i++) {
			msaIdx1 = random.nextInt(msa1.length());
			msaIdx2 = random.nextInt(msa2.length());
		}

		return logLikelihoodRatioStr(msa1, msaIdx1, msa2, msaIdx2, false);
	}

	/**
	 * Calculate likelihood using all entries in VCF file (within transcript likelihood)
	 */
	public void likelihoodVcf(String vcfFile, PdbGenomeMsas pdbGenomeMsas) {
		// Build interval forest
		msas.buildForest();

		// Pre-calculate matrix exponentials
		Timer.showStdErr("Pre-calculating matrix exponentials");
		precalcExps();

		// Read VCF file
		Timer.showStdErr("Reading VCF file: " + vcfFile);
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		ArrayList<VcfEntry> ves = new ArrayList<VcfEntry>();
		for (VcfEntry ve : vcf)
			ves.add(ve);
		Timer.showStdErr("Done. Added " + ves.size() + " entries.");

		// Process entries 
		ves.stream() //
				.parallel() // 
				.filter(ve -> !msas.query(ve).isEmpty()) // Only use entries that can be mapped
				.map(ve -> new GenotypePos(ve)) // Convert to genotyping position
				.filter(gp -> gp.map2MsaAa(pdbGenomeMsas)) // Successfully mapped to MSA ?
				.forEach(gp -> logLikelihoodGenomicPosVsTranscript(gp)) // Calculate likelihood 
		;
	}

	void loadTree(String phyloFileName) {
		Timer.showStdErr("Loading phylogenetic tree from " + phyloFileName);
		tree = new LikelihoodTreeAa();
		tree.load(phyloFileName);
	}

	/**
	 * Calculate log-likelihood ratio for all MSAs matching gene names
	 * Output format:
	 * 		msa1.id [msaIdx1] \t msa2.id[msaIdx2] \t logLikRatio \t likNull \t likAlt \t seqsStr
	 */
	void logLikelihoodGenes(String gene1, String gene2, Set<MultipleSequenceAlignment> msasGene1, Set<MultipleSequenceAlignment> msasGene2, String genesDir) {
		if (msasGene1 == null || msasGene2 == null) return;

		// Create output directory
		String dir = genesDir + "/likelihood_genes/" + gene1;
		(new File(dir)).mkdirs();

		// Create output file
		String outFile = dir + "/" + gene2 + ".txt";
		String tmpFile = dir + "/" + gene2 + ".tmp";
		if (Gpr.exists(outFile)) {
			Timer.showStdErr("Likelihood genes '" + gene1 + "' and '" + gene2 + "', output file exists, skipping ('" + outFile + "')");
			return;
		}

		try {
			BufferedWriter tmp = new BufferedWriter(new FileWriter(tmpFile));
			Timer.showStdErr("Likelihood genes '" + gene1 + "' and '" + gene2 + "', tmp file: '" + tmpFile + "'");

			// Compare all MSA combinations
			for (MultipleSequenceAlignment msa1 : msasGene1)
				for (MultipleSequenceAlignment msa2 : msasGene2) {
					String key = msa1.getId() + "\t" + msa2.getId();

					// Already calculated? Ignore
					if (!done.contains(key)) {
						done.add(key);
						String lout = logLikelihoodRatio(msa1, msa2, false);
						if (!lout.isEmpty()) tmp.write(lout + "\n");
					}
				}

			tmp.close();

			// We are done, move tmpFile to outFile
			(new File(tmpFile)).renameTo(new File(outFile));
			Timer.showStdErr("Finished likelihood genes '" + gene1 + "' and '" + gene2 + "', output file: '" + outFile + "'");
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Calculate all likelihoods between a genomic position 'gp' and the rest of AAs in the transcript
	 */
	void logLikelihoodGenomicPosVsTranscript(GenotypePos gp) {
		// Find transcript
		MultipleSequenceAlignment msaVcf = msas.getMsa(gp.getMsaId());
		int aaIdxVcf = gp.getAaIdx();
		String trId = msaVcf.getTranscriptId();

		// Get all MSAs for this transcript
		List<MultipleSequenceAlignment> msasTr = msas.getMsasByTrId(trId);

		// Compare position from VCF against ALL other possitions in the transcript
		for (MultipleSequenceAlignment msaTr : msasTr) {
			for (int aaIdxTr = 0; aaIdxTr < msaTr.length(); aaIdxTr++) {
				// Skip same position
				if (aaIdxVcf == aaIdxTr && msaTr.getId().equals(msaVcf.getId())) continue;

				String res = logLikelihoodRatioStr(msaVcf, aaIdxVcf, msaTr, aaIdxTr, true);
				System.out.println(gp.getId() + "\t" + res);
			}
		}
	}

	/**
	 * Calculate likelihood ratio for these entries
	 */
	double logLikelihoodRatio(MultipleSequenceAlignment msa1, int msaIdx1, MultipleSequenceAlignment msa2, int msaIdx2, boolean brief) {
		// Since we execute in parallel, we need one tree per thread
		LikelihoodTreeAa treeNull = getTreeNull();
		LikelihoodTreeAa treeAlt = getTreeAlt();

		// Calculate likelihoods for the null and alternative model
		double likNull = likelihoodNullModel(treeNull, msa1, msaIdx1, msa2, msaIdx2);
		double likAlt = likelihoodAltModel(treeAlt, msa1, msaIdx1, msa2, msaIdx2);
		double logLikRatio = -2.0 * (Math.log(likNull) - Math.log(likAlt));
		return logLikRatio;
	}

	/**
	 * Calculate likelihood for all possible AA pairs within these two MSAs
	 */
	String logLikelihoodRatio(MultipleSequenceAlignment msa1, MultipleSequenceAlignment msa2, boolean brief) {
		System.err.println(msa1.getId() + " (len: " + msa1.length() + "), " + msa2.getId() + " (len: " + msa2.length() + ") = " + (msa1.length() * msa2.length()));

		StringBuilder sb = new StringBuilder();

		for (int i1 = 0; i1 < msa1.length(); i1++) {
			if (msa1.isSkip(i1)) continue;

			for (int i2 = 0; i2 < msa2.length(); i2++) {
				if (msa2.isSkip(i2)) continue;

				String res = logLikelihoodRatioStr(msa1, i1, msa2, i2, brief);
				if (res != null) {
					if (debug) System.err.println(res);
					sb.append(res + "\n");
				}
			}
		}

		return sb.toString();
	}

	public double logLikelihoodRatio(String msaId1, int msaIdx1, String msaId2, int msaIdx2, boolean brief) {
		MultipleSequenceAlignment msa1 = msas.getMsa(msaId1);
		if (msa1 == null) return 0;

		MultipleSequenceAlignment msa2 = msas.getMsa(msaId2);
		if (msa2 == null) return 0;

		return logLikelihoodRatio(msa1, msaIdx1, msa2, msaIdx2, brief);
	}

	/**
	 * Calculate likelihood ratio for these entries
	 */
	String logLikelihoodRatioStr(MultipleSequenceAlignment msa1, int msaIdx1, MultipleSequenceAlignment msa2, int msaIdx2, boolean brief) {
		// Since we execute in parallel, we need one tree per thread
		LikelihoodTreeAa treeNull = getTreeNull();
		LikelihoodTreeAa treeAlt = getTreeAlt();

		// Calculate likelihoods for the null and alternative model
		double likNull = likelihoodNullModel(treeNull, msa1, msaIdx1, msa2, msaIdx2);
		double likAlt = likelihoodAltModel(treeAlt, msa1, msaIdx1, msa2, msaIdx2);
		double logLikRatio = -2.0 * (Math.log(likNull) - Math.log(likAlt));

		// Return results
		String seqsStr = "";
		if (!brief) {
			String seq1 = msa1.getColumnString(msaIdx1);
			String seq2 = msa2.getColumnString(msaIdx2);
			seqsStr = "\t" + seq1 + "\t" + seq2;
		}

		return msa1.getId() + "[" + msaIdx1 + "]\t" + msa2.getId() + "[" + msaIdx2 + "]"//
				+ "\t" + logLikRatio //
				+ "\t" + likNull //
				+ "\t" + likAlt //
				+ seqsStr //
		;
	}

	public String logLikelihoodRatioStr(String msaId1, int msaIdx1, String msaId2, int msaIdx2, boolean brief) {
		MultipleSequenceAlignment msa1 = msas.getMsa(msaId1);
		if (msa1 == null) return null;

		MultipleSequenceAlignment msa2 = msas.getMsa(msaId2);
		if (msa2 == null) return null;

		return logLikelihoodRatioStr(msa1, msaIdx1, msa2, msaIdx2, brief);
	}

	/**
	 * Pre-calculate matrix exponential
	 */
	public void precalcExps() {

		// Find all times in the tree
		HashSet<Double> times = new HashSet<>();
		tree.times(times);

		// Pre-calculate Q's exponentials
		times.parallelStream() //
				.peek(t -> System.err.println("Matrix\tdim:" + Q.getRowDimension() + "x" + Q.getColumnDimension() + "\tExp(" + t + ")")) //
				.forEach(t -> Q.matrix(t)) //
		;

		// Calculate all gene-gene
		(cpus == 1 ? times.stream() : times.parallelStream()) //
				.peek(t -> System.err.println("Matrix\tdim:" + Q2.getRowDimension() + "x" + Q2.getColumnDimension() + "\tExp(" + t + ")")) //
				.forEach(t -> Q2.matrix(t)) //
		;
	}

	/**
	 * Create a copy of 'seq1' having gaps present in either seq1 or seq2
	 */
	int[] sequenceGaps(byte seq1[], byte seq2[]) {
		int c[] = new int[seq1.length];

		for (int i = 0; i < seq1.length; i++)
			if (seq1[i] < 0 || seq2[i] < 0) c[i] = -1;
			else c[i] = seq1[i];

		return c;
	}

}
