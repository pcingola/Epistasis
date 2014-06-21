package ca.mcgill.pcingola.epistasis;

import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

import org.apache.commons.math3.util.Pair;

import ca.mcgill.mcb.pcingola.snpEffect.commandLine.CommandLine;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.entropy.Entropy;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTree;
import ca.mcgill.pcingola.epistasis.phylotree.MaxLikelihoodTm;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix;

/**
 * Main command line
 *
 * @author pcingola
 */
public class Epistasis implements CommandLine {

	public static int MIN_DISTANCE = 1000000;
	public static boolean debug = true;

	public static void main(String[] args) {
		Epistasis epistasis = new Epistasis(args);
		epistasis.run();
	}

	boolean nextProt;
	String[] args;
	String cmd;
	String aaContactFile, configFile, genome, idMapFile, multAlignFile, pdbDir, qMatrixFile, treeFile;
	LikelihoodTree tree;
	DistanceResults aaContacts;
	TransitionMatrix Q;
	MultipleSequenceAlignmentSet msas;
	MaxLikelihoodTm mltm;
	IdMapper idMapper;
	PdbGenome pdbGenome;

	public Epistasis(String[] args) {
		this.args = args;
	}

	@Override
	public String[] getArgs() {
		return args;
	}

	/**
	 * Is this sequence fully conserved?
	 */
	boolean isFullyConserved(String seq) {
		char prevBase = '-';
		for (int i = 0; i < seq.length(); i++) {
			char base = seq.charAt(i);
			if (base == '-') continue;
			if (prevBase != '-' && prevBase != base) return false;
			prevBase = base;
		}
		return true;
	}

	/**
	 * Are all sequences fully conserved?
	 */
	boolean isFullyConserved(String[] seqs) {
		for (int i = 0; i < seqs.length; i++)
			if (seqs[i] != null && !isFullyConserved(seqs[i])) return false;
		return true;
	}

	/**
	 * Load files
	 */
	void load() {
		if (treeFile != null) loadTree(treeFile);
		if (multAlignFile != null) loadMsas(multAlignFile);
		if (idMapFile != null) loadIdMap(idMapFile);
		if (aaContactFile != null) loadAaContact(aaContactFile);

		if (genome != null) {
			pdbGenome = new PdbGenome(configFile, genome, pdbDir);
			pdbGenome.setDebug(debug);
			pdbGenome.setIdMapper(idMapper);
			pdbGenome.setMsas(msas);
			pdbGenome.setTree(tree);
			pdbGenome.setNextProt(nextProt);
			pdbGenome.initialize();
		}

		if (qMatrixFile != null) loadQ(qMatrixFile);
	}

	/**
	 * Load AA contact list
	 */
	List<DistanceResult> loadAaContact(String aaContactFile) {
		Timer.showStdErr("Loading AA contact information " + aaContactFile);
		aaContacts = new DistanceResults();
		aaContacts.load(aaContactFile);
		return aaContacts;
	}

	void loadIdMap(String idMapFile) {
		Timer.showStdErr("Loading id maps " + idMapFile);
		idMapper = new IdMapper(idMapFile);
	}

	/**
	 * Load MSAs
	 */
	void loadMsas(String multAlign) {
		int numAligns = tree.childNames().size();

		// Load: MSA
		Timer.showStdErr("Loading " + numAligns + " way multiple alignment from " + multAlign);
		msas = new MultipleSequenceAlignmentSet(multAlign, numAligns);
		msas.load();
	}

	/**
	 * Load Q matrix
	 */
	void loadQ(String qMatrixFile) {
		Timer.showStdErr("Loading Q matrix  from file '" + qMatrixFile);
		if (!Gpr.exists(qMatrixFile)) {
			System.err.println("Matrix file doesn't exists, nothing done");
			return;
		}

		MaxLikelihoodTm mltm = new MaxLikelihoodTm(tree, msas);
		Q = mltm.loadTransitionMatrix(qMatrixFile);
	}

	/**
	 * Load a tree from a file
	 */
	void loadTree(String phyloFileName) {
		Timer.showStdErr("Loading phylogenetic tree from " + phyloFileName);
		tree = new LikelihoodTree();
		tree.load(phyloFileName);
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length < 1) usage("Missing command");
		cmd = args[0];

		int argNum = 1;
		String type = "";

		switch (cmd.toLowerCase()) {
		case "addmsaseqs":
			configFile = args[argNum++];
			genome = args[argNum++];
			treeFile = args[argNum++];
			multAlignFile = args[argNum++];
			idMapFile = args[argNum++];
			aaContactFile = args[argNum++];
			runAddMsaSeqs();
			break;

		case "aacontactstats":
			type = args[argNum++];
			aaContactFile = args[argNum++];
			runAaContactStats(type);
			break;

		case "background":
			type = args[argNum++];
			int numBases = Gpr.parseIntSafe(args[argNum++]);
			int numSamples = Gpr.parseIntSafe(args[argNum++]);
			treeFile = args[argNum++];
			multAlignFile = args[argNum++];
			if (numBases <= 0) usage("number of bases must be positive number");
			if (numSamples <= 0) usage("number of samples must be positive number");
			runBackground(type, numBases, numSamples);
			break;

		case "conservation":
			treeFile = args[argNum++];
			multAlignFile = args[argNum++];
			aaContactFile = args[argNum++];
			runConservation();
			break;

		case "mappdbgenome":
			configFile = args[argNum++];
			genome = args[argNum++];
			pdbDir = args[argNum++];
			idMapFile = args[argNum++];
			runMapPdbGenome();
			break;

		case "nextprot":
			configFile = args[argNum++];
			genome = args[argNum++];
			idMapFile = args[argNum++];
			aaContactFile = args[argNum++];
			runNextProt();
			break;

		case "pdbdist":
			// Parse command line
			double distThreshold = Gpr.parseDoubleSafe(args[argNum++]);
			if (distThreshold <= 0) usage("Distance must be a positive number: '" + args[argNum - 1] + "'");
			int aaMinSeparation = Gpr.parseIntSafe(args[argNum++]);
			pdbDir = args[argNum++];
			idMapFile = args[argNum++];
			runPdbDist(distThreshold, aaMinSeparation);
			break;

		case "qhat":
			if (args.length < 4) usage("Missing arguments for command '" + cmd + "'");
			treeFile = args[argNum++];
			multAlignFile = args[argNum++];
			qMatrixFile = args[argNum++];
			runQhat();
			break;

		case "test":
			treeFile = args[argNum++];
			multAlignFile = args[argNum++];
			aaContactFile = args[argNum++];
			runTest();
			break;

		default:
			throw new RuntimeException("Unknown command: '" + cmd + "'");
		}

		Timer.showStdErr("Done!");
	}

	/**
	 * Calculate or load transition matrix
	 */
	TransitionMatrix qHat(String qMatrixFile) {
		MaxLikelihoodTm mltm = new MaxLikelihoodTm(tree, msas);
		Timer.showStdErr("Q matrix file '" + qMatrixFile + "' not found, calculating matrix");
		Q = mltm.estimateTransitionMatrix();
		Timer.showStdErr("Saving Q matrix to file '" + qMatrixFile + "'");
		Q.save(qMatrixFile);

		System.out.println("Q matrix:\n" + Q.toString());
		mltm.showEienQ();

		return Q;
	}

	/**
	 * Parse and Dispatch to right command
	 */
	@Override
	public boolean run() {
		Timer.showStdErr("Start");
		parseArgs(args);
		Timer.showStdErr("End");
		return true;
	}

	/**
	 * Calculate MI and other statistics
	 */
	void runAaContactStats(String type) {
		// Select which function to use
		Function<DistanceResult, Double> f = selectFunction(type);
		load();

		//---
		// Group by genomic position
		//---
		Timer.showStdErr("Sort by position");
		DistanceResults aaContactsUniq = new DistanceResults();
		aaContacts.stream() //
				.filter(d -> !d.aaSeq1.isEmpty() && !d.aaSeq2.isEmpty()) // Filter out empty sequences
				.forEach(d -> aaContactsUniq.collectMin(d, d.toStringPos()));
		aaContactsUniq.addMins(); // Move 'best' results from hash to list

		//---
		// Show MI and conservation
		//---
		aaContactsUniq.stream() //
				.forEach( //
						d -> System.out.printf("%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" //
								, d //
								, Entropy.mutualInformation(d.aaSeq1, d.aaSeq2) //
								, Entropy.entropy(d.aaSeq1, d.aaSeq2) //
								, Entropy.variationOfInformation(d.aaSeq1, d.aaSeq2) //
								, Entropy.condEntropy(d.aaSeq1, d.aaSeq2) //
								, Entropy.condEntropy(d.aaSeq2, d.aaSeq1) //
								, Entropy.entropy(d.aaSeq1) //
								, Entropy.entropy(d.aaSeq2) //
								, MsaSimilarity.conservation(d.aaSeq1) //
								, MsaSimilarity.conservation(d.aaSeq2) //
								) //
				);

		//---
		// Count first 'AA' (all)
		//---
		CountByType countFirstAaAll = new CountByType();
		aaContactsUniq.stream().forEach(d -> countFirstAaAll.inc(d.getAaPair()));
		System.err.println("Count fist AA (all):\n" + countFirstAaAll.toStringSort());

		//---
		// Count first 'AA' (not-fully conserved)
		//---
		CountByType countFirstAa = new CountByType();
		aaContactsUniq.stream() //
				.filter(d -> MsaSimilarity.conservation(d.aaSeq1) < 1.0 && MsaSimilarity.conservation(d.aaSeq2) < 1.0) // Do not calculate on fully conserved sequences (entropy is zero)
				.forEach(d -> countFirstAa.addScore(d.getAaPair(), f.apply(d))) //
		;
		System.err.println("Count fist AA (non-fully conserved) " + type + " :\n" + countFirstAa.toStringSort());

		//---
		// Count first 'AA' with annotations (all)
		//---
		CountByType countFirstAaAnnAll = new CountByType();
		aaContactsUniq.stream() //
				.filter(d -> !d.annotations1.isEmpty() && !d.annotations2.isEmpty()) // Only entries having annotations
				.forEach( //
						d -> d.getAaPairAnnotations().forEach(ap -> countFirstAaAnnAll.inc(ap)) //
				) //
		;
		System.err.println("Count fist AA with annotations (all):\n" + countFirstAaAnnAll.toStringSort());

		//---
		// Count first 'AA' with annotations (not-fully conserved)
		//---
		CountByType countFirstAaAnn = new CountByType();
		aaContactsUniq.stream() //
				.filter(d -> !d.annotations1.isEmpty() && !d.annotations2.isEmpty()) // Only entries having annotations
				.filter(d -> MsaSimilarity.conservation(d.aaSeq1) < 1.0 && MsaSimilarity.conservation(d.aaSeq2) < 1.0) // Do not calculate on fully conserved sequences (entropy is zero)
				.forEach( //
						d -> d.getAaPairAnnotations().forEach( // Add to all annotation pairs
								ap -> countFirstAaAnn.addScore(ap, f.apply(d)) //
								) //
				) //
		;
		System.err.println("Count fist AA with annotations (non-fully conserved), " + type + " :\n" + countFirstAaAnn.toStringSort());

	}

	/**
	 * Add MSA sequences to 'AA contact' data
	 */
	void runAddMsaSeqs() {
		load();

		// Sanity check: Make sure MSA protein sequences match genome's protein data
		Timer.showStdErr("Checking MSA proteing sequences vs. genome protein sequences");
		pdbGenome.checkSequenceMsaTr();
		System.err.println("Totals:\n" + pdbGenome.countMatch);
		pdbGenome.resetStats();

		// Add MSA sequences to 'AA contact' entries
		Timer.showStdErr("Adding MSA sequences");
		aaContacts.forEach(d -> pdbGenome.mapToMsa(msas, d));
		System.err.println("Totals:\n" + pdbGenome.countMatch);

		System.err.println("Mapped AA sequences:\n");
		aaContacts.stream().filter(d -> d.aaSeq1 != null).forEach(System.out::println);
	}

	/**
	 * Run Mutual Information
	 */
	void runBackground(String type, int numBases, int numSamples) {
		// Load
		load();

		// Select which function to use
		MsaSimilarity sim = null;
		switch (type) {
		case "mi":
			sim = new MsaSimilarityMutInf(msas);
			break;

		case "miN":
			sim = new MsaSimilarityMutInfN(msas, numBases);
			break;

		case "varInf":
			sim = new MsaDistanceVarInf(msas);
			break;

		default:
			throw new RuntimeException("Unknown type '" + type + "'");
		}

		// Run
		sim.backgroundDistribution(numSamples);

		// Show distribution
		System.err.println(sim);
	}

	void runBayes() {
		//	// Calculate Bayes Factor
		//
		//	// For each MSA...
		//	mltm.calcPi();
		//	for (MultipleSequenceAlignment msa1 : msas) {
		//		for (MultipleSequenceAlignment msa2 : msas) {
		//			// Make sure the MSAs are far away from each other
		//			if (msa2.compareTo(msa1) < MIN_DISTANCE) continue;
		//
		//			// All pair of positions
		//			for (int pos1 = 0; pos1 < msa1.length(); pos1++) {
		//				if (msa1.isSkip(pos1)) continue;
		//
		//				// Log likelihood at position 1
		//				double loklik1 = mltm.logLikelyhood(msa1, pos1);
		//
		//				for (int pos2 = 0; pos2 < msa2.length(); pos2++) {
		//					if (msa2.isSkip(pos2)) continue;
		//
		//					// Log likelihood at position 2
		//					double loklik2 = mltm.logLikelyhood(msa2, pos2);
		//					System.out.println(msa1.getId() + "\t" + msa2.getId() + "\t" + pos1 + "\t" + pos2 + "\t" + loklik1 + "\t" + loklik2);
		//
		//					// Combined log likelihood
		//					mltm.logLikelyhood(msa1, pos1, msa2, pos2);
		//				}
		//			}
		//		}
		//	}
	}

	/**
	 * Conservation statistics
	 */
	void runConservation() {
		load();

		// N bases
		System.out.println("Window\tTotal_bases\tConserved\tConserved%");
		for (int n = 1; n < 10; n++) {
			int num = n;
			Counter total = new Counter();
			Counter conserved = new Counter();

			// Count for all MSAs
			msas.getMsas().stream() //
					.filter(msa -> msa.length() > num) //
					.forEach(msa -> IntStream.range(0, msa.length() - num) //
							.peek(i -> total.inc()) //
							.filter(i -> msa.isFullyConserved(i, num)) //
							.forEach(i -> conserved.inc()) //
					);

			// Count number of 'fully conserved' AA windows of n-columns around 'AA in contact' 
			Counter totalIc = new Counter();
			Counter conservedIc = new Counter();
			aaContacts.stream() //
					.filter(d -> !d.msa1.isEmpty() && !d.msa2.isEmpty() && msas.getMsa(d.msa1) != null) //
					.peek(d -> totalIc.inc()) //
					.map(d -> new Pair<String[], String[]>(msas.findColSequences(d.msa1, d.msaIdx1, num), msas.findColSequences(d.msa2, d.msaIdx2, num))) //
					.filter(p -> isFullyConserved(p.getFirst()) && isFullyConserved(p.getSecond())) //
					.forEach(d -> conservedIc.inc()) //
			;

			System.out.println(num //
					+ "\t" + total //
					+ "\t" + conserved //
					+ "\t" + (100.0 * conserved.get()) / total.get() + " %" //
					+ "\t" + totalIc //
					+ "\t" + conservedIc //
			);
		}

	}

	void runMapPdbGenome() {
		load();
		pdbGenome.checkSequencePdbTr();
	}

	void runNextProt() {
		nextProt = true;
		load();

		// Add nextprot annotations
		aaContacts.forEach(d -> pdbGenome.nextProt(d));
	}

	/**
	 * Run Pdb distance
	 */
	void runPdbDist(double distThreshold, int aaMinSeparation) {
		load();

		// Run analysis
		PdbDistanceAnalysis pdban = new PdbDistanceAnalysis(pdbDir, distThreshold, aaMinSeparation, idMapper);
		pdban.run();
		System.out.println(pdban);

		// Write results
		String outFile = "pdb_distance_by_AA_pos.txt";
		Gpr.toFile(outFile, pdban);
		System.err.println("Distance metrics file written to: " + outFile);
	}

	/**
	 * Estimate Q matrix from MSA and Phylogenetic-Tree
	 */
	void runQhat() {
		load();
		sanityCheck(tree, msas); // Sanity check: Make sure that the alignment and the tree match
		qHat(qMatrixFile); // Calculate Qhat
	}

	/**
	 * Test
	 */
	void runTest() {
		load();

		// Count number of 'fully conserved' AA windows of n-columns around 'AA in contact' 
		Counter total = new Counter();
		Counter conserved = new Counter();
		int n = 1; // Number of surrounding columns
		aaContacts.stream()//
				.filter(d -> !d.msa1.isEmpty() && !d.msa2.isEmpty() && msas.getMsa(d.msa1) != null) //
				.peek(System.out::println) //
				.peek(d -> total.inc()) //
				.map(d -> new Pair<String[], String[]>(msas.findColSequences(d.msa1, d.msaIdx1, n), msas.findColSequences(d.msa2, d.msaIdx2, n))) //
				.filter(p -> isFullyConserved(p.getFirst()) && isFullyConserved(p.getSecond())) //
				.peek(d -> conserved.inc()) //
		;

	}

	/**
	 * Check consistency between MSA and tree
	 */
	void sanityCheck(LikelihoodTree tree, MultipleSequenceAlignmentSet msas) {
		// Sanity check: Make sure that the alignment and the tree match
		StringBuilder sbMsa = new StringBuilder();
		for (String sp : msas.getSpecies())
			sbMsa.append(sp + " ");

		StringBuilder sbTree = new StringBuilder();
		for (String sp : tree.childNames())
			sbTree.append(sp + " ");

		if (!sbTree.toString().equals(sbMsa.toString())) throw new RuntimeException("Species form MSA and Tree do not match:\n\tMSA : " + sbMsa + "\n\tTree: " + sbTree);

		System.out.println("\nSpecies [" + tree.childNames().size() + "]: " + sbTree);
	}

	// Select which function to use
	Function<DistanceResult, Double> selectFunction(String type) {
		switch (type) {
		case "mi":
			return d -> Entropy.variationOfInformation(d.aaSeq1, d.aaSeq2);

		case "varInf":
			return d -> Entropy.mutualInformation(d.aaSeq1, d.aaSeq2);

		default:
			throw new RuntimeException("Unknown type '" + type + "'");
		}
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: " + this.getClass().getSimpleName() + " cmd options");

		System.err.println("Command 'aaContactStats' : " + this.getClass().getSimpleName() + " aaContactStats type aa_contact.nextprot.txt ");
		System.err.println("Command 'addMsaSeqs'     : " + this.getClass().getSimpleName() + " addMsaSeqs snpeff.config genome phylo.nh multiple_alignment_file.fa id_map.txt aa_contact.txt ");
		System.err.println("Command 'background'     : " + this.getClass().getSimpleName() + " background number_of_bases number_of_samples phylo.nh multiple_alignment_file.fa");
		System.err.println("Command 'conservation'   : " + this.getClass().getSimpleName() + " conservation phylo.nh multiple_alignment_file.fa");
		System.err.println("Command 'corr'           : " + this.getClass().getSimpleName() + " corr phylo.nh multiple_alignment_file.fa");
		System.err.println("Command 'mapPdbGenome'   : " + this.getClass().getSimpleName() + " mapPdbGenome snpeff.config genome pdbDir idMapFile");
		System.err.println("Command 'pdbdist'        : " + this.getClass().getSimpleName() + " pdbdist distanceThreshold aaMinSeparation path/to/pdb/dir id_map.txt");
		System.err.println("Command 'qhat'           : " + this.getClass().getSimpleName() + " qhat phylo.nh multiple_sequence_alignment.fa transition_matrix.txt");
		System.exit(-1);
	}

}
