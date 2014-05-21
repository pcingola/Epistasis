package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.snpEffect.commandLine.CommandLine;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTree;
import ca.mcgill.pcingola.epistasis.phylotree.MaxLikelihoodTm;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrixMarkov;

/**
 * Main command line
 *
 * @author pcingola
 */
public class Epistasis implements CommandLine {

	public static int MIN_DISTANCE = 1000000;

	public static void main(String[] args) {
		Epistasis epistasis = new Epistasis(args);
		epistasis.run();
	}

	String[] args;

	public Epistasis(String[] args) {
		this.args = args;
	}

	@Override
	public String[] getArgs() {
		return args;
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length < 1) usage("Missing command");
		String cmd = args[0];

		switch (cmd) {
		case "corr":
			int numAligns = Gpr.parseIntSafe(args[1]);
			String multAlignFile = args[2];
			if (numAligns <= 0) usage("number of alignments must be positive number");
			runMsaCorr(numAligns, multAlignFile);
			break;

		case "mi":
			numAligns = Gpr.parseIntSafe(args[1]);
			int numBases = Gpr.parseIntSafe(args[2]);
			multAlignFile = args[3];
			if (numAligns <= 0) usage("number of alignments must be positive number");
			runMsaMi(numAligns, numBases, multAlignFile);
			break;

		case "pdbdist":
			// Parse command line
			double distThreshold = Gpr.parseDoubleSafe(args[1]);
			if (distThreshold <= 0) usage("Distance must be a positive number: '" + args[1] + "'");
			String pdbDir = args[2];
			String idMapFile = args[3];
			String treeFile = args[4];
			multAlignFile = args[5];
			runPdbDist(distThreshold, pdbDir, idMapFile, treeFile, multAlignFile);
			break;

		case "phylo":
			if (args.length < 4) usage("Missing arguments for command '" + cmd + "'");
			treeFile = args[1];
			multAlignFile = args[2];
			String qMatrixFile = args[3];
			runPhylo(treeFile, multAlignFile, qMatrixFile);
			break;

		case "test":
			runTest(args);
			break;

		default:
			throw new RuntimeException("Unknown command: '" + cmd + "'");
		}
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
	 * Run correlation
	 * @param numAligns
	 * @param multAlign
	 */
	void runMsaCorr(int numAligns, String multAlign) {
		// Load
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(multAlign, numAligns);
		msas.load();

		// Run similarity
		MsaSimilarity sim = new MsaSimilarity(msas);
		sim.similarity();

		System.out.println("Score statistics:\n" + sim);
	}

	/**
	 * Run Mutual Information
	 *
	 * @param numAligns
	 * @param numBases
	 * @param multAlign
	 */
	void runMsaMi(int numAligns, int numBases, String multAlign) {
		// Load
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(multAlign, numAligns);
		msas.load();

		// Run similarity
		MsaSimilarity sim = numBases > 1 ? new MsaSimilarityMutInfN(msas, numBases) : new MsaSimilarityMutInf(msas);
		sim.similarity();

		System.out.println("Score statistics:\n" + sim);
	}

	/**
	 * Run Pdb distance
	 * @param distThreshold
	 * @param pdbDir
	 * @param idMapFile
	 * @param treeFile
	 * @param multAlignFile
	 */
	void runPdbDist(double distThreshold, String pdbDir, String idMapFile, String phyloFileName, String multAlignFile) {
		// Load: tree
		Timer.showStdErr("Loading phylogenetic tree from " + phyloFileName);
		LikelihoodTree tree = new LikelihoodTree();
		tree.load(phyloFileName);
		int numAligns = tree.childNames().size();

		// Load: MSA
		Timer.showStdErr("Loading " + numAligns + " way multiple alignment from " + multAlignFile);
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(multAlignFile, numAligns);
		msas.load();

		Timer.showStdErr("Loading id maps " + idMapFile);
		IdMapper idMapper = new IdMapper(idMapFile);

		// Run analysis
		PdbDistanceAnalysis pdban = new PdbDistanceAnalysis(pdbDir, distThreshold, idMapper);
		pdban.run();
		System.out.println(pdban);

		// Write results
		String outFile = "pdb_distance_by_AA_pos.txt";
		Gpr.toFile(outFile, pdban);
		System.err.println("Distance metrics file written to: " + outFile);
	}

	/**
	 * Run phylogenetic analysis
	 * @param phyloFileName
	 * @param numAligns
	 * @param multAlign
	 * @return
	 */
	public boolean runPhylo(String phyloFileName, String multAlign, String qMatrixFile) {
		//---
		// Initialize
		//---

		// Load: tree
		Timer.showStdErr("Loading phylogenetic tree from " + phyloFileName);
		LikelihoodTree tree = new LikelihoodTree();
		tree.load(phyloFileName);
		int numAligns = tree.childNames().size();

		// Load: MSA
		Timer.showStdErr("Loading " + numAligns + " way multiple alignment from " + multAlign);
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(multAlign, numAligns);
		msas.load();

		// Sanity check: Make sure that the alignment and the tree match
		StringBuilder sbMsa = new StringBuilder();
		for (String sp : msas.getSpecies())
			sbMsa.append(sp + " ");

		StringBuilder sbTree = new StringBuilder();
		for (String sp : msas.getSpecies())
			sbTree.append(sp + " ");

		if (!sbTree.toString().equals(sbMsa.toString())) throw new RuntimeException("Species form MSA and Tree do not match:\n\tMSA : " + sbMsa + "\n\tTree: " + sbTree);

		System.out.println("\nSpecies: " + sbTree);

		// Load: Q (calculate if not available)
		MaxLikelihoodTm mltm = new MaxLikelihoodTm(tree, msas);
		TransitionMatrix Q;

		// Load or calculate transition matrix
		if (Gpr.canRead(qMatrixFile)) {
			Timer.showStdErr("Loading Q matrix  from file '" + qMatrixFile);
			Q = mltm.loadTransitionMatrix(qMatrixFile);
		} else {
			Timer.showStdErr("Q matrix file '" + qMatrixFile + "' not found, calculating matrix");
			Q = mltm.estimateTransitionMatrix();
			Timer.showStdErr("Saving Q matrix to file '" + qMatrixFile + "'");
			Q.save(qMatrixFile);
		}
		System.out.println("Q matrix:\n" + Q.toString());
		mltm.showEienQ();

		// Calculate Bayes Factor

		// For each MSA...
		mltm.calcPi();
		for (MultipleSequenceAlignment msa1 : msas) {
			for (MultipleSequenceAlignment msa2 : msas) {
				// Make sure the MSAs are far away from each other
				if (msa2.compareTo(msa1) < MIN_DISTANCE) continue;

				// All pair of positions
				for (int pos1 = 0; pos1 < msa1.length(); pos1++) {
					if (msa1.isSkip(pos1)) continue;

					// Log likelihood at position 1
					double loklik1 = mltm.logLikelyhood(msa1, pos1);

					for (int pos2 = 0; pos2 < msa2.length(); pos2++) {
						if (msa2.isSkip(pos2)) continue;

						// Log likelihood at position 2
						double loklik2 = mltm.logLikelyhood(msa2, pos2);
						System.out.println(msa1.getId() + "\t" + msa2.getId() + "\t" + pos1 + "\t" + pos2 + "\t" + loklik1 + "\t" + loklik2);

						// Combined log likelihood
						mltm.logLikelyhood(msa1, pos1, msa2, pos2);
					}
				}
			}
		}

		return true;
	}

	public boolean runTest(String args[]) {
		int N = 400;
		double d[][] = new double[N][N];

		// Rand matrix
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				d[i][j] = Math.random();

		TransitionMatrixMarkov Q = new TransitionMatrixMarkov(d);
		Timer.show("Start Exp");
		Q.exp(1.0);
		Timer.show("End Exp");

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: " + this.getClass().getSimpleName() + " cmd options");
		System.err.println("Command 'corr'      : " + this.getClass().getSimpleName() + " corr number_of_aligns multiple_alignment_file.fa");
		System.err.println("Command 'mi'        : " + this.getClass().getSimpleName() + " mi number_of_bases number_of_aligns multiple_alignment_file.fa");
		System.err.println("Command 'pdbdist'   : " + this.getClass().getSimpleName() + " pdbdist distanceThreshold path/to/pdb/dir id_map.txt phylo.nh multiple_sequence_alignment.fa");
		System.err.println("Command 'phylo'     : " + this.getClass().getSimpleName() + " phylo phylo.nh multiple_sequence_alignment.fa transition_matrix.txt");
		System.err.println("Command 'test'      : " + this.getClass().getSimpleName() + " ...");
		System.exit(-1);
	}

}
