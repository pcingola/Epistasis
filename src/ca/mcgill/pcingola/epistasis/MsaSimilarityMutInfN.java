package ca.mcgill.pcingola.epistasis;

import java.util.List;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.entropy.Entropy;

/**
 * Implement a 'similarity' by mutual information using N bases around the target
 *
 * @author pcingola
 */
public class MsaSimilarityMutInfN extends MsaSimilarity {

	public MsaSimilarityMutInfN(MultipleSequenceAlignmentSet msas, int numBases) {
		super(msas);
		this.numBases = numBases;
		double n = GprSeq.AMINO_ACIDS.length;
		double p = Math.pow(1.0 / n, numBases);
		maxScore = -Math.log(p) / Math.log(2.0); // Maximum possible entropy
	}

	/**
	 * Calculate one sample of a random distribution
	 */
	@Override
	void backgroundDistribution() {
		while (true) {
			// Select one MSA and position randomly
			MultipleSequenceAlignment msai = msas.rand(random);
			int posi = msai.randomColumnNumber(random);
			if (msai.isSkip(posi)) continue;

			// Select another MSA and position randomly
			String trId = msai.getTranscriptId();
			List<MultipleSequenceAlignment> msasTr = msas.getMsas(trId);
			MultipleSequenceAlignment msaj = msasTr.get(random.nextInt(msasTr.size()));
			int posj = msaj.randomColumnNumber(random);
			if (msaj.isSkip(posj)) continue;

			// Same MSA and position? Find another random
			if (posi == posj && msai.getId().equals(msaj.getId())) continue;

			// Calculate
			double calc = calc(msai.getId(), msaj.getId(), posi, posj);
			if (debug) System.err.println(calc + "\t" + showSeqs(msai, msaj, posi, posj));
			else if (verbose) System.out.printf("%.6e\n", calc);
			else Gpr.showMark(count++, SHOW_EVERY);
			return;
		}
	}

	@Override
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		throw new RuntimeException("Do not use this calc method!");
	}

	public double calc(String msaiId, String msajId, int posi, int posj) {
		// Find column sequences
		String seqsi[] = msas.findColSequences(msaiId, posi, numBases);
		String seqsj[] = msas.findColSequences(msajId, posj, numBases);

		System.out.println(msaiId + "\t" + posi);
		for (int i = 0; i < seqsi.length; i++)
			System.out.println("\t" + (seqsi[i] != null ? seqsi[i] : ""));

		System.out.println(msajId + "\t" + posj);
		for (int i = 0; i < seqsi.length; i++)
			System.out.println("\t" + (seqsj[i] != null ? seqsj[i] : ""));

		Entropy entropy = new Entropy();
		entropy.calc(seqsi, seqsj);
		return entropy.getMi();
	}

	//	/**
	//	 * Measure similarity: Correlation between two loci
	//	 */
	//	@Override
	//	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
	//		// Sanity check
	//		if ((posi + numBases) >= msai.getSeqLen()) return Double.NaN;
	//		if ((posj + numBases) >= msaj.getSeqLen()) return Double.NaN;
	//
	//		int numAligns = msas.getNumAligns();
	//
	//		//---
	//		// Initialize counters
	//		//---
	//		Entropy entropy = new Entropy(numBases);
	//
	//		//---
	//		// Count matching bases
	//		//---
	//		for (int i = 0; i < numAligns; i++) {
	//			long basesI = msai.getCodeLong(i, posi, numBases);
	//			if (basesI == ALIGN_GAP) continue;
	//
	//			long basesJ = msaj.getCodeInt(i, posj, numBases);
	//			if (basesJ == ALIGN_GAP) continue;
	//
	//			// Count
	//			entropy.inc(basesI, basesJ);
	//		}
	//
	//		// Only a few organisms align? Filter out
	//		if (entropy.getCount() < minCount) return Double.NaN;
	//
	//		double mutInf = entropy.mi();
	//
	//		// Results
	//		incScore(mutInf);
	//		if (mutInf > threshold) return mutInf;
	//		return Double.NaN;
	//	}

}
