package ca.mcgill.pcingola.epistasis.msa;

import java.util.List;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.entropy.EntropySeq;
import ca.mcgill.pcingola.epistasis.entropy.EntropySeq.InformationFunction;
import ca.mcgill.pcingola.epistasis.pdb.DistanceResult;

/**
 * Implement a 'similarity' by using N bases around the target
 *
 * @author pcingola
 */
public class MsaSimilarityN extends MsaSimilarity {

	public MsaSimilarityN(MultipleSequenceAlignmentSet msas, int numBases, InformationFunction func) {
		super(msas, func);
		this.numBases = numBases;
		double n = GprSeq.AMINO_ACIDS.length;
		double p = Math.pow(1.0 / n, 2 * numBases + 1);
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
			List<MultipleSequenceAlignment> msasTr = msas.getMsasByTrId(trId);
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
	public double calc(DistanceResult d) {
		return calc(d.msa1, d.msa2, d.msaIdx1, d.msaIdx2);
	}

	@Override
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		throw new RuntimeException("Do not use this calc method!");
	}

	double calc(String msaiId, String msajId, int posi, int posj) {
		// Find column sequences
		String seqsi[] = msas.colSequences(msaiId, posi, numBases);
		String seqsj[] = msas.colSequences(msajId, posj, numBases);

		if (debug) {
			System.out.println(msaiId + "\t" + posi);
			for (int i = 0; i < seqsi.length; i++)
				System.out.println("\t" + (seqsi[i] != null ? seqsi[i] : ""));

			System.out.println(msajId + "\t" + posj);
			for (int i = 0; i < seqsi.length; i++)
				System.out.println("\t" + (seqsj[i] != null ? seqsj[i] : ""));
		}

		// Calculate
		double score = 0;
		EntropySeq entropy = new EntropySeq();
		entropy.calc(seqsi, seqsj);
		switch (func) {
		case MI:
			score = entropy.getMi();
			break;

		case VARINF:
			score = entropy.getVarInf();
			break;

		case HXY:
			score = entropy.getHxy();
			break;

		default:
			throw new RuntimeException("Unknown function type " + func);
		}

		incScore(score);
		return score;

	}

}
