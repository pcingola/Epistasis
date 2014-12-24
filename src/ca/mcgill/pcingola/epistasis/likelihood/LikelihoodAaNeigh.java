package ca.mcgill.pcingola.epistasis.likelihood;

import ca.mcgill.pcingola.epistasis.gwas.GwasResult;

/**
 * Results form neighbouring LL(MSA) calculations
 *
 * @author pcingola
 */
public class LikelihoodAaNeigh {

	int count;
	double logLikelihoodRatioMiddle;
	double logLikelihoodRatio;
	double logLikelihoodAlt;
	double logLikelihoodNull;
	String direction;

	public LikelihoodAaNeigh() {
	}

	public LikelihoodAaNeigh(String direction, GwasResult gwasResult) {
		this.direction = direction;
		update(gwasResult);
		logLikelihoodRatioMiddle = gwasResult.logLikelihoodRatioMsa;
	}

	public double getAvgLogLikAlt() {
		if (count == 0) return 0.0;
		return logLikelihoodAlt / count;
	}

	public double getAvgLogLikNull() {
		if (count == 0) return 0.0;
		return logLikelihoodNull / count;
	}

	public double getAvgLogLikRatio() {
		if (count == 0) return 0.0;
		return logLikelihoodRatio / count;
	}

	@Override
	public String toString() {
		return direction //
				+ "\t" + logLikelihoodRatioMiddle //
				+ "\t" + logLikelihoodRatio //
				+ "\t" + getAvgLogLikRatio() //
				+ "\t" + getAvgLogLikNull() //
				+ "\t" + getAvgLogLikAlt() //
				+ "\t" + count //
		;
	}

	public void update(GwasResult gwasRes) {
		logLikelihoodRatio += gwasRes.logLikelihoodRatioMsa;
		logLikelihoodNull += Math.log(gwasRes.likelihoodMsaNull);
		logLikelihoodAlt += Math.log(gwasRes.likelihoodMsaAlt);
		count++;
	}

}
