package ca.mcgill.pcingola.epistasis.likelihood;

import ca.mcgill.mcb.pcingola.interval.Marker;

/**
 * Likelihood between two regions
 */
public class MarkerPairLikelihood {

	Marker marker1, marker2;
	double logLikelihoodRatio;

	public MarkerPairLikelihood(Marker marker1, Marker marker2, double logLikelihoodRatio) {
		// Store markers sorted by genomic region
		if (marker1.compareTo(marker2) <= 0) {
			this.marker1 = marker1;
			this.marker2 = marker2;
		} else {
			this.marker1 = marker2;
			this.marker2 = marker1;
		}

		this.logLikelihoodRatio = logLikelihoodRatio;
	}

	public double getLogLikelihoodRatio() {
		return logLikelihoodRatio;
	}

	public Marker getMarker1() {
		return marker1;
	}

	public Marker getMarker2() {
		return marker2;
	}

	@Override
	public String toString() {
		return marker1.getChromosomeName() + ":" + marker1.getStart() + "-" + marker1.getEnd() //
				+ "\t" + marker2.getChromosomeName() + ":" + marker2.getStart() + "-" + marker2.getEnd() //
				+ "\t" + logLikelihoodRatio;
	}
}
