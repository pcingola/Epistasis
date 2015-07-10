package ca.mcgill.pcingola.epistasis.likelihood;

import org.apache.commons.math3.analysis.function.Gaussian;

/**
 * Parameters distributed as Gaussians
 *
 * @author pcingola
 */
public class Gauss extends ParameterDistribution {

	double mu, sigma;
	Gaussian gaussian;

	public Gauss(double mu, double sigma) {
		this.mu = mu;
		this.sigma = sigma;
		gaussian = new Gaussian(mu, sigma);
	}

	@Override
	public double p(double x) {
		return gaussian.value(x);
		// return Stat.gaussian(mu, sigma, x);
	}

	@Override
	public String toString() {
		return "Gauss(" + mu + ", " + sigma + ")";
	}
}
