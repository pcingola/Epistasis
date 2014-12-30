package ca.mcgill.pcingola.epistasis.likelihood;

/**
 * Parameters distributed as Gaussians
 * 
 * @author pcingola
 */
public class Gauss extends ParameterDistribution {

	double mu, sigma;

	public Gauss(double mu, double sigma) {
		this.mu = mu;
		this.sigma = sigma;
	}

	@Override
	public String toString() {
		return "Gauss(" + mu + ", " + sigma + ")";
	}
}
