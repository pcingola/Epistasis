package ca.mcgill.pcingola.epistasis.likelihood;

import java.util.ArrayList;
import java.util.Iterator;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Model logistic regression parameter's distributions
 *
 * @author pcingola
 */
public class ParameterDistributionModel implements Iterable<ParameterDistribution> {

	String model;
	ArrayList<ParameterDistribution> distributions;

	public ParameterDistributionModel(String fileName, String model) {
		this.model = model;
		parse(fileName, model);
	}

	@Override
	public Iterator<ParameterDistribution> iterator() {
		return distributions.iterator();
	}

	/**
	 * Calculate probability
	 */
	public double p(double theta[]) {
		double p = 1.0;

		for (int i = 0; i < distributions.size(); i++) {
			ParameterDistribution pd = distributions.get(i);
			p *= pd.p(theta[i]);
		}

		return p;
	}

	/**
	 * Parse file
	 */
	void parse(String fileName, String model) {
		distributions = new ArrayList<ParameterDistribution>();
		String lines[] = Gpr.readFile(fileName).split("\n");

		for (String line : lines) {
			// Ignore comment lines
			if (line.startsWith("#")) continue;

			String fields[] = line.split("\t");

			int num = 0;
			ParameterDistribution distribution = null;

			// Parse 'model' lines
			if (fields[num++].equalsIgnoreCase(model)) {
				Gpr.parseIntSafe(fields[num++]);
				String distributionType = fields[num++];

				if (distributionType.equalsIgnoreCase("GAUSS")) {
					// Parse gaussian parameters
					double mu = Gpr.parseDoubleSafe(fields[num++]);
					double sigma = Gpr.parseDoubleSafe(fields[num++]);
					Gauss gauss = new Gauss(mu, sigma);
					distribution = gauss;
				} else if (distributionType.equalsIgnoreCase("MIXED")) {
					// Parse gaussian mixture parameters
					String lambdas[] = fields[num++].split(",");
					String mus[] = fields[num++].split(",");
					String sigmas[] = fields[num++].split(",");

					MixOfGaussians mix = new MixOfGaussians();
					for (int i = 0; i < lambdas.length; i++) {
						double lambda = Gpr.parseDoubleSafe(lambdas[i]);
						double mu = Gpr.parseDoubleSafe(mus[i]);
						double sigma = Gpr.parseDoubleSafe(sigmas[i]);

						Gauss gauss = new Gauss(mu, sigma);
						mix.add(lambda, gauss);
					}
					distribution = mix;

				} else throw new RuntimeException("Unknown distribution type '" + distributionType + "'");

				distributions.add(distribution);
			}
		}
	}

	@Override
	public String toString() {
		return toString(null);
	}

	public String toString(double theta[]) {
		StringBuilder sb = new StringBuilder();

		int paramNum = 0;
		sb.append("Distributions '" + model + "'" + (theta != null ? ": " + p(theta) : "") + "\n");
		for (ParameterDistribution pd : this) {
			sb.append("\t" + paramNum //
					+ (theta != null ? "\t" + theta[paramNum] + "\t" + pd.p(theta[paramNum]) : "") //
					+ "\t" + pd //
					+ "\n");
			paramNum++;
		}

		return sb.toString();
	}
}
