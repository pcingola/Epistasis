package ca.mcgill.pcingola.epistasis.likelihood;

import java.util.ArrayList;

/**
 * Parameters distributed as mixtures of other distributions
 *
 * @author pcingola
 */
public class MixOfGaussians extends ParameterDistribution {

	ArrayList<Double> lambdas;
	ArrayList<ParameterDistribution> disributions;

	public MixOfGaussians() {
		lambdas = new ArrayList<Double>();
		disributions = new ArrayList<ParameterDistribution>();
	}

	public void add(double lambda, ParameterDistribution distr) {
		lambdas.add(lambda);
		disributions.add(distr);
	}

	@Override
	public double p(double x) {
		double p = 0;

		for (int i = 0; i < lambdas.size(); i++)
			p += lambdas.get(i) * disributions.get(i).p(x);

		return p;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < lambdas.size(); i++) {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(lambdas.get(i) + " * " + disributions.get(i));
		}

		return sb.toString();
	}
}
