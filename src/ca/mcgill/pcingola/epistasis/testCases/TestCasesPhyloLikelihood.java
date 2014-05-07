package ca.mcgill.pcingola.epistasis.testCases;

import java.util.List;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeDna;
import ca.mcgill.pcingola.epistasis.phylotree.PhylogeneticTree;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix2Times;

/**
 * Phylogenetic trees
 *
 * @author pcingola
 */
public class TestCasesPhyloLikelihood extends TestCase {

	public void test_01() {
		// Example from: "Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 104

		// Leaf nodes
		LikelihoodTreeDna lt1 = new LikelihoodTreeDna("1");
		LikelihoodTreeDna lt2 = new LikelihoodTreeDna("2");
		LikelihoodTreeDna lt3 = new LikelihoodTreeDna("3");
		LikelihoodTreeDna lt4 = new LikelihoodTreeDna("4");
		LikelihoodTreeDna lt5 = new LikelihoodTreeDna("5");

		// Times
		double t1, t2, t3, t4, t5, t6, t7, t8;
		t1 = t2 = t3 = t4 = t5 = 0.2;
		t6 = t7 = t8 = 0.1;

		// Non-leaf nodes
		LikelihoodTreeDna lt7 = new LikelihoodTreeDna("7", lt1, t1, lt2, t2);
		LikelihoodTreeDna lt8 = new LikelihoodTreeDna("8", lt4, t4, lt5, t5);
		LikelihoodTreeDna lt6 = new LikelihoodTreeDna("6", lt7, t7, lt3, t3);
		LikelihoodTreeDna lt0 = new LikelihoodTreeDna("0", lt6, t6, lt8, t8);

		// Show tree
		System.out.println("Tree: " + lt0);

		List<PhylogeneticTree> leafNodes = lt0.child(true);
		List<PhylogeneticTree> nodes = lt0.child(false);
		System.out.println("Leaf nodes: " + leafNodes);

		// Set a sequence
		lt0.setLeafSequence("TCACC");

		// Matrix for time=0.1 (see page 103)
		// Note: In the book the bases are sorted "T C A G" instead of "A C G T"
		double v1 = 0.906563, v2 = 0.045855, v3 = 0.023791;
		double matrix1[][] = { //
				{ v1, v3, v2, v3 } // A
				, { v3, v1, v3, v2 } // C
				, { v2, v3, v1, v3 } // G
				, { v3, v2, v3, v1 } // T
		};

		// Matrix for time=0.2 (see page 103)
		// Note: In the book the bases are sorted "T C A G" instead of "A C G T"
		v1 = 0.825092;
		v2 = 0.084274;
		v3 = 0.045317;
		double matrix2[][] = { //
				{ v1, v3, v2, v3 } // A
				, { v3, v1, v3, v2 } // C
				, { v2, v3, v1, v3 } // G
				, { v3, v2, v3, v1 } // T
		};

		// Transition matrix provider
		TransitionMatrix2Times m2 = new TransitionMatrix2Times(matrix1, matrix2, 0.11);
		System.out.println(m2);

		// Calculate likelihood
		double pi[] = { 0.25, 0.25, 0.25, 0.25 };
		double likelihood = lt0.likelihood(m2, pi);

		// Show probabilities
		for (PhylogeneticTree lt : nodes)
			System.out.println(((LikelihoodTreeDna) lt).toStringP());

		// Check that likelihoods match the book's example
		// Note: In the book the bases are sorted "T C A G" instead of "A C G T"
		Assert.assertEquals(0.000112, lt0.getP(3), 1e-6);
		Assert.assertEquals(0.001838, lt0.getP(1), 1e-6);
		Assert.assertEquals(0.000075, lt0.getP(0), 1e-6);
		Assert.assertEquals(0.000014, lt0.getP(2), 1e-6);

		// Check likelihood
		Assert.assertEquals(0.00050975, likelihood, 1e-6);
	}
}
