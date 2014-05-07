package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.phylotree.PhylogeneticTree;

/**
 * Phylogenetic trees
 *
 * @author pcingola
 */
public class TestCasesPhylo extends TestCase {

	public void test_01() {
		String phyloFile = "test/hg19.100way.commonNames.nh";
		String phyloTxtFile = "test/hg19.100way.commonNames.txt";

		PhylogeneticTree tree = new PhylogeneticTree();
		tree.load(phyloFile);

		String treeExpected = Gpr.readFile(phyloTxtFile).replace('\n', ' ').trim();

		Assert.assertEquals(treeExpected, tree.toString());
	}

	public void test_02() {
		String phyloFile = "test/hg19.100way.commonNames.nh";

		PhylogeneticTree tree = new PhylogeneticTree();
		tree.load(phyloFile);

		double d = tree.distance("Human", "Chimp");
		Assert.assertEquals(0.006429 + 0.00638, d, 1e-6);
	}

	public void test_03() {
		String phyloFile = "test/hg19.100way.commonNames.nh";

		PhylogeneticTree tree = new PhylogeneticTree();
		tree.load(phyloFile);

		double d = tree.distance("Human", "Gorilla");
		Assert.assertEquals(0.006429 + 0.002176 + 0.008821, d, 1e-6);
	}

}
