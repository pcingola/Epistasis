package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Read and parse phylogenetic tree files (*.nh)
 *
 * @author pcingola
 */
public class PhylogeneticTree {

	public static boolean debug = false;

	public static final int NO_UNIFORM_CODE = -2;

	int sequenceCode = -1;
	int uniformCode = NO_UNIFORM_CODE; // Do all leaf nodes have the same code?
	int level = -1;
	String name; // Leaf nodes only have 'name'
	int id; // All nodes have 'id'
	PhylogeneticTree parent, left, right;
	List<PhylogeneticTree> leafNodes;
	UniformTreeValueCache lcache;
	double distanceLeft, distanceRight;

	/**
	 * Create root node
	 */
	public PhylogeneticTree() {
		parent = left = right = null;
		distanceLeft = distanceRight = 0;
		sequenceCode = -1;
	}

	/**
	 * Create and parse descendant nodes
	 */
	public PhylogeneticTree(PhylogeneticTree parent, String phyloStr, Counter ids) {
		this.parent = parent;
		left = right = null;
		distanceLeft = distanceRight = 0;
		sequenceCode = -1;
		parse(phyloStr, ids);
	}

	public PhylogeneticTree(String name) {
		parent = left = right = null;
		this.name = name;
		distanceLeft = distanceRight = 0;
		sequenceCode = -1;
	}

	public PhylogeneticTree(String name, PhylogeneticTree left, double distanceLeft, PhylogeneticTree right, double distanceRight) {
		this.left = left;
		this.right = right;
		this.name = name;
		this.distanceLeft = distanceLeft;
		this.distanceRight = distanceRight;
		sequenceCode = -1;
	}

	/**
	 * Get a balanced parenthesis string
	 */
	int balancedIdx(String phylo, int start) {
		if (phylo.charAt(start) != '(') return phylo.indexOf(':', start) - 1; // Must be a leaf node, return index right before ':'

		int countParen = 1;
		for (int i = start + 1; i < phylo.length(); i++) {
			char c = phylo.charAt(i);

			if (c == '(') countParen++;
			else if (c == ')') countParen--;

			if (countParen == 0) return i;
		}

		throw new RuntimeException("String does not balance:\n\t" + phylo);
	}

	/**
	 * Names for all child nodes
	 */
	public List<PhylogeneticTree> child(boolean onlyLeaf) {
		ArrayList<PhylogeneticTree> nodes = new ArrayList<PhylogeneticTree>();
		if (left != null) nodes.addAll(left.child(onlyLeaf));
		if (!onlyLeaf || isLeaf()) nodes.add(this);
		if (right != null) nodes.addAll(right.child(onlyLeaf));
		return nodes;
	}

	/**
	 * Names for all child nodes
	 */
	public List<String> childNames() {
		ArrayList<String> names = new ArrayList<String>();
		if (left != null) names.addAll(left.childNames());
		if (name != null) names.add(name);
		if (right != null) names.addAll(right.childNames());
		return names;
	}

	/**
	 * Find distance to node recursing down
	 */
	public double distance(String name) {
		if (name.equals(this.name)) return 0;

		double d = -1;
		if (left != null) d = left.distance(name);
		if (d >= 0) return d + distanceLeft;

		if (right != null) d = right.distance(name);
		if (d >= 0) return d + distanceRight;

		return -1;
	}

	/**
	 * Find distance to node
	 */
	public double distance(String name1, String name2) {
		if (name1.equals(name2)) return 0;

		double dl1 = -1, dl2 = -1;
		if (left != null) {
			dl1 = left.distance(name1);
			dl2 = left.distance(name2);
		}

		// Both nodes on left side: recurse
		if ((dl1 >= 0) && (dl2 >= 0)) return left.distance(name1, name2);

		double dr1 = -1, dr2 = -1;
		if (right != null) {
			dr1 = right.distance(name1);
			dr2 = right.distance(name2);
		}

		// Both nodes on right side: recurse
		if ((dr1 >= 0) && (dr2 >= 0)) return right.distance(name1, name2);

		// Cannot find node
		if (dl1 < 0 && dr1 < 0) return -1;
		if (dl2 < 0 && dr2 < 0) return -1;

		// One node to the left, one node to the right: Add distances
		return (dl1 >= 0 ? dl1 + distanceLeft : dr1 + distanceRight) + (dl2 >= 0 ? dl2 + distanceLeft : dr2 + distanceRight);
	}

	/**
	 * Find node by name
	 */
	public PhylogeneticTree find(String name) {
		if (name.equals(this.name)) return this;

		PhylogeneticTree node = null;
		if (left != null) node = left.find(name);
		if (node != null) return left;
		if (right != null) return right.find(name);
		return null;
	}

	public double getDistanceLeft() {
		return distanceLeft;
	}

	public double getDistanceRight() {
		return distanceRight;
	}

	public int getId() {
		return id;
	}

	public UniformTreeValueCache getLcache() {
		return lcache;
	}

	public PhylogeneticTree getLeft() {
		return left;
	}

	public int getLevel() {
		if (level < 0) {
			if (parent != null) level = parent.getLevel() + 1;
			else level = 0;
		}
		return level;
	}

	public String getName() {
		return name;
	}

	public PhylogeneticTree getParent() {
		return parent;
	}

	public PhylogeneticTree getRight() {
		return right;
	}

	public char getSequence() {
		if (sequenceCode < 0) return ' ';
		return GprSeq.code2aa((byte) sequenceCode);
	}

	public String getSequenceAaPair() {
		if (sequenceCode < 0) return "  ";
		return GprSeq.code2aaPair(sequenceCode);
	}

	public int getUniformCode() {
		return uniformCode;
	}

	/**
	 * Is this a GAP? (a leaf node assigned a GAP as sequence)
	 */
	public boolean isGap() {
		return isLeaf() && sequenceCode < 0;
	}

	/**
	 * Is this a leaf node?
	 */
	public boolean isLeaf() {
		return left == null && right == null;
	}

	/**
	 * Load from file
	 */
	public void load(String phyloFile) {
		String phylo = Gpr.readFile(phyloFile);
		phylo = Gpr.noSpaces(phylo.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ').replace(';', ' ')); // Remove unwanted characters
		parse(phylo, new Counter());
	}

	protected PhylogeneticTree newNode(PhylogeneticTree parent, String phyloStr, Counter ids) {
		return new PhylogeneticTree(parent, phyloStr, ids);
	}

	/**
	 * Get a balanced parenthesis string
	 */
	int numIdx(String phylo, int start) {
		if (phylo.charAt(start) != ':') return -1;

		for (int i = start + 1; i < phylo.length(); i++) {
			char c = phylo.charAt(i);
			if (!Character.isDigit(c) && c != '.') return i;
		}

		return phylo.length();
	}

	/**
	 * Parse a phylogenetic string (from an 'NH' file)
	 */
	void parse(String phylo, Counter ids) {
		if (debug) System.out.println(phylo);

		if (phylo.charAt(0) == '(' && phylo.charAt(phylo.length() - 1) == ')') phylo = phylo.substring(1, phylo.length() - 1); // Remove first and last parenthesis

		// Parse left node
		int leftIdx = balancedIdx(phylo, 0);
		if (leftIdx <= 0) {
			name = phylo;
			return;
		}

		// Parse distance
		int leftNum = numIdx(phylo, leftIdx + 1);
		if (leftNum <= 0) return;

		// Comma
		if (phylo.charAt(leftNum) != ',') return;

		// Parse right node
		int rightIdx = balancedIdx(phylo, leftNum + 1);
		if (rightIdx <= 0) return;

		// Parse distance
		int rightNum = numIdx(phylo, rightIdx + 1);
		if (rightNum <= 0) return;

		String leftStr = phylo.substring(0, leftIdx + 1);
		String leftDistStr = phylo.substring(leftIdx + 2, leftNum);
		distanceLeft = Gpr.parseDoubleSafe(leftDistStr);

		String rightStr = phylo.substring(leftNum + 1, rightIdx + 1);
		String rightDistStr = phylo.substring(rightIdx + 2, rightNum);
		distanceRight = Gpr.parseDoubleSafe(rightDistStr);

		id = (int) ids.inc();
		left = newNode(this, leftStr, ids);
		right = newNode(this, rightStr, ids);

	}

	protected void resetNode(int size) {
	}

	public void setDistanceLeft(double distanceLeft) {
		this.distanceLeft = distanceLeft;
	}

	public void setDistanceRight(double distanceRight) {
		this.distanceRight = distanceRight;
	}

	public void setLcache(UniformTreeValueCache lcache) {
		this.lcache = lcache;
		if (left != null) left.setLcache(lcache);
		if (right != null) right.setLcache(lcache);
	}

	/**
	 * Set a string as leaf node sequences
	 */
	public void setLeafSequence(String sequence) {
		// Get all leafs and set each one
		if (leafNodes == null) leafNodes = child(true);

		// Sanity check
		if (leafNodes.size() != sequence.length()) throw new RuntimeException("Incompatible lengths:\n\tTree leaf nodes: " + leafNodes.size() + "\n\t" + sequence.length());

		// Set sequence
		for (int i = 0; i < sequence.length(); i++)
			leafNodes.get(i).setSequence(sequence.charAt(i));
	}

	public void setLeafSequenceAaPair(byte seq1[], byte seq2[]) {
		if (seq1.length != seq2.length) throw new RuntimeException("Sequence lengths do not match: " + seq1.length + " != " + seq2.length);

		// Get all leafs and set each one
		if (leafNodes == null) leafNodes = child(true);

		// Sanity check
		if (leafNodes.size() != seq1.length) throw new RuntimeException("Incompatible lengths:\n\tTree leaf nodes: " + leafNodes.size() + "\n\t" + seq1.length);

		// Set sequence
		for (int i = 0; i < seq1.length; i++)
			leafNodes.get(i).setSequence(seq1[i], seq2[i]);

	}

	public void setLeafSequenceAaPair(String seq1, String seq2) {
		if (seq1.length() != seq2.length()) throw new RuntimeException("Sequence lengths do not match: " + seq1.length() + " != " + seq2.length());

		// Get all leafs and set each one
		if (leafNodes == null) leafNodes = child(true);

		// Sanity check
		if (leafNodes.size() != seq1.length()) throw new RuntimeException("Incompatible lengths:\n\tTree leaf nodes: " + leafNodes.size() + "\n\t" + seq1.length());

		// Set sequence
		for (int i = 0; i < seq1.length(); i++)
			leafNodes.get(i).setSequence(seq1.charAt(i), seq2.charAt(i));
	}

	/**
	 * Set a string as leaf node sequences
	 */
	public void setLeafSequenceCode(int sequenceCode[]) {
		// Get all leafs and set each one
		if (leafNodes == null) leafNodes = child(true);

		// Sanity check
		if (leafNodes.size() != sequenceCode.length) throw new RuntimeException("Incompatible lengths:\n\tTree leaf nodes: " + leafNodes.size() + "\n\t" + sequenceCode.length);

		// Set sequence
		for (int i = 0; i < sequenceCode.length; i++)
			leafNodes.get(i).setSequenceCode(sequenceCode[i]);
	}

	public void setLeft(PhylogeneticTree left) {
		this.left = left;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setRight(PhylogeneticTree right) {
		this.right = right;
	}

	public void setSequence(byte aa1, byte aa2) {
		sequenceCode = GprSeq.aaPairCode(aa1, aa2);
	}

	public void setSequence(char aa) {
		sequenceCode = GprSeq.aa2Code(aa);
	}

	public void setSequence(char aa1, char aa2) {
		sequenceCode = GprSeq.aaPairCode(aa1, aa2);
	}

	public void setSequenceCode(int seqCode) {
		sequenceCode = seqCode;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		if (left != null && right != null) {
			sb.append("(");
			if (left != null) sb.append(left + ":" + distanceLeft);
			sb.append(",");
			if (right != null) sb.append(right + ":" + distanceRight);
			sb.append(")");
		} else {
			return name;
		}

		return sb.toString();
	}

	/**
	 * Do all leaf nodes have the same code?
	 */
	public int uniformCode() {
		if (isLeaf()) {
			uniformCode = sequenceCode;
			return sequenceCode;
		}

		int codeLeft = left.uniformCode();
		int codeRight = right.uniformCode();

		if (codeLeft != codeRight //
				|| codeLeft == NO_UNIFORM_CODE //
				|| codeRight == NO_UNIFORM_CODE //
				) uniformCode = NO_UNIFORM_CODE;
		else uniformCode = codeLeft;

		return uniformCode;
	}
}
