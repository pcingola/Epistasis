package ca.mcgill.pcingola.epistasis;

/**
 * For results with 3 parameters
 */
public class Triplet<T1, T2, T3> {

	public T1 a;
	public T2 b;
	public T3 c;

	public Triplet(T1 a, T2 b, T3 c) {
		this.a = a;
		this.b = b;
		this.c = c;
	}
}
