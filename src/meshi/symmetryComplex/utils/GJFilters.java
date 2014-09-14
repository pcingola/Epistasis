package meshi.symmetryComplex.utils;

import meshi.util.filters.Filter;
import meshi.util.MeshiException;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.symmetryComplex.topologyMap.BoundariesMap;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:48:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class GJFilters {

    /**
         * Accepts residues and atoms whose residue numbers are within the ranges
         * given in the constructor. Passing valid arguments to the constructors is
         * not checked and is within the user's responsibility.
         *
         * throws ClassCastException
         *             if an object of a different type is checked for acception.
         * @author Oren Wolfshtat, Tatiana Maximova
         */
    public static class ResidueRangeFilter implements Filter {
            public final int[][] ranges;
        /*	public ResidueRangeFilter(final int... flatRanges) {
                ranges = new int[flatRanges.length / 2][];
                for (int i = 0; i < ranges.length; i++)
                    ranges[i] = new int[] { flatRanges[2*i], flatRanges[2*i + 1] };
            }
            public ResidueRangeFilter(final int from, final int to) {
                this(new int[] { from, to });
            }

            public ResidueRangeFilter(final int[]... ranges) {
                this.ranges = ranges;
            }
        */
            public ResidueRangeFilter(final int[][] ranges) {
                this.ranges = ranges;
            }

            public ResidueRangeFilter(final int[] flatRanges) {
                ranges = new int[flatRanges.length / 2][];
                for (int i = 0; i < ranges.length; i++)
                    ranges[i] = new int[] { flatRanges[2*i], flatRanges[2*i + 1] };
            }

            public ResidueRangeFilter(final int from, final int to) {
                this(new int[] { from, to });
            }

            public boolean accept(final Object obj) {
                final int residueNumber = (obj instanceof Residue) ? ((Residue) obj)
                        .number() : ((Atom) obj).residueNumber();
                for (final int[] range : ranges)
                    if (residueNumber >= range[0] & residueNumber <= range[1])
                        return true;
                return false;
            }
        }

/*
public static PartialCasFilter partialCasFilter = new PartialCasFilter();

public static PartialBackboneFilter partialBackboneFilter = new PartialBackboneFilter();

public static class PartialCasFilter implements Filter {
public boolean accept(final Object obj) {
final Atom atom = (Atom) obj;
if (!atom.name().equals("CA"))
return false;
final int residueNumber = atom.residueNumber();
for (final int[] tm : tmResiduesRanges)
if (tm[0] <= residueNumber & residueNumber <= tm[1])
return true;
return false;
}

} // class PartialCasFilter
public static class PartialBackboneFilter implements Filter {
public boolean accept(final Object obj) {
if (!(obj instanceof Atom))
return false;
final Atom atom = (Atom) obj;
final String name = atom.name();
if (!name.equals("N") & !name.equals("H") & !name.equals("CA")
& !name.equals("CB") & !name.equals("C")
& !name.equals("O"))
return false;
final int residueNumber = atom.residueNumber();
for (final int[] tm : tmResiduesRanges)
if (tm[0] <= residueNumber & residueNumber <= tm[1])
return true;
return false;
}

} // class PartialBackboneFilter
*/
        /*
              public static void freezeTMCas(final Protein protein) {
                  protein.atoms().filter(partialCasFilter).freeze();
              }
          /*
              private static void setTMHelicesSS(final Protein protein) {
                  final AtomList tmCas = protein.atoms().filter(partialCasFilter);
                  final Iterator tmCasIter = tmCas.iterator();
                  Atom tmCa;
                  while (tmCasIter.hasNext()) {
                      tmCa = (Atom) tmCasIter.next();
                    tmCa.residue().setSecondaryStructure('H');
                }
            }
        */

        public static class AndFilter implements Filter {
            private Filter[] filters;
            public AndFilter(Filter... filters) {
                this.filters = filters;
            }
            public boolean accept(Object obj) {
                for (Filter filter : filters)
                    if (! filter.accept(obj))
                        return false;
                return true;
            }
        }

        /**
         * Accepts atoms whose chain is identical to the one given to the constructor.
         */
        public static class AcceptChainFilter implements Filter {
            public final String chain;

            public AcceptChainFilter(String chain) {
                if (chain == null || chain.length() != 1)
                    throw new MeshiException
                        ("ChainFilter must receive a one-letter String");

                this.chain = chain;
            }

            public boolean accept(Object obj) {
                return (obj instanceof Atom) &&
                        chain.equals(((Atom) obj).chain());
            }
        }

        /**
         * Accepts only atoms whose chain is different than the one given to the constructor.
         */
        public static class RejectChainFilter implements Filter {
            public final String chain;

            public RejectChainFilter(String chain) {
                if (chain == null || chain.length() != 1)
                    throw new MeshiException
                        ("ChainFilter must receive a one-letter String");

                this.chain = chain;
            }

            public boolean accept(Object obj) {
                return (obj instanceof Atom) &&
                        ! chain.equals(((Atom) obj).chain());
            }
        }

    public static final Filter chainA_Ca_loops_Atom_Filter = new Filter() {
        private final Filter loopsFilter = new ResidueRangeFilter(BoundariesMap.loopMap());

        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            return (atom.name.equals("CA") && atom.chain().equals("A") && loopsFilter.accept(atom));
        }
    };
/*
    public static final Filter chainA_Ca_loops_ADD_Filter = new Filter() {
        private final Filter loopsFilter = new ResidueRangeFilter(BoundariesMap.loopMap());

        public boolean accept(Object obj) {
            AtomDensityData aDD = (AtomDensityData) obj;
            return (aDD.atomName.equals("CA") && aDD.chain.equals("A") && loopsFilter.accept(new DummyResidue(aDD.residueNumber)));
        }
    };
  */
    public static final Filter frozenResiduesFilter = new ResidueRangeFilter(BoundariesMap.frozenMap());
    public static final Filter absentResiduesFilter = new ResidueRangeFilter(BoundariesMap.absentMap());
    public static final Filter loopsResiduesFilter = new ResidueRangeFilter(BoundariesMap.loopMap());


}
