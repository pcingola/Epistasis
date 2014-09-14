/*
 * Created on 26/12/2004
 */
package meshi.energy.hydrogenBondsPairs;
import meshi.energy.Parameters;
import meshi.util.*;
import java.util.StringTokenizer;

/**
 * @author amilev
 *
 * Let di be the residue number of a hydrogen-donor at hydrogen-bond i,
 *  and ai the residue number of a hydrogen-acceptors in the same hydrogen-bond.
 *  Thus, we can uniquely represent each pattern <hbi ,hbk> by sextuplet:  <|ai-di| ,|ak-dk|,|di-dk|,  |di-ak|,|ai-dk|,|ai-ak|>.
 *  To reduce the number of possible sextuplets, |ai-di| is bounded by 6, and|ai-dk| is bounded by 10.

 **/
public class HydrogenBondsPairsParameters implements Comparable, Parameters  {

    //------------------------------------ data ------------------------------------
    /**
     * Key values:
     * firstHBSeqGap = the gap between the two atoms in the first HB   (O1-H1)
     * secondHBSeqGap = the gap between the two atoms in the second HB   (O2-H2)
     * h1h2SeqGap = the gap between the two Hydrogens
     * o1o2SeqGap = the gap between the two Oxygens
     * h1o2SeqGap = the gap between the H of the first HB and the O of the second HB
     * o1h2SeqGap = the gap between the O of the first HB and the H of the second HB
     * value is proportionaly to the friquency of such pair of HBs in the data base
     **/
	public final int firstHBSeqGap,secondHBSeqGap,h1h2SeqGap,o1o2SeqGap,h1o2SeqGap,o1h2SeqGap;

    private int value=0;
    public final int value(){return value;}
    
    //------------------------------------ constructors ----------------------------

    /*
     * create an element from line that contains all the key values and the friquncy value
     */
    public HydrogenBondsPairsParameters(String line) {
		this(new StringTokenizer(line,": "));
    }
	    
    private HydrogenBondsPairsParameters(StringTokenizer line) {
        line.nextToken();
        firstHBSeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        secondHBSeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        h1h2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        o1o2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        h1o2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        o1h2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        value = Utils.toInt(line.nextToken());
    }
	    
    public HydrogenBondsPairsParameters(int firstHBGap,int secondHBGap,int h1h2Gap,int o1o2Gap,int h1o2Gap,int o1h2Gap){
        firstHBSeqGap =firstHBGap;
        secondHBSeqGap = secondHBGap;
        h1h2SeqGap = h1h2Gap;
        o1o2SeqGap = o1o2Gap;
        h1o2SeqGap = h1o2Gap;
        o1h2SeqGap = o1h2Gap;
    }

    public HydrogenBondsPairsParameters(int firstHBGap,int secondHBGap,int h1h2Gap,int o1o2Gap,int h1o2Gap,int o1h2Gap,int value){
        this(firstHBGap,secondHBGap,h1h2Gap,o1o2Gap,h1o2Gap,o1h2Gap);
        this.value = value;
    }

    //--------------------------------- methods ------------------------------------------
    
    public String toString()
    {
        return "firstBond: "+firstHBSeqGap+"   secondBond: "+secondHBSeqGap+"    h_dis: "+h1h2SeqGap+"    o_dis: "+o1o2SeqGap+"    ho-dis: "+h1o2SeqGap+"    oh_dis: "+o1h2SeqGap+" value: "+value;
    }
	    
	   
    public int compareTo(Object obj) {
        HydrogenBondsPairsParameters other = (HydrogenBondsPairsParameters)obj;
		if (firstHBSeqGap == other.firstHBSeqGap) {
		    if (secondHBSeqGap == other.secondHBSeqGap) {
                if (h1h2SeqGap == other.h1h2SeqGap) {
                    if (o1o2SeqGap == other.o1o2SeqGap) {
                        if (h1o2SeqGap == other.h1o2SeqGap) {
                            if (o1h2SeqGap == other.o1h2SeqGap) return 0;
                            else {
                                if (o1h2SeqGap > other.o1h2SeqGap) return 1;
                                return -1;
                            }
                        }
                        else {
                            if (h1o2SeqGap > other.h1o2SeqGap) return 1;
                            return -1;
                        }
                    }
                    else {
                        if (o1o2SeqGap > other.o1o2SeqGap) return 1;
                        return -1;
                    }
                }
                else {
                    if (h1h2SeqGap > other.h1h2SeqGap) return 1;
                    return -1;
                }
		    }
		    else {
                if (secondHBSeqGap > other.secondHBSeqGap) return 1;
                return -1;
		    }
		}
		else {
		    if (firstHBSeqGap > other.firstHBSeqGap) return 1;
		    return -1;
		}
    }


    // 	    public int compareTo(Object other) {
    // 		//  A comment on comment - instanceof is a very heavy opperation so I commented out the 
    // 		// three lines below.
    // 		// 		if (! (other instanceof  HydrogenBondsPairsParameters))
    // 		// 		    throw new RuntimeException("Weird argument to "+
    // 		// 					       " HydrogenBondsPairsParameters.compairTo(Object other)");
    // 	        long l1,l2,c256;
    // 	        HydrogenBondsPairsParameters obj2;
    // 	        obj2 = (HydrogenBondsPairsParameters)other;
    // 	        c256 = 256000;
    // 	        c256 *= 1000000;
    // 	        l1= (val6+20)
    // 	            +(val5+20)*40
    // 	            +(val4+20)*1600
    // 	            +(val3+20)*64000
    // 	            +(val2+1000)*((long)128000000)
    // 	            +(val1+1000)*c256;
    // 	        l2= (obj2.val6+20)
    // 	            +(obj2.val5+20)*40
    // 	            +(obj2.val4+20)*1600
    // 	            +(obj2.val3+20)*64000
    // 	            +(obj2.val2+1000)*((long)128000000)
    // 	            +(obj2.val1+1000)*c256;
    // 	        if(l1 < l2)
    // 	            return -1;
    // 	        else if (l1 == l2)
    // 	            return 0;
    // 	        else return 1;
    // 	    }
	    
	    
}
