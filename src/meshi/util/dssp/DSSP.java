package meshi.util.dssp;
import meshi.molecularElements.Residue;
import meshi.sequences.*;
import meshi.parameters.*;
import java.io.BufferedReader;
import java.io.FileReader;

public class DSSP {
    protected String fileName;
    private char[] aa;
    private int[] resNum;
    private char[] ss;
    private double[] solvACC;
    private double[] fullSol = {115,135,150,190,210,
				75 ,195,175,200,170,
				185,160,145,180,225,
				115,140,155,255,230};
    private double hbPrecentage;
    private double parallelPrecentage;
    private double antyparallelPrecentage;

    public void printPresentage(){
	System.out.println("Total hb: "+hbPrecentage);
	System.out .println("parallel: "+parallelPrecentage );
	System.out .println("AntiParallel: "+antyparallelPrecentage );

	System.out.println("Number Of Betta Residues: ");
    }

    public DSSP(String dsspFileName)
    {
	this.fileName = dsspFileName;
	readDSSP();
    }
	
    /*
     *This function returns true if the DSSP assignment of residue number 'res' 
     *appears in the given char list. Else false is returned.
     */
    public boolean inList(int res, char[] SSlist) {
	char ch;
	int i;
		
	ch = SSofRes(res);
	for (i=0 ; i<SSlist.length ; i++) {
	    if (ch == SSlist[i])
		return true;
	}
	return false;
    }
	
	
    /*
     *This function returns true if the residue given is not in the edge of a structure. 
     *i.e. the residues +/-1 has the same SS as it has.
     */
    public boolean notInEdge(int res) {
	char chp1=' ',chm1=' ',ch=' ';
	int i;
	for (i=0; (i<resNum.length) ; i++) {
	    if (resNum[i] == (res-1))
		chm1 = ss[i];
	    if (resNum[i] == (res))
		ch = ss[i];
	    if (resNum[i] == (res+1)) 
		chp1 = ss[i];
	}
        return (chp1 == ch) && (ch == chm1);
    }	 	
	
    /*
     *Returns the SS structure of residue number 'res'.
     *If this residue number doesn't exist - 'X' is returned.
     */
    public char SSofRes(int res) {
	int i;
	for (i=0; i<resNum.length ; i++) {
	    if (resNum[i] == res)
		return ss[i];
	}
	return 'X';
    } 


    /*
     *Returns the AA of residue number 'res'.
     *If this residue number doesn't exist - 'X' is returned.
     */
    public char AAofRes(int res) {
	int i;
	for (i=0; i<resNum.length ; i++) {
	    if (resNum[i] == res)
		return aa[i];
	}
	return 'X';
    } 
	
    /*
     *Returns the solvate accessibilty of residue number 'res'.
     *If this residue number doesn't exist - (-1) is returned.
     */
    public double ACCofRes(int res) {
	int i;
	for (i=0; i<resNum.length ; i++) {
	    if (resNum[i] == res)
		return solvACC[i];
	}
	return -1;
    } 


    /*
     *Returns the relative solvate accessibilty of residue number 'res'.
     *If this residue number doesn't exist - (-1) is returned.
     */
    public double relACCofRes(int res) {
	int i;
	for (i=0; i<resNum.length ; i++) {
	    if (resNum[i] == res)
		return solvACC[i]/fullSol[ResidueType.type(aa[i]).ordinal()];
	}
	return -1;
    }

    public String getACC(){
        String ans = "";
        for (int i=0; i<resNum.length ; i++)
	    if (resNum[i] > -999)
		ans = ans+  solvACC[i]+";";
	return ans;
    }

    public String getRelativeACC(){
        String ans = "";
        for (int i=0; i<resNum.length ; i++)
	    if (resNum[i] > -999){
		if ((int)((solvACC[i]/fullSol[ResidueType.type(aa[i]).ordinal()])+0.5) >= 1)
		    ans = ans+ "A" ;
		else
		    ans = ans+"B";
            }
	return ans;
    }

    public String getAA(){
        String ans = "";
        for (int i=0; i<resNum.length ; i++){
	    if (resNum[i] > -999)
		ans = ans+  aa[i];
            else
                throw new RuntimeException("res "+i+" is -999"+"\n"+
				           ans);
        }
	return ans;
    }

    public void printSS() {
	for (int i=0; i<resNum.length ; i++) {
	    if (resNum[i] > -999)
		System.out.println(ss[i] );
	    //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
	}
    }

    public String getSSInOneLine() {
        String ans = "";
        for (int i=0; i<resNum.length ; i++) {
	    if (resNum[i] > -999)
		ans = ans+  ss[i];
	    /*if(ss[i] == 'B' | ss[i] =='E')
	      ans = ans+'E';
	      else if(ss[i] == 'H')
	      ans = ans+'H';
	      else if(ss[i] == 'S' | ss[i] == 'G' | ss[i] == 'I' | ss[i] == 'T')
	      ans = ans+'A';
	      else
	      ans = ans+'C';*/

	    //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
	}
        return ans;
    }

    public String getSSInOneLineConventional() {
        String ans = "";
        for (int i=0; i<resNum.length ; i++) {
	    if (resNum[i] > -999)

		if(ss[i] == 'B' | ss[i] =='E')
		    ans = ans+'E';
		else if(ss[i] == 'H')
		    ans = ans+'H';
	    //else if(ss[i] == 'S' | ss[i] == 'G' | ss[i] == 'I' | ss[i] == 'T')
	    //  ans = ans+'A';//TODO - add this !
		else
		    ans = ans+'C';

	    //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
	}
        return ans;
    }

    public int getNumberOfBetaResidues(){
        int ans =0;
        for (int i=0; i<resNum.length ; i++) {
	          if(resNum[i] > -999 && ss[i] =='E')  //TODO(ss[i] == 'B' |
		                ans++;
          }
        return ans;
    }

    public void readDSSP()
    {
	int i=0;
	boolean cont = true;
	String line;
	try {
	    FileReader fr = new FileReader(fileName);
	    BufferedReader fdssp = new BufferedReader(fr);
	    line = fdssp.readLine();
	    while ((line != null) && (cont))
		{
		    if (line.length() > 25)
			if (line.substring(0,25).compareTo("  #  RESIDUE AA STRUCTURE") == 0)
	    		    {
    	    			cont = false;
			    }
			else{
			    String word = line.substring(50,65);
			    if(word.equalsIgnoreCase(" Parallel BRIDG"))
				parallelPrecentage  = Double.parseDouble(line.substring(6,10));
			    else if(word.equalsIgnoreCase("iParallel BRIDG"))
				antyparallelPrecentage = Double.parseDouble(line.substring(6,10));
			    else if(word.equalsIgnoreCase("E O(I)-->H-N(J)") )
				hbPrecentage = Double.parseDouble(line.substring(6,10));
			}
		    line = fdssp.readLine();
		} 
	    while (line != null)
		{
		    i++;
		    line = fdssp.readLine();
		    if((line != null) && (line.indexOf("!")) > -1) i--;
		} 		
	    fdssp.close();
	    fr = new FileReader(fileName);
	    fdssp = new BufferedReader(fr);
	    resNum = new int[i];
	    ss = new char[i];
	    aa = new char[i];
	    solvACC = new double[i];
	    for (int cc=0 ; cc<resNum.length ; cc++)
		resNum[cc] = -999;
	    i=0;	
	    cont = true;
	    line = fdssp.readLine();
	    while ((line != null) && (cont))
		{
		    if (line.length() > 25)
			if (line.substring(0,25).compareTo("  #  RESIDUE AA STRUCTURE") == 0)
			    {
				cont = false;
			    }
		    line = fdssp.readLine();
		} 
	    while (line != null)
		{
		    try {
			resNum[i] = (Integer.valueOf(line.substring(5, 10).trim()));
			aa[i] = line.charAt(13);
			ss[i] = line.charAt(16);
			if (ss[i] == ' ')
			    ss[i] = 'C';
			solvACC[i] = (Double.valueOf(line.substring(35, 38).trim()));
		    }
		    catch (Exception e) {
			i--;
		    }	
		    i++;
		    line = fdssp.readLine();
		} 		
	    fdssp.close();					
	} // of the try
	catch (Exception e) {
	    throw new RuntimeException(e);
	} // of the catch
    }
    
    public SequenceList sequenceList() {
	SequenceList seqList = new SequenceList();
	Sequence aaSeq = new ResidueSequence(getAA(),"aa sequence of "+fileName);
	seqList.add(aaSeq);
	Sequence ssSeq = new SecondaryStructureSequence(getSSInOneLineConventional() ,"ss sequence of "+fileName);
	seqList.add(ssSeq);
	Sequence accSeq = new AccesibilitySequence(getRelativeACC() ,"accesibility sequence of "+fileName);                seqList.add(accSeq);
	seqList.add(accSeq);
	return seqList;
    }
}


/*
 * The dssp line:
 1         2         3         4         5         6         7         8         9         0
 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
 18   23 A I  T <4 S+     0   0   64     -3,-3.4    -2,-0.2    -4,-0.2     7,-0.2   0.722 140.9  39.8-114.4 -49.1   11.0
   
*/
