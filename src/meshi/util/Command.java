package meshi.util;
import java.util.*;
import meshi.util.string.*;
public class Command implements KeyWords{
    private String line;
    private String commandsComment;
    private StringList words;

    public Command() {}

    public Command(String line, String commandsComment) {
	this.line = line;
	this.commandsComment = commandsComment;
	words = new StringList();
	StringTokenizer tokens = new StringTokenizer(line);
	while (tokens.hasMoreTokens())
	    words.add(tokens.nextToken());
    }

    public String firstWord() {
	if (words.size() < 1) throw new RuntimeException("No first word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	return words.get(0);
    }

    public String secondWord() {
	if (words.size() < 2) throw new RuntimeException("No second word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	return words.get(1);
    }

    public boolean secondWordExists() {
	return (words.size() >= 2);
    }

    public String thirdWord() {
	if (words.size() < 3) throw new RuntimeException("No third word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	return words.get(2);
    }

    public int secondWordInt() {
	if (words.size() < 2) throw new RuntimeException("No second word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return new Integer(words.get(1));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse int in second word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public int thirdWordInt() {
	if (words.size() < 3)
        throw new RuntimeException("No third word in "+line+"\n"+"comandsFile "+commandsComment);
	try {
	    return new Integer(words.get(2));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse int in third word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public int thirdWordPositiveInt() {
	if (words.size() < 3) return -1;
	try {
	    return new Integer(words.get(2));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse int in third word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public int fourthWordPositiveInt() {
	if (words.size() < 4) return -1;
	try {
	    return new Integer(words.get(3));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse int in third word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public double secondWordDouble() {
	if (words.size() < 2) throw new RuntimeException("No second word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return new Double(words.get(1));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse double in second word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public double thirdWordDouble() {
	if (words.size() < 3) throw new RuntimeException("No third word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return new Double(words.get(2));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse double in third word of "+line+"\n"+
							 words.get(2)+"\n"+
							 "in "+commandsComment);}
    }

    public double fourthWordDouble() {
	if (words.size() < 4) throw new RuntimeException("No third word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return new Double(words.get(3));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse double in third word of "+line+"\n"+
							 words.get(3)+"\n"+
							 "in "+commandsComment);}
    }

    public double fifthWordDouble() {
	if (words.size() < 5) throw new RuntimeException("No third word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return new Double(words.get(4));
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse double in third word of "+line+"\n"+
							 words.get(4)+"\n"+
							 "in "+commandsComment);}
    }

    
    public boolean secondWordOn() {
	if (words.size() < 2) throw new RuntimeException("No second word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return ( words.get(1)).equals(ON.key);
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse boolean in second word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public boolean thirdWordOn() {
	if (words.size() < 3) throw new RuntimeException("No third word in "+line+"\n"+
							 "comandsFile "+commandsComment);
	try {
	    return (words.get(2)).equals(ON.key);
	}
	catch (Exception e) { throw new RuntimeException("Failed to parse boolean in third word of "+line+"+\n"+
							 "in "+commandsComment);}
    }

    public void printLine(String prompt) {
	System.out.println(prompt+"  "+line);
    }
    
    public String getLine(){return line;}
    public boolean hasSecondWord(){
        return words.size()>=2;
    }
    
    public String toString() {return line;}
}

    
