package meshi.util;
public class Terminator {
    private boolean dead;
    public boolean dead() { return dead;}
    private String message;
    public String message() {return message;}
    public Terminator() {
	dead = false;
    }
    public void kill(String message) {
	this.message = message;
	dead = true;
    }
    
    public void reset() {
	message = "";
	dead = false;
    }
}
