package meshi.optimizers;

public class OptimizationTerminator {

	private boolean dead;
	private String message;

	public OptimizationTerminator() {
		dead = false;
	}

	public boolean dead() {
		return dead;
	}

	public void kill(String message) {
		this.message = message;
		dead = true;
	}

	public String message() {
		return message;
	}

	public void reset() {
		message = "";
		dead = false;
	}
}
