
									TO DO
									-----
									
- Theory:

	- Explain: Non-additive models (min, log, etc.)
	
	- Explain: Missing heritability
	
	- Bayesian factors: Why are these immune to multiple testing
	

- R model: Check if model CANNOT detect some models (e.g. |x1 - x2|, or  max(x1-x2,0), etc. )
	- Create cases and controls based on an ODDs ratio
	- Sample N/2 cases and N/2 controls
	- Perform logistic and LR test

- How are statistical models reflecting pathways? how about biological known interactions?	

- Looking for known interactions to validate the model. We need a test dataset.

- Assumptions: Loci 'i' and 'j' are in HWE (not in LD)

- Logistic regression using dummy variables

- BIC to select the right model: Additive, multiplicative, min, log? (max?)

- Power calculation: When can I really detect? 
	- Given and odds ratio / penetrance, which sample size do we need?
	- Power calculation standard is to find the loci at a 5% significance level, with an 80% probability

- Warning: Likelihood ratio test may not be Chi-Squared (it should assimptoticaly, but it may not)

- Typical test:
				mu = b0 + b1 x1 + b2 vx + b3 x1 x2

	WARNING: Should x1 and x2 be normalized?

- Flaws on the typical test?

				mu = b0 + b1 x1 + b2 vx + b3 x1 x2
				
				Is there a combination of x1 / x2 that this test cannot detect?
				I.e.: Similar to the XOR and linear unit problem (cannot be separated)
				 
- Another test strategy:

				mu = b0 + b1 x1 + b2 vx + b3 x1 x2
					
				Instead of testing for b3 = 0
				you could test b2 = b3 = 0

- Strategies for limiting the number of tests:

	- Use only coding variants: Impact moderate or higher
	- Use allele frequency threshold (AF or AC)
	- Use "known" interactions: PPI, pathways, gene sets
	- Interaction between Mendelian loci (this is known to exists, see nancy Cox's paper)
	 					 
- New methods: Bayesian statistics? MCMC? Other stuff?
