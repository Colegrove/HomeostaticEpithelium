# Epidermis_Project_Final
Epidermis_Project with the final framework

Steps To Finish Epidermis Skin Model:

1) Neutral Story
	- Different Biopsy Sizes Need to be ran (1mm^2 for all right now, need 0.75mm^2 to 3.14mm^2)
	- Different mutation rates.

2) Non-Neutral Dynamics
	- 1mm^2 Biopsy Sizes
		- How does FIM and dN/dS compare with the Neutral results with Non-Neutral Dynamics
	- Heiko's suggestion (parameterization around fitness values.
	
3) How robust is the homeostasis to changes in parameters. MAYBE!!!
	- How large of a distribution can different cellular population take for parameters
	before loss of homeostasis
	
	
	
	
To Run Model:
Ages: 55, 58, 65, 73
(Currently Commented out) cellPositions.replicate_
2) lineage.replicate_
0) parents.replicate_
1) popSizes.replicate_
3) rLambda.replicate_
4) Size.
5) Time.
6) Mutation Rate Set (Default=0, must be present though)
7) SunDaysFreqency (Default=0, must be switched on in source files)
8) SunDaysDeathProb (Default=0, (theta) must be switched on in source files)
9) Replicate (Used for random seed along with Record time)

####Datasets:

- Neutral:
    - MutSet1 Dataset (Older dataset)
    - MutSet2-X Datasets for different Mutation Rates
    - Size difference = 0.75mm^2 to 1.0mm^2
- Selection:
    - TP53:
		- ThreeDNeutralTP53test.5302018.* are the death probability off, but for fishplots
		- threeDTP53deathproboff.5222018.* is the case with the death probability off for TP53
		- threeDTP53.Xtheta.Xsetfreq.6072018.* is the different values of theta and Sun Day Sets 30xDim
			- theta:
				- 0.001
				- 0.002
				- 0.006
				- 0.014
				- 0.033
				- 0.079
			- setFreq:
				- 5
				- 50
		- threeDTP53.Xtheta.Xsetfreq.6122018.* is the above dataset but for 100xDim
			- theta:
				- 0.014
				- 0.079
			- setFreq:
				- 5
				- 50
				
				
				
