Modeling epithelial competition in Fanconi Anemia gene therapy: https://doi.org/10.1101/2025.02.26.640284

Oral epithelial model adapted from Schenk et al (2022): https://doi.org/10.1073/pnas.2006487119


![Pipeline schematic](model_pipeline_schematic.png)


The computational pipeline uses a Snakemake workflow to manage model simulation runs. The Snakemake file takes as input experimental specific parameters such as the number of simulation years, number of replicates, persistence coefficients for corrected and non-corrected cells, and mutation rates for each genotype. These parameters are passed to the model itself, a Java-based .jar artifact to execute each model simulation. Once model simulations are run, output data is processed using analysis scripts.
