Anti-Viral Antibody Response Deconvolution Algorithm (AVARDA) v1.1
Date: February 9th, 2018
@Architects: Ben Larman, Daniel Monaco, Sanjay Kottapalli
@Developers: Sanjay Kottapalli, Tiezheng Yuan
Language: Python version 3.6
For bugs or help, contact Sanjay: skottap3@jhu.edu

-----------------------------------------------------------------------------

PURPOSE: This pipeline takes VirScan peptide enrichment Z-scores and 
calculates confidence of infection for ~500 viruses known to infect humans.

REQUIRED MODULES/PACKAGES:
pandas (v0.19.2 or later)
NumPy (v1.11.3 or later)
NetworkX (v1.11 or later)
SciPy (v0.18.1 or later)
StatsModels (v0.6.1 or later)

INSTALLATION/SETUP:
Clone the GitHub repository to the directory of your choice. The /bin/ folder
contains the python code, the /ref_seq/ folder contains internal VirScan
files read by the pipeline, the /input/ folder contains the user parameters
file ('variables_virus.txt') and the Z-scores file(s) to be analyzed, and
finally the /results/ folder, generated after the first run, contains the 
results of the run (contained within a newly generated, timestamped folder). 
Avoid editing any files in the /ref_seq/ or /bin/ folders as they are 
specific to the pipeline and should remain static.

USER DIRECTIONS:
Open the 'variables_virus.txt' file within the input directory. Change the
value of the 'dir_home' variable to reflect your local absolute path to the
cloned 'AVARDA' directory. Likewise, change the name of the zscore file
to be analyzed. The optimized 'Z_threshold', 'p_threshold', and 'x_threshold'
values have been entered but can be changed as desired. The 'bh_threshold'
value determines the threshold for adjusted p-value for which an infection 
will show up in the 'results_summary.txt' file, and can also be adjusted.

After the variables_virus.txt file has been set up, simply execute the 
mainfile 'main.py' in the /bin/ directory. NOTE: Do not execute multiple runs
of the pipeline in the same minute or else output files may be overwritten.

INTERPRETING RESULTS:
For each run, a date- and time-stamped subfolder within /results/ will be 
generated. The subfolder /sample_networks/ will contain '.graphml' file 
representations of enriched peptide 7 amino acid-overlap networks for each 
sample, numbered starting from 0 according to the header of the Z-scores file. 
These files can be opened with graph software such as Cytoscape, yEd,
Gephi, etc. The remainder of the files for the run will also be available in
the aforementioned datestamped folder. 'results_summary.txt' is the most 
relevant, and contains a tab-delimited list of all samples and their detected
viral infections.

-------------------------------------------------------------------------------