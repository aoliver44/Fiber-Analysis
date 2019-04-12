## Fiber analysis
## 
This markdown is here to help you make sense of all the files.

Lets begin with the raw data. The raw data off the sequencer took birth as four copies of R1 and R2
for each sample, which is how many lanes the illumina Nextseq has.

**STEP 1:** Filter the reads, decontaminate the reads, and apply taxonomy. All of these answers
youll find in the Full-fiber-basic.sh script. This script will take you from raw reads to 
decontaiminated and quality filtered reads. It will also run metaphlan and midas on the reads
too as an initial look at taxonomy. 

* In order to vizualize the taxonomy data i would refer you to the Fiber_project.md. Specifically
check out the Metaphlan section for the shiny app to get rel abundance graphs. For plotting the
MIDAS data, find the section called Relative abundance Species.


**STEP 2:** Next make a megahit directory and run the Fiber-script-2.sh

* This will merge all the metaphlan and midas data, and make a cross-assembly with all the decontaminated reads. 

**STEP 3:** 