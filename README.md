## Fiber analysis
## 
This markdown is here to help you make sense of all the files.

Lets begin with the raw data. The raw data off the sequencer took birth as four copies of R1 and R2
for each sample, which is how many lanes the illumina Nextseq has.

**STEP 1:** Filter the reads, decontaminate the reads, and apply taxonomy. All of these answers
youll find in the Full-fiber-basic.sh script

**STEP 2:** Next i had to merge some files, like the metaphlan ones and the Midas outputs

**STEP 3:** Stats on those files, which are found in the Fiber_project.Rmd 