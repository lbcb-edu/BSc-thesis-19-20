# BSc-thesis-19-20 (computer science - 2019/2020)
Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić and assistant prof. Krešimir Križanović, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Abstract
The main concern of this thesis is the clustering of modified nucleotides in nanopore sequenced RNA reads. The data for this thesis are reads and the reference from the tetrahymena ribozyme, a group I intron from Tetrahymena. The goal was to integrate existing tools such as Graphmap, Minimap and Minimap2, Ram, and Racon and link them with python scripts into a pipeline. To do this, all of the alignment tools were tested. A different approach for alignment was also tried out, aligning the reads in the signal domain. The pipeline was implemented in three parts: read selection, polishing and clustering of modified nucleotides.

## Usage
There are two versions of the pipeline. One is using [Ram](https://github.com/lbcb-sci/ram) as the main tool in the pipeline and other is using [Graphmap](https://github.com/isovic/graphmap). 
First one has better execution time and can be used when only there is only one transcript as a reference, and the second can be used even if there are more transcripts.
The main script is set to leave all the files, but if you want to delete unnecessary files and directories created through the process, you can change 6th parameter when calling reads_selection.sh from no to yes and 5th parameter when calling polish.sh. 
To speed up the process, reads can be split into smaller batches and clustering can be called in the end when results from the batches are connected. 

In order to start the pipeline, use the following command:
```bash
./main.sh <reads> <reference> <output_directory> <create_transcript_bins:yes/no>
```

## Disclaimer
Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.

## Acknowledgement
The data for the research was obtained from the Genome Institute of Singapore (GIS). The RNA used in the research is the tetrahymena ribozyme, a group I intron from Tetrahymena. The modifications on the RNA was done by Dr Wan Yue and Ashley Ghut Jong Aw. Special thanks to both of them.

I would also like to acknowledge the contribution of Dominik Batić for his help and start of this research while he was on GIS. 
