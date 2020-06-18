# BSc-thesis-19-20 (computer science - 2019/2020)
Final BSc thesis is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the sixth semester of the undergraduate study. The main focus is to apply knowledge and skills obtained from Software Design Project course to recreate or improve existing methods which are widely used in bioinformatics. Under the supervision of prof. Mile Šikić and assistant prof. Krešimir Križanović, students will implement one such algorithm, thoroughly test it on simulated and real data, and formally encapsulate the whole process by writing and defending a thesis. Each student will have access to a private branch of this repository, on which this README will be updated with the specific task.

## Disclaimer
Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.

## Genome Scaffolding Based on HERA algorithm
HERA algorithm implementation based on published paper about HERA algorithm (https://www.biorxiv.org/content/early/2018/06/13/345983).

Implementation is heavily based on ScaRa tool(https://github.com/lbcb-sci/ScaRa) with minor changes in code.

## Instalation

```bash
git clone https://github.com/lbcb-edu/BSc-thesis-19-20.git
```
switch to branch ```matjuric```

To setup ```hera``` tool:
```bash
mkdir build && cd build
cmake ..
make
```

While in ```build``` folder, to run ```hera``` tool:
```bash
./hera contigs.fasta reads.fasta ovlpContRead.paf ovlpReadRead.paf
```

*Reads* file can be changed to FASTQ format with modifying code in ```src/gradi.cpp```.

Example test files for running ```hera``` tool are located in *```example```* folder.

Expected output assembled sequence is generated in ```build``` folder in file named *```output.fasta```*.

