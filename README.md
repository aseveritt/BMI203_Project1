# Project 1 - Sequence Alignment
## Due 01/27/2021

![BuildStatus](https://github.com/aseveritt/BMI203_Project1/workflows/HW1/badge.svg?event=push)

In this assignment, you will implement two classical alignment algorithms and then evaluate each algorithm’s performance with a range of parameters. There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating alignments

### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m align
```

### testing
Testing is as simple as running
```
python -m pytest test/*
```
from the root directory of this project.


################
### Amanda additions:

- All functions/classes are contained within the single script align/algs.py
- An alignment can be run using the format:

import align as algs
sw = algs.SmithWaterman(gap_open=-4, gap_extension = -1, substitutionMatrix = "BLOSUM62")
b = sw.align(algs.FastaRecord("sequences/prot-0004.fa"), algs.FastaRecord("sequences/prot-0022.fa"), print_flag=True)
print(b.raw_score())
print(b.number_gaps())
nw = algs.NeedlemanWunsch(gap_open=-4, gap_extension = -0.5, substitutionMatrix = "BLOSUM50")
nw.align(algs.FastaRecord("sequences/prot-0004.fa"), algs.FastaRecord("sequences/prot-0022.fa"), print_flag=True)


- Script documentation is done through docstrings rendered with Sphinx:
open aseveritt/BMI203_Project1/docs/build/html/index.html

- Unit testing is all self contained within the single script as well that calls sequences already in sequences/ (directory not uploaded here)

- Plotting and images located in jupyter notebook:
open aseveritt/BMI203_Project1/Amanda_Everitt_BMI203_HW1.ipynb



