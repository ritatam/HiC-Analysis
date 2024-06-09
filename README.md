## HiC Analysis
HiC analysis script to visualize HiC contact links via circos plot and provide bar plots per contig.

## Dependencies
1. python (>=3.12)
2. [pycircos](https://github.com/ponnhide/pyCircos)
3. [numpy](https://numpy.org)
4. [matplotlib](https://matplotlib.org)
5. [seaborn](https://seaborn.pydata.org)

## Usage
```sh
$ python src/hic_analysis.py
Generate HiC statistics and circos plot
src/hic_analysis.py <hap1.lst> <hap2.lst> <na.lst> <bed> <matrix> <outdir> <w>
<hap1.lst> path to a file consists the list of contig id prefixed with > from haplotype 1, one per line
<hap2.lst> path to a file consists the list of contig id prefixed with > from haplotype 2, one per line
<na.lst>   path to a file consists the list of contig id prefixed with > from unclassified, one per line
<bed>      path to HiC-Pro bed file
<matrix>   path to HiC-Pro matrix file
<outdir>   output directory
<w>        HiC-Pro window size
```

## Contact
john.luo@anu.edu.au
