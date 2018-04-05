# ImmediateAncestry

This project finds immediate ancestors using an HMM by looking at the population frequencies of different populations at different SNPs.

This project requires a slightly modified form of the tmhmm package https://github.com/dansondergaard/tmhmm.py/tree/master/tmhmm. This modified form is not available online (yet).

The input is three types of file. Here is a small example of how they look like
```bash
$ cat seq1.txt
Chimpanzee_1
0 2 1 2 2 3 0
Chimpanzee_2
0 2 2 0 1 1 3
...
Chimpanzee_n
1 2 2 1 1 0 0

$ cat rho_chr1.txt
0.001 0.001 0.002 0.0011 0.0004 0.001

$cat chr1_freqs.txt
Chimpanzee_population_1
0.5 0.1 0.43 0.1 0.3 0.1 0.7
Chimpanzee_population_2
0.6 0.2 0.53 0.6 0.2 0.65 0.3
...
Chimpanzee_population_k
0.2 0.3 0.45 0.4 0.1 0.04 0.9
```
The program is run with
```bash
$ python ggrandparents-model.py --seq_files seq1.txt seq2.txt \
        --recomb_map rho_chr1.txt rho_chr2.txt \
        --allele_frequencies chr1_freqs.txt chr2_freqs.txt [OTHER_ARGS]
```
