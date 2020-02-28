### Idea about the deconvolution based on MHL

In my previous report, if I mix MHL between WBC and CCT with different
recombination, the deconvolution effect is quite good. However, in the
reality, it is not mixture of MHL but mixture of haplotype. Here, I want
to observe the change of the MHL with the different propotion of the
"CCCC" reads in the plasma samples population ("TTTT") in 1% to 100%. In
the beginning, we have worried that the haplotype mixture will not
provide linear change to MHL, the influence to the deconvolution will be
huge. and I want to find a function to transfer this
non-linear/generalized linear to linear change so that deconvolution can
conducted with the transform data.

In order to make all the analysis reproducible, I start to use Rmd to do
the analysis from now on. It will be easy to repeat the analysis
anytime, anywhere and by anybody. I will creat a perl script to
reconstructed Rmd so that I could be easy to input to our wiki.

### Simulation analysis by mixture of 4-C-Haplotype

This figure show the principle, how to mix the CCT haplotype (black)
with WBC haplotype (white) by different propotion.
![](Supplementary._MHL_simulation_based_sampled_methylation_and_unmethyation_reads_files/figure-markdown_strict/haplotypeMixture-1.png)
Relationship between MHL and 5mC.We can obsever the relationship between
haplotype distribution, MHL and 5mC with the following table

    ##     PROPC PROPT  amC MHL
    ## 1       0   100 0.00 0.0
    ## 2       1    99 0.01 0.5
    ## 3       2    98 0.02 0.5
    ## 4       3    97 0.03 0.5
    ## 5       4    96 0.04 0.5
    ## 6       5    95 0.05 0.5
    ## 7       6    94 0.06 0.5
    ## 8       7    93 0.07 0.5
    ## 9       8    92 0.08 0.5
    ## 10      9    91 0.09 0.5
    ## 11     10    90 0.10 0.5
    ## 12     11    89 0.11 0.5
    ## 13     12    88 0.12 0.5
    ## 14     13    87 0.13 0.5
    ## 15     14    86 0.14 0.5
    ## 16     15    85 0.15 0.5
    ## 17     16    84 0.16 0.5
    ## 18     17    83 0.17 0.5
    ## 19     18    82 0.18 0.5
    ## 20     19    81 0.19 0.5
    ## 21     20    80 0.20 0.5
    ## 22     21    79 0.21 0.5
    ## 23     22    78 0.22 0.5
    ## 24     23    77 0.23 0.5
    ## 25     24    76 0.24 0.5
    ## 26     25    75 0.25 0.5
    ## 27     26    74 0.26 0.5
    ## 28     27    73 0.27 0.5
    ## 29     28    72 0.28 0.5
    ## 30     29    71 0.29 0.5
    ## 31     30    70 0.30 0.5
    ## 32     31    69 0.31 0.5
    ## 33     32    68 0.32 0.5
    ## 34     33    67 0.33 0.5
    ## 35     34    66 0.34 0.5
    ## 36     35    65 0.35 0.5
    ## 37     36    64 0.36 0.5
    ## 38     37    63 0.37 0.5
    ## 39     38    62 0.38 0.5
    ## 40     39    61 0.39 0.5
    ## 41     40    60 0.40 0.5
    ## 42     41    59 0.41 0.5
    ## 43     42    58 0.42 0.5
    ## 44     43    57 0.43 0.5
    ## 45     44    56 0.44 0.5
    ## 46     45    55 0.45 0.5
    ## 47     46    54 0.46 0.5
    ## 48     47    53 0.47 0.5
    ## 49     48    52 0.48 0.5
    ## 50     49    51 0.49 0.5
    ## 51     50    50 0.50 0.5
    ## 52     51    49 0.51 0.5
    ## 53     52    48 0.52 0.5
    ## 54     53    47 0.53 0.5
    ## 55     54    46 0.54 0.5
    ## 56     55    45 0.55 0.5
    ## 57     56    44 0.56 0.5
    ## 58     57    43 0.57 0.5
    ## 59     58    42 0.58 0.5
    ## 60     59    41 0.59 0.5
    ## 61     60    40 0.60 0.5
    ## 62     61    39 0.61 0.5
    ## 63     62    38 0.62 0.5
    ## 64     63    37 0.63 0.5
    ## 65     64    36 0.64 0.5
    ## 66     65    35 0.65 0.5
    ## 67     66    34 0.66 0.5
    ## 68     67    33 0.67 0.5
    ## 69     68    32 0.68 0.5
    ## 70     69    31 0.69 0.5
    ## 71     70    30 0.70 0.5
    ## 72     71    29 0.71 0.5
    ## 73     72    28 0.72 0.5
    ## 74     73    27 0.73 0.5
    ## 75     74    26 0.74 0.5
    ## 76     75    25 0.75 0.5
    ## 77     76    24 0.76 0.5
    ## 78     77    23 0.77 0.5
    ## 79     78    22 0.78 0.5
    ## 80     79    21 0.79 0.5
    ## 81     80    20 0.80 0.5
    ## 82     81    19 0.81 0.5
    ## 83     82    18 0.82 0.5
    ## 84     83    17 0.83 0.5
    ## 85     84    16 0.84 0.5
    ## 86     85    15 0.85 0.5
    ## 87     86    14 0.86 0.5
    ## 88     87    13 0.87 0.5
    ## 89     88    12 0.88 0.5
    ## 90     89    11 0.89 0.5
    ## 91     90    10 0.90 0.5
    ## 92     91     9 0.91 0.5
    ## 93     92     8 0.92 0.5
    ## 94     93     7 0.93 0.5
    ## 95     94     6 0.94 0.5
    ## 96     95     5 0.95 0.5
    ## 97     96     4 0.96 0.5
    ## 98     97     3 0.97 0.5
    ## 99     98     2 0.98 0.5
    ## 100    99     1 0.99 0.5
    ## 101   100     0 1.00 1.0

And you can check the Figur as the following:

![](Supplementary._MHL_simulation_based_sampled_methylation_and_unmethyation_reads_files/figure-markdown_strict/MFvsMHL-1.png)

Our MHL metric is difficult to be used in deconvolution, since the value
is none linear or generalzied linear increasement as the propotion of
the methylated haplotypes.

### Real data analysis by mixture of haplotype

for the real mixture data, some regions will have linear change and
majority them are chaotic. I have tried many ways to find a function to
transfer them as a generalized linear way. However, it seems it is hard
to do that since in the definition of the MHL, there is a unique
operation for the haplotype therefore the relatiohsip will be
unexpectable.  
![](Supplementary._MHL_simulation_based_sampled_methylation_and_unmethyation_reads_files/figure-markdown_strict/unnamed-chunk-4-1.png)

### Conclusion

Again, really data result also show our MHL metric is difficult to be
used in deconvolution, since the value is none linear or generalzied
linear increasement as the propotion of the methylated haplotypes.

pandoc("")
