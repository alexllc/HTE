
## extracting results from object

    aldex = readRDS("./dat/ALDEx2_result.rds.gz")
    sig = filter(aldex, abs(effect) > 1 & abs(overlap) < 0.01 ) # default filter setting
    head(sig)[, c(1:3, 1219:1226)]
    sig = filter(aldex, abs(effect) > 1 & abs(overlap) < 0.015  # more generous setting
    # BRCA 414 sig genes vs 5 genes in default filter

The recommended $|\Delta_{Q_0}^{A}| \leq 0.01$ and $|\Delta_R| \geq 1.5$ requirement only resulted in 5 DE genes in the TCGA-BRCA dataset:

>As a theoretical minimum, the criterion that $|\Delta_{R}^{50}| \geq 1$ can be used to select genes where expected between-condition changes are of the same order or greater than the expected within-condition changes. However, this figure shows that $|\Delta_{R}^{50}| \geq 1$ is too inclusive.

                    rab.all rab.win.NT rab.win.TP  diff.btw diff.win    effect
    ENSG00000099953 10.384380   4.481162  10.653495  6.020681 1.999537  2.844156
    ENSG00000154736  6.431999   9.662481   6.260838 -3.382345 1.434293 -2.276067
    ENSG00000161649  2.665743   8.773723   2.331235 -6.499735 3.108106 -2.000147
    ENSG00000165197  1.742301   7.848220   1.480719 -6.147180 2.340500 -2.511602
    ENSG00000168497  5.288417   9.742912   5.053344 -4.682565 1.986490 -2.249140
                        overlap         we.ep        we.eBH        wi.ep
    ENSG00000099953 0.005798912 2.920229e-127 3.443559e-124 6.162514e-67
    ENSG00000154736 0.007000059 2.098791e-116 1.600772e-113 1.357598e-66
    ENSG00000161649 0.007398576 8.006541e-129 9.690166e-126 6.354361e-67
    ENSG00000165197 0.005798912 8.286394e-134 1.330798e-130 2.371859e-67
    ENSG00000168497 0.006598743 3.020443e-110 1.831472e-107 6.228782e-67
                        wi.eBH
    ENSG00000099953 9.150922e-63
    ENSG00000154736 1.492571e-62
    ENSG00000161649 9.151790e-63
    ENSG00000165197 9.150922e-63
    ENSG00000168497 9.150922e-63

