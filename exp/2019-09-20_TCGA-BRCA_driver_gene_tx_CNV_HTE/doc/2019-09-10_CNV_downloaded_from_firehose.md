10/09/2019 21:57

CNV data downloaded from Firehose has the following structure

             Sample Chromosome     Start       End Num_Probes Segment_Mean
    68 TCGA-3C-AAAU          1   3218610  63469503      33369       0.1791
    69 TCGA-3C-AAAU          1  63471492  63472103          3      -0.8257
    70 TCGA-3C-AAAU          1  63472868  85632596      13663       0.2994
    71 TCGA-3C-AAAU          1  85635413  85878136        146       0.6498
    72 TCGA-3C-AAAU          1  85881278 149890533      20618       0.2781
    73 TCGA-3C-AAAU          1 149898951 150333087        165       2.2844
    
Teh segment_Mean here is already the log(tumor:)

[Biostar tutorial with GDC firehose](https://www.biostars.org/p/311199/#311746)


[Can I simply use the segment mean as a value for CNA?](https://www.biostars.org/p/257175/)

TCGA uses a normal reference set for tangent normalization of raw intensity values followed by segmentation using CBS algorithm. These reference normal set are different from matched normal (Blood derived sample from the patient) as well as Normal samples (derived from healthy individual). Find more here

Tumor samples, their corresponding matched-normal samples (Blood derived Normal) as well as Normal sample (Solid Tissue Normal), all are normalized against reference sample set.


CNV GISTIC data is retreived using [tutorial of RTCGAToolbox](http://bioconductor.org/packages/release/bioc/vignettes/RTCGAToolbox/inst/doc/RTCGAToolbox-vignette.html)

    firehose <- getFirehoseData("BRCA", runDate = "20160128", GISTIC = TRUE)
    CN <- getData(firehose, "GISTIC", "AllByGene") 
