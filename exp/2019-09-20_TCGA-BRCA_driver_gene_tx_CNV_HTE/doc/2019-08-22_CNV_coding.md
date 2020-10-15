22/08/2019 10:52

Use [Davoli et al.](http://www.sciencedirect.com/science/article/pii/S0092867413012877) list of TSG and OG.
We start with the TSG and OG by PAN cancer calculation (estimated to cover about ~70% of all the TSG and OG) first, we can move on to do tissue specific TSG and OG.

***

Sort by column and then select top n
head(TSG[order(TSG$TUSON_q_value_TSG, decreasing= F),], n = 300)

What is the number after gene symbol?
`[1] "ENSG00000008128.21" "ENSG00000008130.14" "ENSG00000067606.14" "ENSG00000078369.16"
   [5] "ENSG00000078808.15" "ENSG00000107404.16" "ENSG00000116151.12" "ENSG00000127054.17"
   [9] "ENSG00000131584.17" "ENSG00000131591.16" "ENSG00000142606.14" "ENSG00000142609.16"`
   

The number after ensembl ID is the version number.

But I'm going to just try to convert directlty into HUGO symbols according to a [biostar post](https://www.biostars.org/p/200810/)

***

From [biostar post](https://www.biostars.org/p/343204/)

"You just took the tumour samples, merged them together, and then summarised copy number per gene? There's nothing wrong with that; however, if you try to publish it, a reviewer would ask if these copy number events are somatic or germline. Without subtracting the normal copy number profiles, you cannot know whether or not these are true somatic events."


By using the *copy number estimation* it's previous step of generating *Masked Copy Number Segment* has already subtracted the normal sample (see [GDC documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/))

-> so all the "CNV" I have here is already *somatic*, it should be called **sCNA** to be proper


***

Expresion data

- only selecting primary solid tumors (-01)
- Convert ensembl ID to hugo symbol


***

Ensembl ID unique but gene ID not unique
19729 unique ensembl ID, but only 18729 unique gene names
-> So let's translate them after getting a result


***

Do I have to disable the execution in source script?


wifi password: cuhksbs624a 
