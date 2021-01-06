
# Demographic comparison between two breast cancer cohorts

# Potential items to include for comparison

1. Age (mean, and range >/<50)
2. Gender
3. ER status
4. HER2 status
5. menopause status
6. Breast cancer subtype
7. Stage (I, IIA, IIB, IIIA, IIIB, IIIC)

# Issues
- TCGA HER2 and ER status requires extraction from GDC legacy file
- Chemotherapy and Other pharmacological therapy no longer considered, because this is not part of the main focus on demographic background checks
- include lymphnodes status?

# Report

| item                 | TCGA # | TCGA % | METABRIC # | METABRIC % |
| -------------------- | ------ | ------ | ---------- | ---------- |
| **gender**           |        |        |            |            |
| female               | 1085   |        | 2509       | 100.0      |
| male                 | 12     |        |            |            |
| age group            |        |        |            |            |
| [20,30)              | 9      | 0.83   | 17         | 0.68       |
| [30,40)              | 65     | 6.013  | 146        | 5.84       |
| [40,50)              | 214    | 19.8   | 404        | 16.17      |
| [50,60)              | 279    | 25.8   | 587        | 23.50      |
| [60,70)              | 285    | 26.36  | 714        | 28.58      |
| [70,80)              | 156    | 14.43  | 483        | 19.34      |
| [80,90)              | 64     | 5.92   | 140        | 5.60       |
| [90,100)             | 9      | 0.83   | 7          | 0.28       |
| **menopause**        |        |        |            |            |
| pre                  | 229    | 22.74  | 424        | 21.41      |
| post                 | 705    | 70.01  | 1556       | 78.59      |
| intermediate or peri | 73     | 7.1    | -          | -          |
| **ER status**        |        |        |            |            |
| positive             | 808    | 89.68  | 1817       | 72.42      |
| negative             | 238    | 26.42  | 609        | 24.27      |