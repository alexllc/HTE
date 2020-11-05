
# 03/30 HPV status of the CD8 T cell HTE

03/30/2020 19:14

## Is the treatment effect of CD8 T cell proportion due to HPV status?

(A very brief proof)

Among the patients with known HPV status, we checked whether the extreme tau (>/<2) is due to HPV status.

We performed the following chi square test:

|  |  | tau |z value|
|  |  |>2|<2|
|--|--|--|--|
|HPV |+|13|17|
|status |-|26|56|

result:
```
Pearson's Chi-squared test with Yates' continuity correction

data:  a
X-squared = 0.84596, df = 1, p-value = 0.3577
```

So we can say that the difference in the observed number of extreme tau values are not due to HPV status, so this is not the underlying reason behind these extreme values. What are the factors? They will be in the variable importance, and you will have to dig into those.

