
# Brainstorming for performing an HTE with clinically meaningful predictors

## Motivation

- forest of 50k is too large for interpretation
- people would probably want something that is more clinically oriented

## Ideas

- Use molecular indices as covariates and use tumor progression/metastatic status as outcome, some expressions are fine, such as HER2 or ERBB
- Preferably use real drugs/targeted therapy drugs as treatment
- Inferring patient outcome using existing data is not impossible, but that signal had to be strong enough in order for us to detect it.
- Using drug pathways as the `pathway_edges` in [PROPS](https://rdrr.io/bioc/PROPS/man/props.html)
