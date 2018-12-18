Classification of CSES cases -- TODO

1. Read up on KNN algorithm
2. Sort CSES labels such that B is always the centrist
3. Run 3d scatterplots (bow tie result)
4. Fit regression plane: mAB and mBA predicting mCB (hypothesis: $beta_1$ = 1 and $beta_2$ = 0)
5. Compute "outlyingness" from plane
6. Add "ideal cases" for single-peaked, neutral, and divided majority
7. Run classification of CSES cases based on three "ideal cases" as training set
8. Check results. If meaningful, look at whether they are good predictors of patterns of optimal vote choice (level 0). 