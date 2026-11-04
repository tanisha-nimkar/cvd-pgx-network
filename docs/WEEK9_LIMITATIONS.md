# Week 9 Methodological Limitations

1. **Discrete Distance Granularity for Single-Gene Queries:**
   For a single pharmacogene ($|T|=1$), the shortest path distance $d_c$ is an integer ($0, 1, 2, \dots$). When testing against large disease modules ($>1,000$ genes), the random null distribution is heavily concentrated at distance 1. Consequently, standard deviations are narrow ($\sigma pprox 0.22 - 0.26$), causing modest integer shifts to produce large $z$-scores ($|z| > 3.5$) while true empirical tail probabilities remain non-significant ($p_{\text{emp}} pprox 0.06$).

2. **Direct Overlap Inflation:**
   Genes included directly within disease annotation sets (such as *CYP2C19* in CAD) achieve an artificial distance of $0.0000$. This reflects annotation co-membership rather than true network-mediated proximity across intermediate interactors.

3. **Literature Degree Inflation in STRING:**
   Highly studied cytochrome P450 enzymes (*CYP2C9* degree 89, *CYP2C19* degree 53) have elevated network connectivity driven by historical publication volume. Although degree-binning partially mitigates this bias, literature-rich nodes have smaller random search radii.

4. **Pharmacogenes Do Not Predict Clinical Drug Efficacy:**
   Network proximity between a pharmacokinetic enzyme/transporter and a disease module measures molecular contextualization. It does not imply that genetic variants in these enzymes modulate disease risk itself or that the drug possesses off-target efficacy.

## 5. Methodological Constraints Mandated by Supervisor Sign-Off (W9.10)
1. **Open Targets Scores Are Heuristic, Not Probabilities:**
   Open Targets association scores aggregate disparate data streams (GWAS, ClinVar, text mining, animal models) using heuristic weighting rules. They represent relative evidence ranking, NOT calibrated biological or statistical confidence probabilities that an association is true.
2. **CYP2C19 Is Not Statistically Significant Despite Negative z-Score:**
   Although *CYP2C19* has an observed distance of $d_c = 0.0000$ and a Gaussian approximation $z = -2.8520$, its empirical permutation $p$-value is $p_{\text{emp}} = 0.099900$ (with FDR $q = 0.2402$). Because $\sim 10\%$ of random degree-matched interactome nodes fall directly into the 1,006-gene CAD module by chance, *CYP2C19* fails the pre-registered significance threshold ($p < 0.05$). Under no circumstances may *CYP2C19* be declared statistically significant based on $z$-score alone.
