# Week 9 Methodological Limitations

1. **Discrete Distance Granularity for Single-Gene Queries:**
   For a single pharmacogene ($|T|=1$), the shortest path distance $d_c$ is an integer ($0, 1, 2, \dots$). When testing against large disease modules ($>1,000$ genes), the random null distribution is heavily concentrated at distance 1. Consequently, standard deviations are narrow ($\sigma pprox 0.22 - 0.26$), causing modest integer shifts to produce large $z$-scores ($|z| > 3.5$) while true empirical tail probabilities remain non-significant ($p_{\text{emp}} pprox 0.06$).

2. **Direct Overlap Inflation:**
   Genes included directly within disease annotation sets (such as *CYP2C19* in CAD) achieve an artificial distance of $0.0000$. This reflects annotation co-membership rather than true network-mediated proximity across intermediate interactors.

3. **Literature Degree Inflation in STRING:**
   Highly studied cytochrome P450 enzymes (*CYP2C9* degree 89, *CYP2C19* degree 53) have elevated network connectivity driven by historical publication volume. Although degree-binning partially mitigates this bias, literature-rich nodes have smaller random search radii.

4. **Pharmacogenes Do Not Predict Clinical Drug Efficacy:**
   Network proximity between a pharmacokinetic enzyme/transporter and a disease module measures molecular contextualization. It does not imply that genetic variants in these enzymes modulate disease risk itself or that the drug possesses off-target efficacy.
