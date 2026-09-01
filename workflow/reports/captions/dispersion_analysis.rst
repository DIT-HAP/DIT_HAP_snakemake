DESeq2 dispersion estimates
=============================

Per-insertion dispersion estimates from DESeq2, plotted against the mean of
normalized counts (both axes log-scaled). Black points are the raw genewise
(maximum-likelihood) estimates, red is the fitted trend, and blue is the final
MAP estimate after shrinkage toward that trend. A healthy fit shows dispersion
decreasing with mean count and blue points pulled from the black cloud toward
the red curve; noisy insertions shrink the most. Genewise estimates far above
the trend that stay high after shrinkage are dispersion outliers and are left
unshrunk by design.
