Morphic data analysis 


This analysis only covers the cell pellets and supernatant ran in March of 2023. For media experiments, contact Amnah Siddiqi.

This is the workflow for all the experiments:

1. Convert .raw to .mzML
2. Analyze .mzML to feature tables using Asari
3. Examine QCQA metrics on Asari feature tables
4. Mark samples that are not experimental samples (e.g., pooled QC, blanks, etc.)
5. Examine QCQA metrics on Asari feature tables using just the experimental sample subset
6. Drop any sample that has a missing feature Z-score greater than 2.5 std. 
7. Preprocess the feature table by performing TIC normalization on features in >80% of samples, dropping features in less than 80% of samples, and imputing missing feature values by half of the minimum value across the experiment. 
8. Perform PCA and t-SNE on the features to examine overall patterns. PCA assumes linearity, t-SNE is more complex but does not assume linearity. 
9. Next analyze each feature table using a linear model. Construct the clustermaps using the omnibus p-values and the per-term p-values. Correct these p-values in a bonferroni-esque method.
10. Using a t-test, compare each KO, media type, differentiation against the appropriate control. Find any features that are statistically significant after FDR correction (Benjamini-Hochberg)
11. Perform differential abundance analysis for lipid categories (obviously this only applies in lipidomics experiments). Examine if a lipid category is up or down changed for the given comparison (same groups as the t-test) using a Fischer exact test, two-sided. Correct raw p-values using Benjamini-Hochberg and report. 

The results from these files are located as follows:

Figures from 8 and 9 are stored in ./figures.
Significant features from step 10 are located in /annotated_significant_features. These values are also used to generate the intermediates in ./mummichog_data
Differential abundant analysis of lipid categories is located in ./lipidomics_differential_abudance_analysis/

Annotation Details

Annotation was performed in Metabolomics experiments using the Human Metabolome Database (HMDB) and using the LipidMaps Structure Database (LMSD) for Lipidomics experiments. Matching was done using m/z value alone and a 5ppm tolerance 

!!! MS1 matching only is not reliable for compound identity. Formulas are likely to be correct but the compound assignment may not be. This is the nature of mass spectrometry. It is fine to report them as hypotheses, but be careful making any definitive statement claiming a compound occurs or not. Always use language such as, we had annotations to compounds in the blah blah pathway or an annotation to Vitamin D3 was found to be differenitally abundant in... !!! you get the point. 


TODO:

Finish Mummichog analysis
Write methods

Analyze DmPA
Analyze DnHZ and DnCl data.

FUTURE:

The 80% cutoff rule may eliminate features that are only in a subset of samples but also does a great job eliminating noise. It is a trade off, but some form of group filtering may be better suited for finding significant metabolites.