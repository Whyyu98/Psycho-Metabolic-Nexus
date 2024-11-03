# Shared Genetic Architecture and Bidirectional Clinical Risks within the Psycho-Metabolic Nexus
## 1. Clinical Risk Assessment (epi)
### interaction_univariate_logistic.r
To investigate the bidirectional clinical risks, psychiatric or metabolic profiles in UKBB were focused on. We analysed 6 available psychiatric disorders and 24 metabolic profiles, treating them as either exposures or outcomes. Logistic and linear regression models with FDR correction were employed to estimate these clinical risks. Confounders including the top 10 principal components, age, sex, and average total household income before tax were adjusted.
## 2.GWAS META
### lambda_meta.R
We ensured minimal sample overlap between participating studies and homogeneity of effect sizes by calculating λmeta.
### metal.txt
For meta-analyses without sample overlap, we used the inverse variance weighted (IVW) fixed-effects method of the METAL software.
## 3. Genetic Correlation
### LDSC_for.sh ,LDSC_result.R, and GNOVA_for.sh
Genome-wide genetic correlations between psychiatric disorders and metabolic dysregulations were examined using linkage disequilibrium score regression (LDSC) and the Genetic Covariance Analyzer (GNOVA).
### s-LDSC.sh
To explore correlations across various cell-and-tissue types, stratified LD score regression (s-LDSC) was performed. This analysis included data from 5 datasets: 39 brain tissue cells from Zeisel et al., 292 immune cells from ImmGen, 152 cell types and tissue types from Franke et al., and 53 cell types from GTEx V8.
## 4. Local Genetic Correlation
### LAVA_for.R
For local genetic correlation among 2495 regions, Local Analysis of Variance Annotation (LAVA) was used to define the regions with significant heritability in both trait and correlation.
### GWAS-PW_for.sh
Pairwise GWAS (GWAS-PW) was used to validate LAVA results, applying a threshold of posterior probability of A3/4 (PP.A3/4) > 0.8.
## 5.Pleiotropic Loci
### MTAG_for.sh, and CPASSOC_for.R
To identify pleiotropic loci shared between traits, Multi-Trait Analysis of GWAS (MTAG) was performed, complemented by Cross-Phenotype Association (CPASSOC) for sensitivity analysis.
### clump.sh
Independent genome-wide significant pleiotropic loci were identified using PLINK1.9 (r2 < .02 and 500-kb window).
### MTAG_CPASSOC.R and SNV+GWAS.R
Significant pleiotropic SNVs were defined as P < 1 × 10−4 in GWAS studies and P < 5 × 10−8 in the meta-analysis, ensuring the meta-analysis P value was smaller than the original value.
### COLOC.R, susie.R, and hyprcoloc.R
For these potential pleiotropic loci (markers in the ±100 kb range), bayesian colocalisation (COLOC), SuSiE, and HyPrColoc determined the probability of shared causal variants.
## 6. GARFIELD
### garfield_input_for.R, Garfield, and garfield_output_for.R
The functional impact of shared genetic signalling was investigated using GARFIELD, focusing on SNPs identified in MTAG. GARFIELD compared annotated markers from specific cell types and adjusted for various genetic factors.  Under Bonferroni correction based on the number of annotations from each tissue, relevant tissues were grouped into 26 categories. Proportions of significant annotations in each category were examined with Pearson correlation test.
## 7. Identification of Pleiotropic Genes
### MAGMA_for.sh, mBAT.sh, POPS_for.sh, SMR.sh, and TWAS.sh
To identify associations at loci identified in cross-trait meta-analyses for each trait pair, we conducted Summary-data-based Mendelian Randomization (SMR), functional Summary-Based Imputation (FUSION), Polygenic Priority Score (PoPS), and mBAT-combo (mBAT) analyses.
## 8. Mendelian Randomization (MR)
### MR_for.R, and MRlap_for.R
To infer causality between psychiatric disorders and metabolic profiles, MR was performed. IVW with FDR correction was employed as the primary method. As the sensitive tests, MR-Egger, weighted median, weighted mode, and MRlap were applied to validate the results.
## 9. Polygenic Risk Score (PRS)
### get_prs_validation.r
To validate the genetic role in comorbidity, PRS was constructed using Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM) based on GWAS statistics for psychiatric disorders and metabolic profiles. Continuous traits were modelled with Gaussian linear models and binary traits with generalised logistic regression, testing PRS effectiveness for each trait. Adjustments were made for age, sex, the top 10 genetic principal components, and income. Models with two PRS were assessed for improved predictive performance beyond outcome-specific PRS. Pearson's R² and McFadden's pseudo-R² evaluated model performance with FDR correction. Bootstrap calculated 95% confidence intervals (CI) for Pearson's R² and area under the curve for binary traits.  
DBSLMM: https://github.com/biostat0903/DBSLMM
## Other Web Tools
GWAShug: http://www.gwashug.com/  
PGSFusion: http://www.pgsfusion.net/#/  
PGS-Depot: http://www.pgsdepot.net/  


## Q&A
If you have any questions, please contact the author, Yu Feng (fengyu19981121@gmail.com).
