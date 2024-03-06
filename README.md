The present study hypothesizes that a detailed analysis of TCGA data-specific
to Cholangiocarcinoma will reveal further insights into its molecular subtypes, prognostic biomarkers, and
potential therapeutic targets. This hypothesis is grounded in the belief that a more nuanced understanding
of the genomic alterations in Cholangiocarcinoma can lead to better diagnosis, prognostic, and treatment
strategies, ultimately improving patient outcomes.



The Deseq2 library facilitated the idenƟficaƟon of Differentially Expressed Genes (DEGs) with an adjusted
p-value threshold of 0.05 and an absolute Fold Change (FC) greater than 1.2. Employing these DEGs, we
computed gene co-expression networks for the two conditions (cancerous and normal) using the
Spearman measure of similarity. We then constructed binary adjacency matrices for both cancerous and
normal samples based on Spearman similarity results. Thereafter, we identified network hubs (top 5% of
nodes by degree values) and compared the differential co-expression network (Cancer vs. Normal) by
computing the degree index. Finally, we analyzed the Patient Similarity Network using cancer gene
expression profiles and performed community detection.
