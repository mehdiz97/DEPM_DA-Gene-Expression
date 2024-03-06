The present study hypothesizes that a detailed analysis of TCGA data-specific
to Cholangiocarcinoma will reveal further insights into its molecular subtypes, prognosƟc biomarkers, and
potenƟal therapeuƟc targets. This hypothesis is grounded in the belief that a more nuanced understanding
of the genomic alteraƟons in Cholangiocarcinoma can lead to better diagnosis, prognosƟc, and treatment
strategies, ulƟmately improving paƟent outcomes.



The Deseq2 library facilitated the idenƟficaƟon of DifferenƟally Expressed Genes (DEGs) with an adjusted
p-value threshold of 0.05 and an absolute Fold Change (FC) greater than 1.2. Employing these DEGs, we
computed gene co-expression networks for the two conditions (cancerous and normal) using the
Spearman measure of similarity. We then constructed binary adjacency matrices for both cancerous and
normal samples based on Spearman similarity results. ThereaŌer, we identified network hubs (top 5% of
nodes by degree values) and compared the differenƟal co-expression network (Cancer vs. Normal) by
compuƟng the degree index. Finally, we analyzed the PaƟent Similarity Network using cancer gene
expression profiles and performed community detecƟon.
