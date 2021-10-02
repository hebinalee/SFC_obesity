## SFC_obesity - functions used in the main analysis ##
These functions are defined for convinience of analysis since are repeatedly used in main analysis.<br /><br />

- **binarize_conn.m**　　　-----　To binarize weighted connectivity matrix with threshold of 95%</br>
　　　　　　　　　　　　　SFC analysis requires binarized connectivity matrix as input, thus this should be preceded.<br />
- **compute_sfc.m**　　　-----　To calculate stepwise functional connectivity matrices from binarized connectivity matrix</br>
　　　　　　　　　　　　　Binarized connectivity matrix is considered as direct connections and indirect connections are inferred.</br>
- **mean_subcortical.m**　-----　To average feature values of 36 subcortical regions in BNA into 7 sub-structures of subcortex.</br>
　　　　　　　　　　　　　(amygdala, hippocampus, globus pallidus, nucleus accumbens, putamen, caudate, and thalamus)<br />
　　　　　　　　　　　　　The function gets input of feature vector/matrix with length of 246 (BNA) and returns with length of 14.<br /><br />

**copyrightⓒ 2021 All rights reserved by Hyebin Lee<br /><br />**
