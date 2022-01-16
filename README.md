## SFC_obesity ##
:large_blue_diamond: This is code for the paper **"Disrupted stepwise functional brain organization in overweight individuals"**<br />
:large_blue_diamond: **Paper link:** https://www.nature.com/articles/s42003-021-02957-7<br />
:large_blue_diamond: If you use this code, please cite the article.<br />
　　: Lee, Hyebin, et al. "Disrupted stepwise functional brain organization in overweight individuals." Communications Biology 5.1 (2022): 1-9.<br /><br />

✔ We obtained data from **[Enhanced Nathan Kline Institute-Rockland Sample (eNKI) database](http://fcon_1000.projects.nitrc.org/indi/enhanced/access.html)**.<br />
✔ Imaging data were preprocessed using **[FuNP pipeline](https://gitlab.com/by9433/funp)**.<br />
✔ **[Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/)** was used for SFC analysis.<br />
✔ **[BrainSpace Toolbox](https://brainspace.readthedocs.io/en/latest/#)** was used for visualization of cortex/subcortex features.<br /><br />

- **meanbold.m**　　　　　-----　to average BOLD signal for each atlas<br />
- **connectivity_seed.m**　　-----　to compute connectivity matrix from mean BOLD signals and define seed
　　　　　　　　　　　　　　regions<br />
- **stepwise_fc.m**　　　　　-----　to construct stepwise connectivity matrix<br />
- **group_analysis.m**　　　-----　to perform group analysis using degree of SFC<br />
- **TFEQ_association.m**　　-----　to associate degree of SFC with TFEQ scores<br />
- **main_figures.m**　　　　-----　to plot the results of group analysis<br /><br />

:pushpin: **copyrightⓒ 2021 All rights reserved by Hyebin Lee<br /><br />**
