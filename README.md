# apogee_rc_dimred

This project contains most of the code for the MSc thesis of Ignacio García Soriano (U Barcelona, June 2020) and the TFG thesis of Jaume Dolcet Monès (U Barcelona, 2021). Some of the code is also used in [Casamiquela et al. (in prep.:)](https://ui.adsabs.harvard.edu/abs/?????). It is inspired by and extends the analysis of [Anders et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A.125A/abstract) using APOGEE RC stars. 

To make the notebooks in py/ work, you will need to download [the DR16 APOGEE RC catalogue](https://data.sdss.org/sas/dr16/apogee/vac/apogee-rc/cat/apogee-rc-DR16.fits) (~260 MB) to the data/ folder. 
The corresponding data model is [here](https://data.sdss.org/datamodel/files/APOGEE_RC/cat/apogee-rc-DR16.html). The reference is [Bovy
et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...790..127B/abstract).

Explanation of the other files:
* data/occam_cluster-DR16.fits   -  Summary table from Donor+2020
* data/occam_member-DR16.fits    -  Cluster members from Donor+2020
* data/occam_member-DR16.fits    -  Cluster members from Donor+2020
* data/rc_clustermembers.txt     -  Crossmatch between occam_member-DR16.fits & apogee-rc-DR16.fits
* data/dimred_results/apogee_rc_dimred_hyperparametertest.fits   - Results of the t-SNE & umap runs for the DR16 RC sample
* data/dimred_results/apogee_rc_H_dimred_hyperparametertest.fits - Results of the t-SNE & umap runs for the DR16 RC sample (using abundances rel to H)
