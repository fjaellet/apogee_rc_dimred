# apogee_rc_dimred

To make the noteboos in py/ work, you will need to download [the DR16 APOGEE RC catalogue](https://data.sdss.org/sas/dr16/apogee/vac/apogee-rc/cat/apogee-rc-DR16.fits) (~260 MB) to the data/ folder. 
The corresponding data model is [here](https://data.sdss.org/datamodel/files/APOGEE_RC/cat/apogee-rc-DR16.html).

Explanation of the other files:
* data/occam_cluster-DR16.fits   -  Summary table from Donor+2020
* data/occam_member-DR16.fits    -  Cluster members from Donor+2020
* data/occam_member-DR16.fits    -  Cluster members from Donor+2020
* data/rc_clustermembers.txt     -  Crossmatch between occam_member-DR16.fits & apogee-rc-DR16.fits
* data/dimred_results/apogee_rc_dimred_hyperparametertest.fits - Results of the t-SNE & umap runs for the DR16 RC sample
