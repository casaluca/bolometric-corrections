Bolometric Corrections and Synthetic Stellar Photometry 
-------------------------------------------------------

``bcutil.py`` contains python functions to interpolate bolometric corrections for Gaia DR2, DR3 and 2MASS at desired input stellar parameters (effective temperature, gravity and metallicity) and reddening (with the possibility of changing R<sub>V</sub> in the adopted extinction law). To use ``bcutil.py`` please download and unpack ``grid.tar.gz`` in the same directory.

``BCcodes``, ``BCtables`` and ``BCtables2`` are not needed to run bcutil.py. They contain codes and tables to compute bolometric corrections in many more photometric systems using a different set of programs. Also in this case, bolometric corrections can be computed at desired input stellar parameters and reddening (but without the possibility of varying R<sub>V</sub>).  

- ``BCcodes.tar.gz`` contains FORTRAN programs and shell scripts to be used with the tables below. 

- ``BCtables.tar.gz`` contains tables for the systems described in Casagrande & VandenBerg (2014): Johnson-Cousins, SDSS, HST, 2MASS.

- ``BCtables2.tar.gz`` contains tables for the systems described in Casagrande & VandenBerg (2018a, 2018b): JWST, PanSTARRS1, SkyMapper, Tycho, Hipparcos and Gaia DR2.

The content of the unpacked BC*.tar files must be in the same directory (i.e., the directories within BCtables and/or BCtables2 must be where files from BCcodes are located). For instructions and examples how to use the programs provided, please **read** INSTRUCTIONS.txt located in BCcodes.tar.gz.

Synthetic colours and reddening coefficients can be derived from bolometric corrections as described in Appendix A of Casagrande & VandenBerg (2014). Comparison of synthetic photometry with observations indicates that performances are overall good across optical and infrared bands (but please, refer to Section 4 of Casagrande & VandenBerg 2014 and 2018a for a discussion of successes and limitations).

Please cite [Casagrande & VandenBerg (2014)](http://adsabs.harvard.edu/abs/2014MNRAS.444..392C) and/or [Casagrande & VandenBerg (2018a,](http://adsabs.harvard.edu/abs/2018MNRAS.475.5023C) [2018b)](http://adsabs.harvard.edu/abs/2018MNRAS.479L.102C) if you find these packages useful for your research. 

If you are looking for empirical colour-effective temperature relations in the Gaia DR2 and DR3 systems, please check [here](https://github.com/casaluca/colte).

| ![My image](https://github.com/casaluca/bolometric-corrections/blob/master/DBC.jpg)
|:--:| 
| *Difference in bolometric corrections derived using Gaia DR2 vs DR3 passbands and zero-points.* |

