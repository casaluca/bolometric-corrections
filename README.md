Bolometric Corrections and Synthetic Stellar Photometry 
-------------------------------------------------------

The programs and tables below allow users to compute bolometric corrections, synthetic colours and reddening coefficients for a variety of photometric systems given an input effective temperature, surface gravity and metallicity. 

- ``BCcodes.tar.gz`` contains FORTRAN programs and shell scripts to be used with the our tables below. 

- ``BCtables.tar.gz`` contains tables for the systems described in Casagrande & VandenBerg (2014): Johnson-Cousins, SDSS, HST, 2MASS.

- ``BCtables2.tar.gz`` contains tables for the systems described in Casagrande & VandenBerg (2018a, 2018b): JWST, PanSTARRS1, SkyMapper, Tycho, Hipparcos and Gaia.

For our programs to work, the content of the unpacked .tar files must be in the same directory (i.e., the directories within BCtables and/or BCtables2 must be where files from BCcodes are located). For instructions and examples how to use the computer programs we provide, please **read** INSTRUCTIONS.txt located in BCcodes.tar.gz.
Comparison of our synthetic photometry with observations indicates that performances are overall good across optical and infrared bands (but please, refer to Section 4 of Casagrande & VandenBerg 2014 and 2018a for a discussion of successes and limitations).

``bcutil.py`` contains python functions to interpolate bolometric corrections for Gaia DR2, DR3 and 2MASS at desired input stellar parameters and reddening (with the possibility of changing R<sub>V</sub> in the adopted extinction law). To use ``bcutil.py`` please download and unpack ``grid.tar.gz`` in the same directory.  

| ![My image](https://github.com/casaluca/bolometric-corrections/blob/master/DBC.jpg)
|:--:| 
| *Difference in bolometric corrections derived using Gaia DR2 vs DR3 passbands and zero-points.* |

Please cite [Casagrande & VandenBerg (2014)](http://adsabs.harvard.edu/abs/2014MNRAS.444..392C) and/or [Casagrande & VandenBerg (2018a,](http://adsabs.harvard.edu/abs/2018MNRAS.475.5023C) [2018b)](http://adsabs.harvard.edu/abs/2018MNRAS.479L.102C) if you find these packages useful for your research. 

If you are looking for empirical colour-effective temperature relations in the Gaia DR2 and DR3 system, please check [here](https://github.com/casaluca/colte).
