"""
Functions to interpolate MARCS tables of bolometric corrections (BCs) at desired
input Teff, [Fe/H], log(g) and E(B-V). 

Akima splines are used to interpolate in Teff, [Fe/H] and log(g), plus linear 
interpolation in E(B-V). These are the same prescriptions used for the Fortran 
routines in BCcodes.tar. Agreement with those routines is typically at millimag
level, expect for models at the coolest Teff and log(g)<2 (currently excluded 
by this routine, see rme tuple below).  

By default, the routine will compute BCs for the Gaia DR3 system and save those 
into a file named bctable_Rv31.csv. In presence of reddening, the extinction 
law of Cardelli/O'Donnell is used with Rv=3.1 (same as in the BCcodes and 
BCtables packages). If BCs cannot be computed for a star, nan values will be 
returned. 

REQUIRED INPUT (scalar or array): 
sid:     star ID/name (list also OK)
teff:    effective temperature
logg:    surface gravity
feh:     metallicity [Fe/H]
ebv:     reddening. Accepted values 0 <= ebv < 0.72

OPTIONAL PARAMETERS:
Rv:      3.1 by default (standard extinction law). Accepted values: 3.1 or 2.5

filters: by default only BCs for Gaia DR3 are computed. Available filters:
         G3     Gaia DR3 G 
         BP3    Gaia DR3 BP
         RP3    Gaia DR3 RP
         G2     Gaia DR2 G 
         BP2    Gaia DR2 BP
         RP2    Gaia DR2 RP
         J2M    2MASS J
         H2M    2MASS H
         K2M    2MASS K

         Gaia DR3 are obtained using updated passband and Vega zero points from 
         https://www.cosmos.esa.int/web/gaia/edr3-passbands
         Gaia DR2 are obtained with revised transmission curves and data 
         processing zero-points. Same as 27 = Gaia_DR2 in selectbc.data in 
         the BCcodes package.

outfile: name for output file. By default bctable_Rv?.?.csv
retarr:  True to return an array with BCs

USAGE:
 import bcutil

EXAMPLE 1: derive BCs for stars of know input parameters in the Gaia DR3 system,
assuming standard extinction law. Results are saved in bctable_Rv3.1.csv:

 bcutil.bcstar(sid,teff,logg,feh,ebv)

EXAMPLE 2: derive BCs for a star having Teff=5300, log(g)=3.89, [Fe/H]=-0.23, 
E(B-V)=0.04 in Gaia DR3 G, RP and 2MASS K with Rv=2.5 extinction law. 
Results are saved in output file 'mystar.csv', and returned into array mybc:

 mybc=bcutil.bcstar('',5300,3.89,-0.23,0.04,Rv=2.5,filters='G3 RP3 K2M',outfile='mystar.csv',retarr=True)

Luca Casagrande Feb 2022
"""

import glob
import numpy as np
from astropy.table import Table 
from scipy.interpolate import Akima1DInterpolator, interp1d

# bracket by +/-nn values over (irregular) grid. If idx True, then indices 
# are returned instead
def bracket(inval,grval,nn,idx=False):
    
    norep = np.sort(np.array(list(dict.fromkeys(list(grval)))))
    
    x1    = np.where(norep<=inval)
    x2    = np.where(norep>inval)
    
    if idx==False:
        lo = norep[x1][-nn::]
        up = norep[x2][0:nn]        
    else:
        lo = x1[0][-nn::]
        up = x2[0][0:nn]
        
    return(lo,up)

# linear interpolation for 2 points, Akima for more. Returns nan if 
# not possible or if extrapolated. The MARCS grid of BC used here is ordered
# such that gridt is monotonic. If not, sorting is necessary.
def mal(val,gridt,gridbc,dset):
    if len(dset[0])>2:
        mfun = Akima1DInterpolator(gridt[dset],gridbc[dset])
        itp  = mfun(val)
    if len(dset[0])==2:
        mfun = interp1d(gridt[dset],gridbc[dset],bounds_error=False) 
        itp  = mfun(val)        
    if len(dset[0])<2:
        itp = np.nan
    return(itp)

# compute Bolometric Corrections for stars of known input parameters
def bcstar(sid,teff,logg,feh,ebv,Rv=False,filters=False,outfile=False,retarr=False):

    # check input data are OK
    if np.isscalar(teff)==True:
        sid  = np.atleast_1d(sid)
        teff = np.atleast_1d(teff)
        logg = np.atleast_1d(logg)
        feh  = np.atleast_1d(feh)
        ebv  = np.atleast_1d(ebv)
    else:
        biglist = [sid,teff,logg,feh,ebv]
        it      = iter(biglist)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            print('** ERROR ** Not all required inputs have same length. Exiting now ...')
            raise SystemExit
            #raise ValueError('Not all required inputs have same length!')
            
        check_arr = isinstance(teff,np.ndarray)+isinstance(logg,np.ndarray)+isinstance(feh,np.ndarray)+isinstance(ebv,np.ndarray)
        if check_arr < 4:
            print('** ERROR **: teff,logg,feh,ebv all need to be arrays. Exiting now ...')
            raise SystemExit
            
    # Rv=3.1 by default, or Rv=2.5. Other values not allowed
    if Rv==False or Rv==3.1:
        el = '3.1'
    elif Rv==2.5:
        el = '2.5'
    else:
        el=input('Wrong Rv, enter either 3.1 or 2.5: ')
    if el!='2.5' and el!='3.1':
        print('Wrong Rv, exiting now ...')
        raise SystemExit

    # keep only requested filters. By default only Gaia DR3 G,BP,RP
    flist  = ['','','','','','','','','']
    rmi    = []
    frange = np.arange(9)
    if type(filters)==str:
        if 'G2'  in filters: flist[0]='BC_G2'
        else: rmi.append(0)

        if 'BP2' in filters: flist[1]='BC_BP2'
        else: rmi.append(1)
        
        if 'RP2' in filters: flist[2]='BC_RP2'
        else: rmi.append(2)
        
        if 'G3'  in filters: flist[3]='BC_G3'
        else: rmi.append(3)
        
        if 'BP3' in filters: flist[4]='BC_BP3'
        else: rmi.append(4)
        
        if 'RP3' in filters: flist[5]='BC_RP3'
        else: rmi.append(5)

        if 'J2M' in filters: flist[6]='BC_J'
        else: rmi.append(6)
        
        if 'H2M' in filters: flist[7]='BC_H'
        else: rmi.append(7)
        
        if 'K2M' in filters: flist[8]='BC_Ks'
        else: rmi.append(8)

        frange = np.delete(frange,rmi)
        flist  = list(np.delete(flist,rmi))

        if len(rmi)==9:
            print('Not suitable filters have been selected.')
            print('Try again by choosing among the following ones:')
            print('G2 BP2 RP2 G3 BP3 RP3 J2M H2M K2M. Exiting now ...')
            raise SystemExit
    else:
        frange = np.arange(3,6) 
        flist  = ['BC_G3','BC_BP3','BC_RP3']
        rmi    = [0,1,2,6,7,8]                 

    # write output file. By default bctable_Rv?.?.dat    
    if type(outfile)==str:
        ouf=open(outfile,'w')
    else: ouf = open('bctable_Rv'+el+'.csv','w')

    header='starID,Teff,feh,logg,ebv,'+','.join(flist)
    ouf.write(header+'\n')

    # read input tables of BCs for several values of E(B-V)
    files  = np.sort(glob.glob('grid/STcol*Rv'+el+'*.dat'))
    gebv   = []
    gri_bc = []
    
    kk=0
    for f in files:

        gebv.append(float(f[-8:-4]))
    
        grid = Table.read(f,format='ascii')
        if kk==0:
            gteff, gfeh, glogg = grid['Teff'],grid['feh'],grid['logg']
        
        bc_g2  = grid['mbol']-grid['G2']
        bc_bp2 = grid['mbol']-grid['BP2']
        bc_rp2 = grid['mbol']-grid['RP2']

        bc_g3  = grid['mbol']-grid['G3']
        bc_bp3 = grid['mbol']-grid['BP3']
        bc_rp3 = grid['mbol']-grid['RP3']

        bc_j   = grid['mbol']-grid['J']
        bc_h   = grid['mbol']-grid['H']
        bc_k   = grid['mbol']-grid['Ks']
    
        tmp = np.transpose([bc_g2,bc_bp2,bc_rp2,bc_g3,bc_bp3,bc_rp3,bc_j,bc_h,bc_k])
        gri_bc.append(tmp)
    
        kk=kk+1
    
    gebv   = np.array(gebv)
    gri_bc = np.array(gri_bc)

    itp_bc = np.zeros(9) + np.nan
    if retarr==True: arr_bc  = np.zeros((len(teff),9)) + np.nan

    # remove entries with E(B-V) outside grid or undetermined. Also remove a
    # few points towards the edge of the grid, where differences wrt Fortran
    # interpolation routines in BCcodes.tar are the largest.    
    rme       = np.where((ebv>=0.72) | (ebv<0.) | (np.isfinite(ebv)==False) | ((logg<2.5) & (teff<3500)))
    # use this instead not to remove grid edges at low log(g) values
    #rme       = np.where((ebv>=0.72) | (ebv<0.))

    
    fold      = feh.copy()
    bold      = ebv.copy()
    
    fold[rme] = -99.
    bold[rme] =   0.
    
    for i in range(len(teff)):
    
        # take +/-3 steps in [Fe/H] grid
        snip = np.concatenate(bracket(fold[i],gfeh,3))
        itp1 = np.zeros((2,len(snip),9))+np.nan

        # take +/-1 step in E(B-V) grid
        eb   = np.concatenate(bracket(bold[i],gebv,1,idx=True))
        
        for k in range(len(snip)):
        
            x0   = np.where((gfeh==snip[k]) & (np.abs(glogg-logg[i])<1.1))
            lg0  = np.array(list(dict.fromkeys(list(glogg[x0]))))
            itp0 = np.zeros((2,len(lg0),9))+np.nan
        
            # at given logg and feh, range of Teff to interpolate across
            for j in range(len(lg0)):
                ok      = np.where((np.abs(gteff-teff[i])<1000) & \
                                   (gfeh==snip[k]) & (glogg==lg0[j]))
                # do it for all selected filters
                for f in frange:
                    itp0[0,j,f] = mal(teff[i],gteff,gri_bc[eb[0],:,f],ok)
                    itp0[1,j,f] = mal(teff[i],gteff,gri_bc[eb[1],:,f],ok)
                
            for f in frange:
                # remove any nan, in case. Either of itp[?,:,:] is enough
                k0 = np.where(np.isnan(itp0[0,:,f])==False)
                # interpolate in logg at correct Teff
                itp1[0,k,f] = mal(logg[i],lg0,itp0[0,:,f],k0)
                itp1[1,k,f] = mal(logg[i],lg0,itp0[1,:,f],k0)

        for f in frange:
            # remove any nan, in case
            k1  = np.where(np.isnan(itp1[0,:,f])==False)
            lor = mal(fold[i],snip,itp1[0,:,f],k1)
            upr = mal(fold[i],snip,itp1[1,:,f],k1)

            # linear interpolate in reddening
            itp_bc[f] = lor + (upr-lor)*(bold[i]-gebv[eb][0])/(gebv[eb][1]-gebv[eb][0])
            if retarr==True: arr_bc[i,f] = itp_bc[f]
            
        ouf.write(','.join([str(sid[i]),str(teff[i]),str(feh[i]),str(logg[i]),str(ebv[i])] + [str(format(itp_bc[f],'.3f')) for f in frange]))
        ouf.write('\n')
            
    ouf.close()

    if retarr==True:
        bctab = np.delete(arr_bc,rmi,1) 
        print('The columns of the output array are '+' '.join(flist))
        return(bctab)
    else:
        return(None)

    del fold
    del bold
    
    print(' ')
    print('These Bolometric Corrections have been computed using:')
    print('- MARCS models with standard [Fe/H] vs [alpha/Fe] enrichment')
    print("- Cardelli/O'Donnel extinction law with Rv="+el)
    print('- MBol_Sun=4.75.')
    print('To derive bolometric flux on the Earth, use Eq (3) of Casagrande & VandenBerg (2018, MNRAS, 475, 5023).\n')
