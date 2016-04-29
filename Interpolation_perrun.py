#! /usr/bin/env python
import numpy as np
from astropy.io import fits
#import pyfits
from astropy.table import Table
from scipy import interpolate
import math
from astropy.io.fits import Column
#from pyfits import Column
import sys
import os
from glob import glob

"""
Script that interpolates the area, edisp and the psf on the zenith and muon efficiency of each run.
You have to give the analysis name, the run number and the dst prod
"""
#./Interpolation_perrun.py 'elm_north_stereo_Prod15_5' 23526 "Prod15_4_stereo"



class Observation:
    """Helper functions to compute file and folder names.
    """

    #filetypes = ['events', 'aeff', 'edisp', 'psf_3gauss']
    filetypes = ['events']

    def __init__(self, obs_id, hap_config=None, telpattern=None):
        self.obs_id = obs_id
        self.hap_config = hap_config
        self.telpattern = telpattern

    @property
    def obs_group(self):
        obs_id_min = self.obs_id - (self.obs_id % 200)
        obs_id_max = obs_id_min + 199
        return obs_id_min, obs_id_max

    @property
    def _obs_group_folder(self):
        return 'run{:06d}-{:06d}'.format(self.obs_group[0], self.obs_group[1])

    @property
    def _obs_folder(self):
        return 'run{:06d}'.format(self.obs_id)

    def folder(self, step=None):
        """Create folder for a given step.
        """
        if step is None:
            return self._obs_group_folder+"/"+self._obs_folder
        else:
            return step+"/"+self._obs_group_folder+"/"+self._obs_folder

    def filename(self, filetype, format='old'):
        if format == 'old':
            TAGS = dict(
                events='events',
                aeff='aeff_2d',
                edisp='edisp_2d',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='bkg_offruns',
            )
        elif format == 'new':
            TAGS = dict(
                events='events',
                aeff='aeff',
                edisp='edisp',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='background',
            )

        tag = TAGS[filetype]
        if(filetype=="events"):
            filename = '{}_{:06d}.fits.gz'.format(tag, self.obs_id)
        else:
            filename = '{}_{:06d}.fits'.format(tag, self.obs_id)
        return self.folder()+"/"+filename
    
    
def gauss(x,sigma, mean):
    f=1/(np.sqrt(2*math.pi)*sigma)*np.exp(-(x-mean)**2/(2*sigma**2))
    return f





#Load the MCs information from the IRF table where is stored for each MC energy, zenith, offset and efficiency, the value of the surface area, the biais and the sigma for the resolution and the s1,s2,s3,A2, A3 for the tripple gauss used to fit the PSF
PathTableIRF=os.path.expandvars('$HESSCONFIG')
PathTablePSF=os.path.expandvars('$HESSCONFIG')


analysis_name=sys.argv[1]
nrun=sys.argv[2]
dstprod=sys.argv[3]
PathListRun = os.path.expandvars('$CALDB')+"/data/hess/"+os.path.expandvars('$HESSVERSION')+"/"+dstprod+'/'+analysis_name

obs = Observation(int(nrun))
informat="old"
namerun=PathListRun +"/"+ obs.filename('events', format=informat)
try:
    table = Table.read(namerun, hdu='EVENTS')
except Exception:
    print "fits corrupted for file "+namerun
else:
    hdurun=fits.open(namerun)
    AltRun=hdurun[1].header["ALT_PNT"]
    ZenRun=90-AltRun
    EffRun=hdurun[1].header["MUONEFF"]*100
    IRF=np.load(PathTableIRF+"/"+analysis_name+"/IRF_"+analysis_name+".npz")
    IRFArea=IRF["TableArea"]
    IRFSigma=IRF["TableSigma"]
    IRFBiais=IRF["TableBiais"]
    enMC=IRF["enMC"]
    lnenMC=IRF["lnenMC"]
    zenMC=IRF["zenMC"]
    effMC=IRF["effMC"]
    offMC=IRF["offMC"]

    PSF=np.load(PathTablePSF+"/"+analysis_name+"/PSF_triplegauss_"+analysis_name+".npz")
    PSFs1=PSF["TableSigma1"]
    PSFs2=PSF["TableSigma2"]
    PSFs3=PSF["TableSigma3"]
    PSFA2=PSF["TableA2"]
    PSFA3=PSF["TableA3"]


    binoffMC=len(offMC)
    binEMC=len(enMC)
    binEreco=50
    bineffarea=len(offMC)*len(enMC)
    bineffresol=len(offMC)*len(enMC)*binEreco

    off_low=offMC
    off_hi=offMC

    #Define Etrue low and up in log
    binlnEMC=lnenMC[1:]-lnenMC[:-1]
    #For the first bin in order to define the low edge we take the width of the first bin
    binlnEMClow=np.insert(binlnEMC,0,binlnEMC[0])
    #For the last bin in order to define the up edge we take the width of the last bin
    binlnEMCup=np.insert(binlnEMC,-1,binlnEMC[-1])
    lnEMClow=lnenMC-binlnEMClow/2
    lnEMCup=lnenMC+binlnEMCup/2
    E_true_low=pow(10,lnEMClow)
    E_true_up=pow(10,lnEMCup)

    #Definition de Etrue/Ereco
    lnEtrue_reco=np.linspace(-1,1,binEreco)
    #For each bin in log the width is the same so we take the width f the first bi
    binlnEtrue_reco=lnEtrue_reco[1]-lnEtrue_reco[0]
    lnE_true_reco_low=lnEtrue_reco-binlnEtrue_reco/2
    lnE_true_reco_up=lnEtrue_reco+binlnEtrue_reco/2
    Etrue_reco=pow(10,lnEtrue_reco)
    E_true_reco_low=pow(10,lnE_true_reco_low)
    E_true_reco_hi=pow(10,lnE_true_reco_up)

    AreaRun=np.zeros((binoffMC,binEMC))
    ResolRun=np.zeros((binoffMC,binEreco,binEMC))

    PSFS1Run=np.zeros((binoffMC,binEMC))
    PSFS2Run=np.zeros((binoffMC,binEMC))
    PSFS3Run=np.zeros((binoffMC,binEMC))
    PSFA2Run=np.zeros((binoffMC,binEMC))
    PSFA3Run=np.zeros((binoffMC,binEMC))

    for (iEMC,EMC) in enumerate(enMC):
        for (ioff, off) in enumerate(offMC):
            InterArea=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFArea[iEMC,ioff,:,:])
            InterBiais=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFBiais[iEMC, ioff,:,:])
            InterSigma=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFSigma[iEMC, ioff,:,:])
            AreaRun[ioff,iEMC]=InterArea(EffRun,np.cos(ZenRun*math.pi/180))
            BiaisRun=InterBiais(EffRun,np.cos(ZenRun*math.pi/180))
            SigmaRun=InterSigma(EffRun,np.cos(ZenRun*math.pi/180))
            ResolRun[ioff, : ,iEMC]=gauss(lnEtrue_reco,SigmaRun,BiaisRun)
            #Resolution is normalized
            norm=np.sum(ResolRun[ioff, : ,iEMC]*(E_true_reco_hi-E_true_reco_low))            

            if(np.isnan(norm)):
                ResolRun[ioff, : ,iEMC]=0
            else:
                ResolRun[ioff, : ,iEMC]=ResolRun[ioff, : ,iEMC]/norm



            ind_zen, ind_eff= np.where(PSFs1[iEMC, ioff, :, :] != -1)
            #If there is at least one simu for this offset and this energy for wich the fit works
            if(len(ind_zen)!=0):
                zensame=np.where(ind_zen != ind_zen[0])
                effsame=np.where(ind_eff != ind_eff[0])
                #In order for the interpolation to work, you have to have at least two values in efficiency and zenith
                if((len(zensame[0])!=0) & (len(effsame[0])!=0)):
                    coord_eff=effMC[ind_eff]
                    coord_zen = zenMC[ind_zen]
                    points= (coord_eff, np.cos(coord_zen * math.pi / 180))

                    PSFS1Run[ioff, iEMC] = interpolate.griddata(points, PSFs1[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='linear')
                    if np.isnan(PSFS1Run[ioff, iEMC]):
                        PSFS1Run[ioff, iEMC] = interpolate.griddata(points, PSFs1[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='nearest')

                    PSFS2Run[ioff, iEMC] = interpolate.griddata(points, PSFs2[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='linear')
                    if np.isnan(PSFS2Run[ioff, iEMC]):
                        PSFS2Run[ioff, iEMC] = interpolate.griddata(points, PSFs2[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='nearest')

                    PSFS3Run[ioff, iEMC] = interpolate.griddata(points, PSFs3[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='linear')
                    if np.isnan(PSFS3Run[ioff, iEMC]):
                        PSFS3Run[ioff, iEMC] = interpolate.griddata(points, PSFs3[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='nearest')

                    PSFA2Run[ioff, iEMC] = interpolate.griddata(points, PSFA2[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='linear')
                    if np.isnan(PSFA2Run[ioff, iEMC]):
                        PSFA2Run[ioff, iEMC] = interpolate.griddata(points, PSFA2[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='nearest')

                    PSFA3Run[ioff, iEMC] = interpolate.griddata(points, PSFA3[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='linear')
                    if np.isnan(PSFA3Run[ioff, iEMC]):
                        PSFS1Run[ioff, iEMC] = interpolate.griddata(points, PSFA3[iEMC, ioff, ind_zen, ind_eff], (EffRun,np.cos(ZenRun * math.pi / 180)), method='nearest')

                else:
                    PSFS1Run[ioff, iEMC] = -1
                    PSFS2Run[ioff, iEMC] = -1
                    PSFS3Run[ioff, iEMC] = -1
                    PSFA2Run[ioff, iEMC] = -1
                    PSFA3Run[ioff, iEMC] = -1
            else:
                PSFS1Run[ioff, iEMC] = -1
                PSFS2Run[ioff, iEMC] = -1
                PSFS3Run[ioff, iEMC] = -1
                PSFA2Run[ioff, iEMC] = -1
                PSFA3Run[ioff, iEMC] = -1


    """
    Writing the fits file for aeff, edisp and psf for the given observation
    """
    outdir =  PathListRun+"/"+obs.folder()
    
    #AEFF FITS FILE
    c1_area = Column(name='ENERG_LO', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_area = Column(name='ENERG_HI', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_area = Column(name='THETA_LO', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_low))
    c4_area = Column(name='THETA_HI', format=str(binoffMC)+'E', unit='def', array=np.atleast_2d(off_hi))
    c5_area = Column(name='EFFAREA', format=str(bineffarea)+'E', unit='TeV', array=np.expand_dims(AreaRun,0))
    c6_area = Column(name='EFFAREA_RECO', format=str(bineffarea)+'E', unit='TeV', array=np.expand_dims(AreaRun,0))
    tbhdu_area = fits.BinTableHDU.from_columns([c1_area,c2_area,c3_area,c4_area,c5_area,c6_area])
    for i in range(1,7):
        tbhdu_area.header.comments['TTYPE'+str(i)]='label for field '+str(i)
        tbhdu_area.header.comments['TFORM'+str(i)]='data format of field: 4-byte REAL'
        tbhdu_area.header.comments['TUNIT'+str(i)]='physical unit of field '

    tbhdu_area.header.set("EXTNAME","EFFECTIVE AREA", "name of this binary table extension ")
    tbhdu_area.header.set("TDIM5","("+str(binEMC)+","+str(binoffMC)+")")
    tbhdu_area.header.set("TDIM6","("+str(binEMC)+","+str(binoffMC)+")")
    tbhdu_area.header.set("LO_THRES",-1,"TeV")
    tbhdu_area.header.set("HI_THRES",-1,"TeV")
    #tbhdu_area.header["EXTNAME"]='EFFECTIVE AREA'
    tbhdu_area.writeto(outdir + '/aeff_2d_0'+str(int(nrun))+'.fits', clobber=True)
    #EDISP FITS FILE
    c1_resol = Column(name='ETRUE_LO', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_resol = Column(name='ETRUE_HI', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_resol = Column(name='MIGRA_LO', format=str(binEreco)+'E', unit='', array=np.atleast_2d(E_true_reco_low))
    c4_resol = Column(name='MIGRA_HI', format=str(binEreco)+'E', unit='', array=np.atleast_2d(E_true_reco_hi))
    c5_resol = Column(name='THETA_LO', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_low))
    c6_resol = Column(name='THETA_HI', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_hi))
    c7_resol = Column(name='MATRIX ', format=str(bineffresol)+'E', unit='TeV', array=np.expand_dims(ResolRun,0))
    tbhdu_resol = fits.BinTableHDU.from_columns([c1_resol,c2_resol,c3_resol,c4_resol,c5_resol,c6_resol,c7_resol])
    for i in range(1,8):
        tbhdu_resol.header.comments['TTYPE'+str(i)]='label for field '+str(i)
        tbhdu_resol.header.comments['TFORM'+str(i)]='data format of field: 4-byte REAL'
    tbhdu_resol.header.set("EXTNAME","EDISP_2D", "name of this binary table extension ")
    tbhdu_resol.header.set("TDIM7","("+str(binEMC)+","+str(binEreco)+","+str(binoffMC)+")")
    tbhdu_resol.writeto(outdir +'/edisp_2d_0'+str(int(nrun))+'.fits', clobber=True)
    
    #PSF FITS FILE
    c1_psf = Column(name='ENERG_LO', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_psf = Column(name='ENERG_HI', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_psf = Column(name='THETA_LO', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_low))
    c4_psf = Column(name='THETA_HI', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_hi))

    norm=2*np.pi*(PSFS1Run**2+PSFA2Run*PSFS2Run**2+PSFA3Run*PSFS3Run**2)
    c5_psf = Column(name='SIGMA_1', format=str(bineffarea)+'E', unit='deg', array=np.expand_dims(PSFS1Run,0))
    c6_psf = Column(name='AMPL_2', format=str(bineffarea)+'E', unit='', array=np.expand_dims(PSFA2Run,0))
    c7_psf = Column(name='SIGMA_2', format=str(bineffarea)+'E', unit='deg', array=np.expand_dims(PSFS2Run,0))
    c8_psf = Column(name='AMPL_3', format=str(bineffarea)+'E', unit='', array=np.expand_dims(PSFA3Run,0))
    c9_psf = Column(name='SIGMA_3', format=str(bineffarea)+'E', unit='deg', array=np.expand_dims(PSFS3Run,0))
    c10_psf = Column(name='SCALE', format=str(bineffarea)+'E', unit='', array=np.expand_dims(1/norm,0))
    tbhdu_psf = fits.BinTableHDU.from_columns([c1_psf, c2_psf, c3_psf, c4_psf, c5_psf, c6_psf, c7_psf, c8_psf, c9_psf, c10_psf]) 
    tbhdu_psf.header.set("EXTNAME","PSF_2D", "name of this binary table extension ")
    tbhdu_psf.writeto(outdir + '/psf_3gauss_0'+str(int(nrun))+'.fits', clobber=True)
    
