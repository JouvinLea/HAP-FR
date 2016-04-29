#! /usr/bin/env python
import time
import sys
import logging
import subprocess
import os
import shutil
import hashlib
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import join as table_join
from astropy.time import Time
from glob import glob


"""
Script to check the runwise fits file
You have to give the runlist you want to check, the dst production and the analysis name
./checkfits.py "GC_Lt2deg_PA.list" "Prod15_4_stereo" "elm_north_stereo_Prod15_5"
"""

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
    

class ListObservations:
    def __init__(self, runlist_file, config):
        self.observations = []
        runlist= np.loadtxt(runlist_file)
        obs_ids=runlist[:,0].astype(int)
        telcodes=runlist[:,1].astype(int)
        for obs_id, telcode in zip(obs_ids, telcodes):
            obs = Observation(obs_id, config, telcode)
            self.observations.append(obs)


def make_checkrun(list_observations, indir, informat, outfile):
    """Check corrupted fits file created during the fits production
    """
   
    missfiles=open(outfile,"w")
    irunbad=0
    irungood=0
    for obs in list_observations.observations:
        events_filename = indir+"/"+obs.filename('events', format=informat)
        try:
            table = Table.read(str(events_filename), hdu='EVENTS')
        except Exception:
            print "fits corrupted for file "+str(events_filename)
            missfiles.write(str(obs.obs_id))
            irunbad=irunbad+1
            continue
        irungood= irungood+1
    missfiles.close()
    print("There are "+str(irunbad)+" corrupted files over "+str(len(list_observations.observations)))
    print("There are "+str(irungood)+" good files over "+str(len(list_observations.observations)))
    


if __name__ == '__main__':
    runlist = sys.argv[1]
    dstprod=sys.argv[2]
    analysis_name=sys.argv[3]
    observation_list = ListObservations(runlist ,analysis_name)
    indir= os.path.expandvars('$CALDB')+"/data/hess/"+os.path.expandvars('$HESSVERSION')+"/"+dstprod+'/'+analysis_name
    make_checkrun(observation_list, indir, "old","missig_file_"+analysis_name+".txt")
    
