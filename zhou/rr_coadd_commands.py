# Get the a list of science exposures for a specficic date and tile
# Print desi_coadd_spectra commands for coadds of single or multiple exposures
# Print redrock commands

from    __future__ import division, print_function

import  ephem
import  fitsio
import  desisurvey.config
import  sys, os, glob, time, warnings
import  numpy as np
import  astropy.units as u
import  desisurvey.utils   as      dutils

from    astropy.time  import Time
from    astropy.table import Table, vstack, hstack


redux_dir     = '/global/cfs/cdirs/desi/spectro/redux/daily'       
output_dir    = '/global/cscratch1/sd/mjwilson/BGS/MINISV/coadds/'

nights        = ['20200225', '20200227', '20200228', '20200229', '20200303']

# number of exposures in a coadded; 1 for single-exposure coadd
ALL           = True # Overrides n_exp.
n_exp         = 1
n_node        = 4

overwrite     = True

# tileid_list = None  # no restriction on tiles
tileid_list   = [70500, 70502, 70510]

petals        = range(10)
# petals      = [0, 3, 6, 7, 9]


################################# Setuop ##################################################
config           = desisurvey.config.Configuration()

mayall           = ephem.Observer()
mayall.lat       = config.location.latitude().to(u.rad).value
mayall.lon       = config.location.longitude().to(u.rad).value
mayall.elevation = config.location.elevation().to(u.m).value

################################## Get list of exposures ##################################
exposure_dir_list = []

for obsdate in nights:
  exposure_dir_list += glob.glob(os.path.join(redux_dir, 'exposures', obsdate, '*'))

cframe_list = []

# Get a list of all science exposures.
for exposure_dir in exposure_dir_list:
    cframe_list_tmp = glob.glob(os.path.join(exposure_dir, 'cframe-*'))

    if len(cframe_list_tmp) > 0:
        if tileid_list is None:
            cframe_list += cframe_list_tmp

        else:
            # only need to check one cframe file in the exposure
            with fitsio.FITS(cframe_list_tmp[0]) as f:
                if f[0].read_header()['TILEID'] in tileid_list:
                    cframe_list += cframe_list_tmp

                    
cframe_list           = sorted(cframe_list)

# Gather exposure/petal information
cframes               = Table()

cframes['cframe']     = np.array(cframe_list)
cframes['night']      = '20200227'
cframes['mjd']        = np.zeros(len(cframes), dtype=np.float) - 1.0 
cframes['lat']        = np.zeros(len(cframes), dtype=np.float)
cframes['lon']        = np.zeros(len(cframes), dtype=np.float)
cframes['elv']        = np.zeros(len(cframes), dtype=np.float)
cframes['tileid']     = np.zeros(len(cframes), dtype=int)
cframes['expid']      = np.zeros(len(cframes), dtype=int)
cframes['exptime']    = np.zeros(len(cframes), dtype=np.float)
cframes['camera']     = ' '
cframes['petal_loc']  = -1 * np.ones(len(cframes), dtype=np.int32)
cframes['spectrograph'] = -1 * np.ones(len(cframes), dtype=np.int32)
cframes['ra']         = np.zeros(len(cframes), dtype=np.float)
cframes['dec']        = np.zeros(len(cframes), dtype=np.float)

##  Celestial.
cframes['MOONALT']   = np.zeros(len(cframes))
cframes['MOONRA']    = np.zeros(len(cframes))
cframes['MOONDEC']   = np.zeros(len(cframes))
cframes['MOONFRAC']  = np.zeros(len(cframes))
cframes['MOONSEP']   = np.zeros(len(cframes))

cframes['SUNALT']    = np.zeros(len(cframes))
cframes['SUNRA']     = np.zeros(len(cframes))
cframes['SUNDEC']    = np.zeros(len(cframes))
cframes['SUNSEP']    = np.zeros(len(cframes))

for index, cframe in enumerate(cframes['cframe']):
    with fitsio.FITS(cframe) as f:
        header                         = f[0].read_header()

        # print(header)
        
        cframes['mjd'][index]          = header['MJD-OBS']
        cframes['night'][index]        = header['NIGHT']
        cframes['tileid'][index]       = header['TILEID']
        cframes['expid'][index]        = header['EXPID']
        cframes['camera'][index]       = header['CAMERA'].strip()[0]
        # cframes['petal_loc'][index]  = int(header['CAMERA'].strip()[1])
        cframes['spectrograph'][index] = int(header['CAMERA'].strip()[1])
        
        cframes['lat'][index]        = header['OBS-LAT']
        cframes['lon'][index] 	     = header['OBS-LONG']
        cframes['elv'][index] 	     = header['OBS-ELEV']
        cframes['exptime'][index]    = header['EXPTIME']
        cframes['ra'][index]         = header['SKYRA']
        cframes['dec'][index]        = header['SKYDEC']

        ##  Celestial.
        mayall.date                 = Time(cframes['mjd'][index], format='mjd').datetime

        _moon                       = ephem.Moon()
        _moon.compute(mayall)

        _sun                        = ephem.Sun()
        _sun.compute(mayall)

        cframes['MOONALT'][index]   = 180./np.pi*_moon.alt
        cframes['MOONRA'][index]    = 180./np.pi*_moon.ra
        cframes['MOONDEC'][index]   = 180./np.pi*_moon.dec
        cframes['MOONFRAC'][index]  = _moon.moon_phase

        cframes['SUNALT'][index]    = 180./np.pi*_sun.alt
        cframes['SUNRA'][index]     = 180./np.pi*_sun.ra
        cframes['SUNDEC'][index]    = 180./np.pi*_sun.dec

        cframes['MOONSEP']          = np.diag(dutils.separation_matrix(cframes['MOONRA'] * u.deg, cframes['MOONDEC'].quantity * u.deg, cframes['ra'] * u.deg, cframes['dec'] * u.deg))
        cframes['SUNSEP']           = np.diag(dutils.separation_matrix(cframes['SUNRA'] * u.deg, cframes['SUNDEC'].quantity * u.deg, cframes['ra'] * u.deg, cframes['dec'] * u.deg))
        
# Sanity check: each petal must have three cframe files
for expid in np.unique(cframes['expid']):
    mask_expid = cframes['expid']==expid
    
    for petal_loc in range(10):
        mask = mask_expid & (cframes['petal_loc']==petal_loc)
        if (np.sum(mask)>0) & (np.sum(mask)!=3):
            raise ValueError('EXPID {} PETAL_LOC {} has only {} cframes files'.format(expid, petal_loc, np.sum(mask)))

print('\n\n')
        
cframes.sort(('tileid', 'petal_loc'))

cframes.pprint(max_lines=-1, max_width=-1)

cframes.write('bgs_allcframes.fits', format='fits', overwrite=True)

exit(0)

## 
output_argument_list = []

if not ALL:
  if (not overwrite) and (os.path.isfile("commands_coadd_nexp_{}.sh".format(n_exp)) | os.path.isfile("commands_rr_nexp_{}.sh".format(n_exp))):
    raise  ValueError('Overwrite=True required to remove exisiting files.')

  output_file   = open("commands_coadd_nexp_{}.sh".format(n_exp), "w")
  output_rrfile = open("commands_rr_nexp_{}.sh".format(n_exp),    "w")
  
else:
  if (not overwrite) and (os.path.isfile("commands_coadd_allexp.sh") | os.path.isfile("commands_rr_allexp.sh")):
    raise  ValueError('Overwrite=True required to remove exisiting files.')
  
  output_file   = open("commands_coadd_allexp.sh", "w")
  output_rrfile = open("commands_rr_allexp.sh",    "w")

##
for tileid in np.unique(cframes['tileid']):
  for night in nights:
    for petal_loc in petals: 
        mask                    = (cframes['tileid'] == tileid) & (cframes['petal_loc'] == petal_loc) & (cframes['night'] == night)

        # Choose one camera for simplicity 
        mask                   &= (cframes['camera'] == 'b')

        cframe1                 = cframes[mask]

        if not ALL:
          if (np.sum(mask) < n_exp):
            print('\n# {} exposures is not enough for TILEID {}, NIGHT {} PETAL_LOC {} for a {} coadd.\n'.format(np.sum(mask), tileid, night, petal_loc, n_exp))
            continue

          # Skip the exposures that do not make the split.
          cframe1               = cframe1[:len(cframe1) - len(cframe1) % (n_exp)]

          nsplit                = len(cframe1)//(n_exp)

          subset_split          = np.split(np.arange(len(cframe1)), nsplit)

        else:
          if (np.sum(mask) == 0):
            print('\n# No exposures for TILEID {}, PETAL_LOC {}.\n'.format(tileid, petal_loc))
            continue
          
          nsplit                = 1
          subset_split          = np.split(np.arange(len(cframe1)), nsplit)
          
        for subset_index in range(len(subset_split)):
            subset              = cframe1[subset_split[subset_index]]
            input_argument      = ''

            for index in range(len(subset)):
                exposure_dir    = os.path.dirname(subset['cframe'][index])
                input_argument += os.path.join(exposure_dir.replace(redux_dir, '$REDUXDIR'), 'cframe-[brz]{}-*.fits ').format(petal_loc)

            if n_exp == 1:
                exposure        = os.path.basename(exposure_dir)
                output_argument = os.path.join('$OUTDIR', 'NEXP{}'.format(n_exp), str(tileid), night, 'coadd-{}-{}-{}.fits'.format(night, petal_loc, exposure))

            elif not ALL:
                output_argument = os.path.join('$OUTDIR', 'NEXP{}'.format(n_exp), str(tileid), night, 'coadd-{}-{}-{}exp-subset-{}.fits'.format(night, petal_loc, n_exp, subset_index))

            else:
                output_argument = os.path.join('$OUTDIR', 'ALL', str(tileid), night, 'coadd-{}-{}-allexp.fits'.format(night, petal_loc))
                
            output_argument_list.append(output_argument)

            if os.path.isfile(output_argument) and (not overwrite):
                print('\nWarninig: {} already exists!\n'.format(output_argument))
                continue

            output_file.write('time desi_coadd_spectra --coadd-cameras -i {} -o {}\n'.format(input_argument, output_argument))

##
output_file.close()
        
for output_argument in output_argument_list:
    rrdesi_argument_redrock = output_argument.replace('coadd', 'redrock').replace('.fits', '.h5')
    rrdesi_argument_zbest   = output_argument.replace('coadd', 'zbest')

    output_rrfile.write('srun -N {} -n {} -c 2 rrdesi_mpi -o {} -z {} {}\n'.format(n_node, 32 * n_node, rrdesi_argument_redrock, rrdesi_argument_zbest, output_argument))

print('\n\nDone.\n\n')
