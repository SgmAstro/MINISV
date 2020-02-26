# configure the desi master environment
source /global/cfs/cdirs/desi/software/desi_environment.sh master

# define some path shortcuts
export REDUXDIR=/global/cfs/cdirs/desi/spectro/redux/daily/
export OUTDIR=$CSCRATCH/BGS/MINISV/rrgen

mkdir -p $OUTDIR

# redrock currently doesn't directly read cframes so reformat into coadds
# even if only a single exposure (<1 min to run)
time desi_coadd_spectra --coadd-cameras \
     -i /global/homes/m/mjwilson/BGS/MINISV/rrgen/70502/cframe-[brz]0-*.fits \
     -o $OUTDIR/coadd-0-00052112.fits

# run redrock on an interactive node
# (~3.5 min on one node for a coadd of one spectrograph)
#
srun -N 1 -n 32 -c 2 -C haswell -t 00:10:00 -q interactive rrdesi_mpi -o $OUTDIR/redrock-0-00051053.h5 -z $OUTDIR/zbest-0-00051053.fits $OUTDIR/coadd-0-00052112.fits
