import os


run       = True

ctype     = 'rr'  ## ['rr', 'coadd']

commands  = open('commands_{}_nexp_1.sh'.format(ctype), 'r')
outdir    = os.environ['OUTDIR']

print('\n\n')

for x in commands:
  fpath = x.split('-o')[-1].split('-z')[0].strip().replace('$OUTDIR', outdir)
    
  print(fpath, os.path.exists(fpath))

  if run & (not os.path.exists(fpath)):
    os.system(x)
        
print('\n\nDone.\n\n')
