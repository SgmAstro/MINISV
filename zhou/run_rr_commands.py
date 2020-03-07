import os


run       =  True
reverse   = False
overwrite = False

ctype     = 'rr'  ## ['rr', 'coadd']

# commands = open('commands_{}_allexp.sh'.format(ctype), 'r')
commands   = open('commands_{}_nexp_1.sh'.format(ctype), 'r')

outdir     = os.environ['OUTDIR']
 
lines      = commands.readlines()

if reverse:
  lines    = reversed(lines)

print('\n\n')

for x in lines:
  fpath = x.split('-o')[-1].split('-z')[0].strip().replace('$OUTDIR', outdir)

  print(fpath, os.path.exists(fpath))

  if (ctype == 'rr') & (not os.path.exists(x.split(' ')[-1].replace('$OUTDIR', outdir).replace('\n',''))):
    cpath = x.split(' ')[-1].replace('$OUTDIR', outdir).replace('\n','')

    print(cpath, os.path.exists(cpath))
    
    continue 
    
  if run & ((not os.path.exists(fpath)) | overwrite):
    os.system(x)
        
print('\n\nDone.\n\n')
