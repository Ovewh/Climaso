import xarray as xr
from pyclim_noresm.aerosol_forcing import calc_cloud_forcing
import copy
import time

import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )

VARS = snakemake.config['variables']

ERFaf = xr.open_dataset(snakemake.input.ERFaf)
ERFafcs = xr.open_dataset(snakemake.input.ERFafcs)

cloudEffect = calc_cloud_forcing(ERFaf[VARS[snakemake.wildcards.vName][0]],
                                    ERFafcs[VARS[snakemake.wildcards.vName][1]])

cloudEffect = cloudEffect.to_dataset(name=snakemake.wildcards.vName)
cloudEffect.attrs=copy.copy(ERFaf.attrs)
cloudEffect.attrs['title']='Aerosol direct radiative effect'
cloudEffect.attrs['history'] = f'@{time.ctime()} Generated by: {__file__} ' + cloudEffect.attrs['history']
cloudEffect.attrs['source'] = ', '.join(snakemake.input) + ', ' + cloudEffect.attrs['source']
cloudEffect.to_netcdf(snakemake.output.outpath)