import xarray as xr
from pyclim_noresm.aerosol_forcing import calc_cloud_forcing
import copy
import time

VARS = snakemake.config['variables']

ERFaf = xr.open_dataset(snakemake.input.ERFaf)
ERFafcs = xr.open_dataset(snakemake.input.ERFafcs)

cloudEffect = calc_cloud_forcing(ERFaf[VARS[snakemake.wildcards.vName][0]],
                                    ERFafcs[VARS[snakemake.wildcards.vName][1]])

cloudEffect = cloudEffect.to_dataset(name=snakemake.wildcards.vName)
cloudEffect.attrs=copy.copy(ERFaf.attrs)
cloudEffect['title']='Aerosol direct radiative effect'
cloudEffect['history'] = f'@{time.ctime()} Generated by: {__file__} ' + cloudEffect.attrs['history']
cloudEffect['source'] = ', '.join(snakemake.input) + ', ' + cloudEffect.attrs['source']
cloudEffect.to_netcdf(snakemake.output.outpath)