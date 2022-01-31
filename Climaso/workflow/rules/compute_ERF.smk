from pyclim_noresm.aerosol_forcing import calc_total_ERF_TOA
import xarray as xr
import os



rule calculate_ERF_TOA:
    input:
        exp_downwelling_SW = 
            lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/{w.model}/{w.experiment}/**/Amon/{VAR[w.vName]}/**/latest/*.nc'),
    output:
        outpath = 'results/{{vName}}_{{experiment}}_{{model}}_{{freq}}.nc'
    run:
        print(input.exp_downwelling_SW)

            
