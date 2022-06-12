

rule plot_ERFs_historical:
    input:
        paths=expand(outdir+'histSST/ERFs/{vName}/{vName}_histSST_{model}_Ayear.nc', 
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM','UKESM1-0-LL',
                      'MIROC6', 'EC-Earth3-AerChem', 'MRI-ESM2-0'], allow_missing=True)
    output:
        outpath=outdir+'figs/forcing_history/ERFs/{vName}_historical_CMIP6.png'
    
    wildcard_constraints:
        vName = 'ERFtsw|ERFtlw|ERFtcs|ERFt'
    params:
        vmin=-2.5,
        vmax=1

    notebook:
        "../notebooks/plot_ERFs_time_series.py.ipynb"

rule plot_cloudEff_historical:
    input:
        paths=expand(outdir+'histSST/ERFs/{vName}/{vName}_histSST_{model}_Ayear.nc', 
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem','UKESM1-0-LL'], allow_missing=True)
    
    output:
        outpath=outdir+'figs/forcing_history/ERFs/{vName}_historical_CMIP6.png'

    wildcard_constraints:
        vName = 'SWCloudEff|LWCloudEff|CloudEff'

    notebook:
        "../notebooks/plot_ERFs_time_series.py.ipynb"
    

rule plot_DirectEff_historical:
    input:
        paths=expand(outdir+'histSST/ERFs/{vName}/{vName}_histSST_{model}_Ayear.nc', 
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem','UKESM1-0-LL'], allow_missing=True)
    
    output:
        outpath=outdir+'figs/forcing_history/ERFs/{vName}_historical_CMIP6.png'

    wildcard_constraints:
        vName = 'SWDirectEff|LWDirectEff|LWDirectEff_cs|SWDirectEff_cs|DirectEff'
    params:
        vmin=-0.5,
        vmax=0.5

    notebook:
        "../notebooks/plot_ERFs_time_series.py.ipynb"

rule convert_od550aer_aerocom:
    input:
        expand(outdir+'{model}_histSST/renamed/converted_CMIP6_aerocom_od550aer.txt'
                ,model=['NorESM2-LM', 'MPI-ESM-1-2-HAM','UKESM1-0-LL',
                      'MIROC6', 'EC-Earth3-AerChem', 'MRI-ESM2-0', 'GISS-E2-1-G'])

rule convert_od550csaer_aerocom:
    input:
        expand(outdir+'{model}_histSST/renamed/converted_CMIP6_aerocom_od550csaer.txt'
                ,model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'MIROC6', 'GISS-E2-1-G'])