

rule plot_ERFs_historical:
    input:
        paths=expand(outdir+'histSST/ERFs/{vName}/{vName}_histSST_{model}_Ayear.nc', 
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem'], allow_missing=True)
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
                      'EC-Earth3-AerChem'], allow_missing=True)
    
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
                      'EC-Earth3-AerChem'], allow_missing=True)
    
    output:
        outpath=outdir+'figs/forcing_history/ERFs/{vName}_historical_CMIP6.png'

    wildcard_constraints:
        vName = 'SWDirectEff|LWDirectEff|LWDirectEff_cs|SWDirectEff_cs|DirectEff'
    params:
        vmin=-0.5,
        vmax=0.5

    notebook:
        "../notebooks/plot_ERFs_time_series.py.ipynb"

rule plot_historical_timeseries:
    input:
        paths=expand(outdir+'{experiment}/{vName}/{vName}_{experiment}_{model}_Ayear.nc',
                   model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem'], allow_missing=True),
    output:
        outpath=outdir+'figs/forcing_history/{vName}_{experiment}_CMIP6_ts.png'

    wildcard_constraints:
        experiment = 'histSST|historical|histSST-piAer'
    notebook:
        '../notebooks/plot_time_series.py.ipynb' 


rule plot_change_historical:
    input:
        paths_hist=expand(outdir+'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                   model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem'], allow_missing=True),
        paths_piaer_hist=expand(outdir+'histSST-piAer/{vName}/{vName}_histSST-piAer_{model}_Ayear.nc',
                    model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'EC-Earth3-AerChem'], allow_missing=True)


    output:
        outpath=outdir+'figs/forcing_history/{vName}_historical_CMIP6_{kind}_diff.png'        

    notebook:
        "../notebooks/plot_change_time_series.py.ipynb"

rule build_historical_ts_figures:
    input:
        expand(rules.plot_change_historical.output.outpath, 
                vName=['clt','clwvi','clivi','od550aer','od550lt1aer'], kind='abs'),
        expand(rules.plot_historical_timeseries.output.outpath, 
                vName=['clt','clwvi','clivi','od550aer','od550lt1aer'], experiment=['histSST', 'histSST-piAer']),
        expand(rules.plot_ERFs_historical.output.outpath,
                vName=['ERFtsw','ERFtlw','ERFt']),
        expand(rules.plot_cloudEff_historical.output.outpath,
                vName=['SWCloudEff','LWCloudEff','CloudEff']),
        expand(rules.plot_DirectEff_historical.output.outpath,
                vName=['SWDirectEff','LWDirectEff','LWDirectEff_cs','SWDirectEff_cs','DirectEff'])
        

rule convert_to_aerocom:
    input:
        expand(outdir+'{model}_histSST/renamed/converted_CMIP6_aerocom_od550aer.txt'
                ,model=['NorESM2-LM', 'MPI-ESM-1-2-HAM','UKESM1-0-LL',
                      'MIROC6', 'EC-Earth3-AerChem', 'MRI-ESM2-0', 'GISS-E2-1-G']),

        expand(outdir+'{model}_histSST/renamed/converted_CMIP6_aerocom_od550csaer.txt'
                ,model=['NorESM2-LM', 'MPI-ESM-1-2-HAM',
                      'MIROC6', 'GISS-E2-1-G']),
        expand(outdir+'{model}_histSST/renamed/converted_CMIP6_aerocom_od550lt1aer.txt',
                model=['NorESM2-LM', 'EC-Earth3-AerChem']),
        
        expand(outdir+'{model}_historical/renamed/converted_CMIP6_aerocom_od550aer.txt'
                ,model=['NorESM2-LM', 'MPI-ESM-1-2-HAM','EC-Earth3-AerChem']),
        expand(outdir+'{model}_historical/renamed/converted_CMIP6_aerocom_od550lt1aer.txt'
                ,model=['NorESM2-LM'])
    