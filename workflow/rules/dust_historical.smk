


rule plot_change_historical_dust:
    input:
        pertubation = expand(outdir+'historical/{vName}/{vName}_historical_{model}_Ayear.nc',
                        model=['NorESM2-LM',
                        'GISS-E2-1-G', 'EC-Earth3-AerChem',
                        'GFDL-ESM4'],allow_missing=True),
        baseline = expand(outdir+'hist-piAer/{vName}/{vName}_hist-piAer_{model}_Ayear.nc',
                        model=['NorESM2-LM', 
                        'GISS-E2-1-G',  'EC-Earth3-AerChem',
                        'GFDL-ESM4'],allow_missing=True),
    
        area_weights = expand('workflow/input_data/gridarea_{model}.nc',model=['NorESM2-LM','CNRM-ESM2-1'])
    output:
        outpath=expand(outdir+'figs/historical/delta_{vName}_historical_ts.png',
            vName=['emidust'])
    
    params:
        regrid=True

    notebook:
        "../notebooks/plot_change_time_series.py.ipynb"




rule plot_effect_of_landuse_change:
    input:
        pertubation = expand(outdir+'historical/{vName}/{vName}_historical_{model}_Ayear.nc',
                        model=['NorESM2-LM','CNRM-ESM2-1','MIROC-ES2L'],allow_missing=True),
        baseline = expand(outdir+'hist-noLu/{vName}/{vName}_hist-noLu_{model}_Ayear.nc',
                        model=['NorESM2-LM','CNRM-ESM2-1','MIROC-ES2L'],allow_missing=True),
    
        area_weights = expand('workflow/input_data/gridarea_{model}.nc',model=['NorESM2-LM','CNRM-ESM2-1'])
    output:
        outpath=outdir+'figs/hist_noLu/delta_{vName}_hist-noLu_ts.png'
    
    params:
        regrid=True

    notebook:
        "../notebooks/plot_change_time_series.py.ipynb"

