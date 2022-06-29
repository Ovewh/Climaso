
rule diagnose_cloud_changes_piClim:
    input:
        path_exp =  expand(rules.calc_experiment_climalogies.output.outpath,freq='Ayear', allow_missing=True),
        path_ctrl = expand(rules.calc_clim_PI_control.output.outpath,freq='Ayear', experiment='piClim-control', 
                        allow_missing=True)
        # ps_ctrl = expand(rules.calc_clim_PI_control.output.outpath,freq='Ayear', experiment='piClim-control', 
                        # variable='ps', allow_missing=True), 
        # ps_exp = expand(rules.calc_experiment_climalogies.output.outpath,freq='Ayear', variable='ps' ,allow_missing=True)
    output:
        outpath = outdir + '{experiment}/delta_{variable}/delta_{variable}_{experiment}_{model}.nc'
    conda:
        "../envs/comp_cat.yaml"
    notebook:
        "../notebooks/diagnose_cloud_changes_piclim.py.ipynb"
