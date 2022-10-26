
rule diagnose_cloud_changes_piClim:
    input:
        path_exp =  expand(rules.get_data_intake.output.outpath,freq='Ayear', allow_missing=True),
        path_ctrl = expand(rules.get_data_intake.output.outpath,freq='Ayear', experiment='piClim-control', 
                        allow_missing=True)
        # ps_ctrl = expand(rules.get_data_intake.output.outpath,freq='Ayear', experiment='piClim-control', 
                        # variable='ps', allow_missing=True), 
        # ps_exp = expand(rules.get_data_intake.output.outpath,freq='Ayear', variable='ps' ,allow_missing=True)
    output:
        outpath = outdir + '{experiment}/delta_{variable}/delta_{variable}_{experiment}_{model}.nc'
    conda:
        "../envs/comp_cat.yaml"
    notebook:
        "../notebooks/diagnose_cloud_changes_piclim.py.ipynb"
