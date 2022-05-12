
rule calc_feedback:
    input:
        emis_exp = expand(rules.calc_experiment_climalogies.output.outpath, freq='Ayear', allow_missing=True),
        emis_ctrl = expand(rules.calc_clim_PI_control.output.outpath, experiment='piClim-control',
                    freq='Ayear',allow_missing=True),
        emis_PIctrl = expand(rules.calc_experiment_climalogies.output.outpath, freq='Ayear', 
                        experiment='piControl',allow_missing=True),
        emis_2xCO2 = expand(rules.calc_experiment_climalogies.output.outpath, freq='Ayear', 
                        experiment='abrupt-4xCO2',allow_missing=True),
        t_ctrl = expand(rules.calc_experiment_climalogies.output.outpath, freq='Ayear', variable='ts',
                        experiment='piControl',allow_missing=True),
        t_exp = expand(rules.calc_experiment_climalogies.output.outpath, freq='Ayear', variable='ts',
                        experiment='abrupt-4xCO2',allow_missing=True),
        area_cello='workflow/input_data/gridarea_{model}.nc',
        forcing =outdir + '{experiment}/ERFs/{erf}/{erf}_{experiment}_{model}_Ayear.nc',


    output:
        outpath=outdir+'{experiment}/Feedback_per_emis/{erf}_{variable}_{experiment}_{model}_Ayear.yaml'
    wildcard_constraints:
        erf='|'.join(config['variables'].keys())

    notebook:
        "../notebooks/calc_feedback_per_emissions.py.ipynb"

#rule calc_forcing_per_burden_change:
