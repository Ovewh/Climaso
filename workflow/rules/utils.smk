
# rule calc_dust_loading:
#     input:
#         mmr = lambda w: get_paths(w,'mmrdust',w.experiment,grid_label=config['default_grid_label']),
#         air_mass = lambda w: get_paths(w,'airmass',w.experiment,grid_label=config['default_grid_label'])
#     output:
#         outpath = expand(output_format['single_variable'],variable=['loaddust'],allow_missing=True)
#     wildcard_constraints:
#         model = 'EC-Earth3-AerChem|MPI-ESM-1-2-HAM',
#         experiment="|".join(EXPERIMENTS)

#     priority: 50
#     notebook:
#         "../notebooks/calc_load.py.ipynb"

# rule calc_dust_loading_ctrl:
#     input:
#         mmr = lambda w: get_control_path(w, 'mmrdust', 'gn'),
#         # air_mass = lambda w: get_control_path(w,'airmass', 'gn')
#     output:
#         outpath = expand(output_format['single_variable'],variable=['loaddust'],allow_missing=True)
#     wildcard_constraints:
#         model = 'EC-Earth3-AerChem|MPI-ESM-1-2-HAM',
#         # experiment="|".join(CONTROL_EXPS)
    
#     priority: 50
#     notebook:
#         "../notebooks/calc_load.py.ipynb"

rule calc_clim_PI_control:
    input:
        pi_clim_var = lambda w: get_control_path(w, w.variable),
    output:
        outpath = output_format['single_variable']
    
    log:
        "logs/calc_clim/{variable}_{model}_{experiment}_{freq}.log"
    wildcard_constraints:
        experiment="|".join(CONTROL_EXPS),
    params:
        accumalative_vars = config['accumalative_vars']
    notebook:
        "../notebooks/calc_clim.py.ipynb"

rule calc_experiment_climalogies:
    input: 
        pi_clim_var = lambda w: get_paths(w,w.variable,w.experiment,grid_label=config['default_grid_label'])
    output:
        outpath = output_format['single_variable']
    log:
        "logs/calc_clim/{variable}_{model}_{experiment}_{freq}.log"
    wildcard_constraints:
        experiment="|".join(EXPERIMENTS),

    params:
        accumalative_vars = config['accumalative_vars']
    notebook:
        "../notebooks/calc_clim.py.ipynb"

rule change_historical_ts:
    input:
        pertubation = outdir + '{hist_pert}/{vName}/{vName}_{hist_pert}_{model}_Ayear.nc',
        baseline = outdir + '{hist_base}/{vName}/{vName}_{hist_base}_{model}_Ayear.nc',
    output:
        outpath=outdir+'figs/{hist_pert}-{hist_base}_ts/{vName}/{vName}_{model}_delta_{hist_pert}_{hist_base}.{ext}'
    wildcard_constraints:
        ext = 'csv|png|pdf'
    notebook:
        "../notebooks/plot_change_time_series.py.ipynb"

rule calc_change_historical:
    input:    
        pertubation = rules.change_historical_ts.input.pertubation,
        baseline = rules.change_historical_ts.input.baseline
    
    output:
        outpath = outdir + '{hist_pert}-{hist_base}/{vName}/{vName}_{model}_delta_{hist_pert}_{hist_base}.nc'

    notebook:
        "../notebooks/calc_difference.py.ipynb"

rule cmip6_to_aerocom_fmt:
    input:
        paths = lambda w: get_paths(w,w.variable,w.experiment,grid_label=config['default_grid_label'])
    output:
        outpath = outdir+'{model}_{experiment}/renamed/converted_CMIP6_aerocom_{variable}.txt'
    
    notebook:
        "../convert_to_aerocom_fmt.py.ipynb"

