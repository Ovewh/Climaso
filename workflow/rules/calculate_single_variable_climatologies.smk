

rule calc_clim_PI_control:
    input:
        pi_clim_var = lambda w: get_control_path(w, w.variable),
    output:
        outpath = output_format['single_variable']
    
    log:
        "logs/calc_clim/{variable}_{model}_{experiment}_{freq}.log"
    wildcard_constraints:
        experiment="|".join(CONTROL_EXPS),
        variable="!".join(config['variables'].keys())

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
        variable="!".join(config['variables'].keys())


    notebook:
        "../notebooks/calc_clim.py.ipynb"
