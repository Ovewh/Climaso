

rule make_local_catalogue:
    output:
        outpath = 'catalogues/{activity}_{source}_CMIP6.csv.gz'
    params:
        root_path = lambda w: config[f"root_{w.source}"] + f'/{w.activity}/',
        depth = 5
    shell:
        "python workflow/scripts/builders/cmip.py --root-path {params.root_path} --csv-filepath {output.outpath} -d {params.depth} -v 6 --pick-latest-version"

rule build_catalogues:
    input:
        expand('catalogues/{activity}_{source}_CMIP6.csv.gz', activity = config['activities'], source = ['betzy', 'noresm'])
    output:
        table='catalogues/merge_CMIP6.csv',
        json='catalogues/merge_CMIP6.json'
    notebook:
        '../notebooks/merge_catalogues.py.ipynb'        

rule get_data_intake:
    input:
        catalog = rules.build_catalogues.output.json
    output:
        outpath = output_format['single_variable']
    params:
        accumalative_vars = config['accumalative_vars']
    
    log:
        "logs/calc_clim/{variable}_{model}_{experiment}_{freq}.log"
    notebook:
        "../notebooks/get_data_intake.py.ipynb"


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

