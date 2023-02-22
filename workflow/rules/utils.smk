

rule make_local_catalogue:
    output:
        outpath = 'catalogues/{activity}_{source}_CMIP6.csv.gz'
    params:
        root_path = lambda w: config[f"root_{w.source}"] + f'/{w.activity}/',
        depth = 9,
    threads: 4
    run:
        import ecgtools
        from ecgtools import Builder
        from ecgtools.parsers.cmip import parse_cmip6_using_directories

        if params.root_path.split('/')[-2] == 'NIRD_noresm':
            exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*']
        else:
            exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*','*/NorESM2-LM/*','*/NorESM2-MM/*']        
        builder = Builder(paths=[params.root_path], depth=params.depth,
                        joblib_parallel_kwargs={'n_jobs': threads, 'verbose':13},
                        exclude_patterns=['*/files/*', '*/latest/*','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*','*/NorESM2-LM/*','*/NorESM2-MM/*'])
        builder.build(parsing_func=parse_cmip6_using_directories)

        builder.clean_dataframe()
        df = builder.df
        df.to_csv(output.outpath, compression='gzip', index=False)
    # shell:
        # "python workflow/scripts/builders/cmip.py --root-path {params.root_path} --csv-filepath {output.outpath} -d {params.depth} -v 6 --pick-latest-version"

rule build_catalogues:
    input:
        expand('catalogues/{activity}_{source}_CMIP6.csv.gz', 
                activity = config['activities'], 
                source = ['betzy', 'noresm']),
        'catalogues/AerChemMIP_noresmdev_CMIP6.csv.gz',
        'catalogues/CMIP_nirdCMIPtemp_CMIP6.csv.gz'
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

def trans_mmr_to_load(mmr:str):
    trans_dict = {'concdust':'mmrdust',
                'concpm1':'mmrpm1',
                'concpm10':'mmrpm10',
                'concpm2p5':'mmrpm2p5',
                'concso4':'mmrso4',
                'concss':'mmrss',
                'concsoa':'mmrsoa',
                'concoa':'mmroa',
                'conch2oaer':'mmraerh2o'}
    return trans_dict[mmr]

rule derive_column_integrated_load_airmass:
    input:
        mmr = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable=trans_mmr_to_load(w.variable)),
        airmass = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable='airmass'),
    output:
        outpath = outdir + '{experiment}/derived_variables/{variable}/{variable}_{model}_{experiment}_{freq}.nc'
    wildcard_constraints:
        model='UKESM1-0-LL',
        variables = 'concdust|concpm1|concpm10|concpm2p5|concso4|concss|concsoa|concoa|conch2oaer'
    notebook:
        "../notebooks/derive_column_integrated_load.py.ipynb"

rule derive_column_integrated_load:
    input:
        mmr = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable=trans_mmr_to_load(w.variable)),
    output:
        outpath = outdir + '{experiment}/derived_variables/{variable}/{variable}_{model}_{experiment}_{freq}.nc'
    wildcard_constraints:
        variables = 'concdust|concpm1|concpm10|concpm2p5|concso4|concss|concsoa|concoa|conch2oaer',
        model="(?!UKESM1-0-LL).*"
    notebook:
        "../notebooks/derive_column_integrated_load.py.ipynb"

rule column_integrate_cdnc:
    input:
        cdnc = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable='cdnc'),
        ta = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable='ta'),
    output:
        outpath = outdir + '{experiment}/derived_variables/cdncvi/cdncvi_{model}_{experiment}_{freq}.nc'
    conda:
        "../envs/comp_cat.yaml"
    
    notebook:
        "../notebooks/derive_column_integrated_cdnc.py.ipynb"

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

