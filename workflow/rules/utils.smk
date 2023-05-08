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
        lustre = config.get('seperate_Noresm', True)
        if wildcards.source == 'noresm':
            exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*']
        elif lustre == False:
            exclude_patterns=['*/files/*']
        
        else:
            exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*','*/NorESM2-LM/*','*/NorESM2-MM/*']    
            
        print(exclude_patterns)
        builder = Builder(paths=[params.root_path], depth=params.depth,
                        joblib_parallel_kwargs={'n_jobs': threads, 'verbose':13},
                        exclude_patterns=exclude_patterns)
        builder.build(parsing_func=parse_cmip6_using_directories)

        builder.clean_dataframe()
        df = builder.df
        df.to_csv(output.outpath, compression='gzip', index=False)

rule make_available_data_tracker:
    input:
        reference = 'config/reference_data_request.yaml',
    
    output:
        outpath = 'config/.data_trackers/{experiment}_{source}_CMIP6.yaml'

    notebook:
        '../notebooks/make_available_data_tracker.py.ipynb'


rule build_catalogues:
    input:
        expand('catalogues/{activity}_{source}_CMIP6.csv.gz', 
                activity = config['activities'], 
                source = config['sources']),
        'catalogues/AerChemMIP_noresm_CMIP6.csv.gz',
        'catalogues/RFMIP_noresm_CMIP6.csv.gz',
        # 'catalogues/CMIP_nirdCMIPtemp_CMIP6.csv.gz'
    output:
        table='catalogues/merge_CMIP6.csv',
        json='catalogues/merge_CMIP6.json'
    notebook:
        '../notebooks/merge_catalogues.py.ipynb'        

rule get_data_intake:
    input:
        catalog = ancient(rules.build_catalogues.output.json),
        data_tracker = ancient('config/.data_trackers/{experiment}_{model}_CMIP6.yaml')
    output:
        outpath = output_format['single_variable']
    params:
        accumalative_vars = config['accumalative_vars']
    
    log:
        "logs/calc_clim/{variable}_{model}_{experiment}_{freq}.log"
    notebook:
        "../notebooks/get_data_intake.py.ipynb"





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

    wildcard_constraints:
        model="(?!UKESM1-0-LL).*"
    
    notebook:
        "../notebooks/derive_column_integrated_cdnc.py.ipynb"


rule column_integrate_cdnc_UKESM:
    input:
        cdnc = lambda w: expand(output_format['single_variable'], model='UKESM1-0-LL', experiment=w.experiment,
                freq=w.freq, variable='cdnc'),
        ta = lambda w: expand(output_format['single_variable'], model='UKESM1-0-LL', experiment=w.experiment,
                freq=w.freq, variable='ta'),
        pfull = lambda w: expand(output_format['single_variable'], model='UKESM1-0-LL', experiment=w.experiment,
                freq=w.freq, variable='pfull'),
    output:
        outpath = outdir + '{experiment}/derived_variables/cdncvi/cdncvi_UKESM1-0-LL_{experiment}_{freq}.nc'
    conda:
        "../envs/comp_cat.yaml"
    
    notebook:
        "../notebooks/derive_column_integrated_cdnc.py.ipynb"




rule derive_column_integrated_load_airmass:
    input:
        mmr = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable=config['burdens_dict'].get(w.variable)),
        airmass = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable='airmass'),
    output:
        outpath = outdir + '{experiment}/derived_variables/{variable}/{variable}_{model}_{experiment}_{freq}.nc'
    wildcard_constraints:
        model='UKESM1-0-LL',
        variable = 'concdust|concpm1|concpm10|concpm2p5|concso4|concss|concsoa|concoa|conch2oaer'
    notebook:
        "../notebooks/derive_column_integrated_load.py.ipynb"

rule derive_column_integrated_load:
    input:
        mmr = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                freq=w.freq, variable=config['burdens_dict'].get(w.variable)),
    output:
        outpath = outdir + '{experiment}/derived_variables/{variable}/{variable}_{model}_{experiment}_{freq}.nc'
    wildcard_constraints:
        variable = 'concdust|concpm1|concpm10|concpm2p5|concso4|concss|concsoa|concoa|conch2oaer',
        model="(?!UKESM1-0-LL).*"
    notebook:
        "../notebooks/derive_column_integrated_load.py.ipynb"



rule derived_windspeed:
    input:
        u = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                    freq=w.freq, variable='ua'),
        v = lambda w: expand(output_format['single_variable'], model=w.model, experiment=w.experiment,
                    freq=w.freq, variable='va'),
    output:
        outpath = outdir + '{experiment}/derived_variables/windspeed/windspeed_{model}_{experiment}_{freq}.nc'
    notebook:
        "../notebooks/derive_windspeed.py.ipynb"

rule cmip6_to_aerocom_fmt:
    input:
        paths = lambda w: get_paths(w,w.variable,w.experiment,grid_label=config['default_grid_label'])
    output:
        outpath = outdir+'{model}_{experiment}/renamed/converted_CMIP6_aerocom_{variable}.txt'
    
    notebook:
        "../convert_to_aerocom_fmt.py.ipynb"

