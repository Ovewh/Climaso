from pyclim_noresm.aerosol_forcing import calc_atm_abs


def glob_control_path(w, var, subdir='Amon'):
    root_path = f'{ROOT_PATH}/{CMIP_VER}/'
    version = config['version'].get(w.model, 'latest')
    wildcards = glob_wildcards(root_path+'RFMIP/{institution}/'+f'{w.model}' +
                        '/piClim-control/{variant}/'+f'{subdir}/{var}/'+'{grid_label}/'+version+'/{fname}',
                        followlinks=True)
    fnames = []

    if wildcards.fname:
    
        if config['default_variant'] in wildcards.variant:
            variant = config['default_variant']
            for f in wildcards.fname:
                if variant in f.split('_'):
                    fnames.append(f)
        else:
            variant = wildcards.variant[0]
            for f in wildcards.fname:
                if variant in f.split('_'):
                    fnames.append(f)
        paths = expand('{root_path}RFMIP/{institution}/{model}/piClim-control/{variant}/{subdir}/{var}/{grid_label}/'+version+'/{fname}'
                    , root_path=[root_path], institution=wildcards.institution[0], model=w.model, variant=variant, subdir=subdir,var =var,
                    grid_label=wildcards.grid_label[0], fname=fnames)
        return paths
    else:
        wildcards = glob_wildcards(root_path+'AerChemMIP/{institution}/'+f'{w.model}' +
                        '/piClim-control/{variant}/'+f'{subdir}/{var}/'+'{grid_label}/'+version+'/{fname}',
                        followlinks=True)
        if config['default_variant'] in wildcards.variant:
            variant = config['default_variant']
            for f in wildcards.fname:
                if variant in f.split('_'):
                    fnames.append(f)
        else:
            variant = wildcards.variant[0]
            for f in wildcards.fname:
                if variant in f.split('_'):
                    fnames.append(f) 
        paths = expand('{root_path}AerChemMIP/{institution}/{model}/piClim-control/{variant}/{subdir}/{var}/{grid_label}/'+version+'/{fname}'
                    , root_path=[root_path], institution=wildcards.institution[0], model=w.model, variant=variant, subdir=subdir, var=var,
                    grid_label=wildcards.grid_label[0], fname=fnames)
        return paths

rule calc_ERF_surf:
    input:
        exp_downwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][1]}/**/latest/*.nc'),
        exp_upwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][0]}/**/latest/*.nc'),
        exp_upwelling_LW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][2]}/**/latest/*.nc'),
        exp_downwelling_LW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][3]}/**/latest/*.nc'),
        ctrl_downwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][1], subdir='Amon'),
        ctrl_upwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][0],subdir='Amon'),
        ctrl_upwelling_LW = lambda w:glob_control_path(w, var=VARS[w.vName][2], subdir='Amon'),
        ctrl_downwelling_LW = lambda w:glob_control_path(w, var=VARS[w.vName][3], subdir='Amon')

    output:
        outpath ='results/{vName}_{experiment}/{vName}_{experiment}_{model}_{freq}.nc'
    
    wildcard_constraints:
        vName = 'ERFsurf|ERFsurfcs'

    log:
        "logs/calc_ERF_toa/{vName}_{model}_{experiment}_{freq}.log"
    
    script:
        "../scripts/compute_ERF_surf.py"

    
    

rule calculate_ERF_TOA:
    input:
        exp_downwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][1]}/**/latest/*.nc'),
        exp_upwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][0]}/**/latest/*.nc'),
        exp_upwelling_LW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][2]}/**/latest/*.nc'),
        ctrl_downwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][1], subdir='Amon'),
        ctrl_upwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][0],subdir='Amon'),
        ctrl_upwelling_LW = lambda w:glob_control_path(w, var=VARS[w.vName][2], subdir='Amon')
    
    output:
        outpath = 'results/{vName}_{experiment}/{vName}_{experiment}_{model}_{freq}.nc'
    
    log:
        "logs/calc_ERF_toa/{vName}_{model}_{experiment}_{freq}.log"

    wildcard_constraints:
        vName = 'ERFt|ERFtcs'
    
    script:
        "../scripts/compute_ERF_TOA.py"



rule calculate_SW_ERF:
    input:
        exp_downwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][1]}/**/latest/*.nc'),
        exp_upwelling_SW = lambda w: glob.glob(f'{ROOT_PATH}/{CMIP_VER}/**/**/{w.model}/' +
                                f'{w.experiment}/**/Amon/{VARS[w.vName][0]}/**/latest/*.nc'),
        ctrl_downwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][1], subdir='Amon'),
        ctrl_upwelling_SW = lambda w:glob_control_path(w, var=VARS[w.vName][0],subdir='Amon')
    
    output:
        outpath = 'results/{vName}_{experiment}/{vName}_{experiment}_{model}_{freq}.nc'
    
    wildcard_constraints:
        vName='ERFtsw|ERFtswcs|ERFsurfsw|ERFsurfswcs'
    log:
        "logs/calc_ERF_SW/{vName}_{model}_{experiment}_{freq}.log"  
    script:
        "../scripts/compute_ERF_SW.py"


rule calc_absorption:
    input:
        delta_rad_surf = lambda w: f'results/{VARS[w.vName][1]}_{w.experiment}/{VARS[w.vName][1]}_{w.experiment}_{w.model}_{w.freq}.nc',
        delta_rad_toa = lambda w: f'results/{VARS[w.vName][0]}_{w.experiment}/{VARS[w.vName][0]}_{w.experiment}_{w.model}_{w.freq}.nc'
    output:
        outpath='results/atm_abs/{vName}_{experiment}_{model}_{freq}.nc'
    wildcard_constraints:
        vName='atmabsSW|atmabs'
    log:
        "logs/calc_absorption/{vName}_{model}_{experiment}_{freq}.log"
    run:
        import xarray as xr
        vName_rad_surf = VARS[wildcards.vName][1]
        vName_rad_toa = VARS[wildcards.vName][0]
        delta_rad_toa = xr.open_dataset(input.delta_rad_toa)
        delta_rad_surf = xr.open_dataset(input.delta_rad_surf)
        atm_abs = calc_atm_abs(delta_rad_surf[vName_rad_surf],delta_rad_toa[vName_rad_toa])

        atm_abs = atm_abs.to_dataset(name=wildcards.vName)
        atm_abs.to_netcdf(output.outpath)