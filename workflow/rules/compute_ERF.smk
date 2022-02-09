
def get_control_path(w, variable,table_id, grid_label):
    try:
        paths = get_paths(w, variable,'piClim-control',table_id, grid_label,activity='RFMIP')
    except KeyError:
        paths = get_paths(w, variable,'piClim-control',table_id, grid_label,activity='AerChemMIP')
    return paths


rule calc_ERF_surf:
    input:
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],w.experiment,'Amon','gn'),
        exp_upwelling_SW = lambda w: get_paths(w,VARS[w.vName][0],w.experiment,'Amon','gn'),
        exp_upwelling_LW = lambda w: get_paths(w,VARS[w.vName][2],w.experiment,'Amon','gn'),
        exp_downwelling_LW = lambda w: get_paths(w,VARS[w.vName][3],w.experiment,'Amon','gn'),
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1], 'Amon','gn'),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0],'Amon','gn'),
        ctrl_upwelling_LW = lambda w:get_control_path(w, VARS[w.vName][2], 'Amon','gn'),
        ctrl_downwelling_LW = lambda w:get_control_path(w, VARS[w.vName][3],'Amon','gn')

    output:
        outpath ='results/{vName}_{experiment}/{vName}_{experiment}_{model}_{freq}.nc'
    
    wildcard_constraints:
        vName = 'ERFsurf|ERFsurfcs'

    log:
        "logs/calc_ERF_surf/{vName}_{model}_{experiment}_{freq}.log"
    
    script:
        "../scripts/compute_ERF_surf.py"

    
    

rule calculate_ERF_TOA:
    input:
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],w.experiment,'Amon','gn'),
        exp_upwelling_SW =lambda w: get_paths(w,VARS[w.vName][0],w.experiment,'Amon','gn'),
        exp_upwelling_LW = lambda w: get_paths(w,VARS[w.vName][2],w.experiment,'Amon','gn'),
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1], 'Amon', 'gn'),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0],'Amon', 'gn'),
        ctrl_upwelling_LW = lambda w:get_control_path(w, VARS[w.vName][2], 'Amon', 'gn')
    
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
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],'Amon','gn'),
        exp_upwelling_SW = lambda w: get_paths(w,VARS[w.vName][0],w.experiment,'Amon','gn'), 
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1], 'Amon', 'gn'),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0],'Amon', 'gn')
    
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
        from pyclim_noresm.aerosol_forcing import calc_atm_abs
        vName_rad_surf = VARS[wildcards.vName][1]
        vName_rad_toa = VARS[wildcards.vName][0]
        delta_rad_toa = xr.open_dataset(input.delta_rad_toa)
        delta_rad_surf = xr.open_dataset(input.delta_rad_surf)
        atm_abs = calc_atm_abs(delta_rad_surf[vName_rad_surf],delta_rad_toa[vName_rad_toa])

        atm_abs = atm_abs.to_dataset(name=wildcards.vName)
        atm_abs.to_netcdf(output.outpath)