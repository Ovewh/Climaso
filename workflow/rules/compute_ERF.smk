

rule calc_ERF_surf:
    input:
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],w.experiment),
        exp_upwelling_SW = lambda w: get_paths(w,VARS[w.vName][0],w.experiment),
        exp_upwelling_LW = lambda w: get_paths(w,VARS[w.vName][2],w.experiment),
        exp_downwelling_LW = lambda w: get_paths(w,VARS[w.vName][3],w.experiment),
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1]),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0]),
        ctrl_upwelling_LW = lambda w:get_control_path(w, VARS[w.vName][2]),
        ctrl_downwelling_LW = lambda w:get_control_path(w, VARS[w.vName][3])

    output:
        outpath = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    
    wildcard_constraints:
        vName = 'ERFsurf|ERFsurfcs'

    log:
        "logs/calc_ERF_surf/{vName}_{model}_{experiment}_{freq}.log"
    
    script:
        "../scripts/compute_ERF_surf.py"


rule calculate_ERF_TOA:
    input:
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],w.experiment),
        exp_upwelling_SW =lambda w: get_paths(w,VARS[w.vName][0],w.experiment),
        exp_upwelling_LW = lambda w: get_paths(w,VARS[w.vName][2],w.experiment),
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1]),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0]),
        ctrl_upwelling_LW = lambda w:get_control_path(w, VARS[w.vName][2])
    
    output:
        outpath = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    
    log:
        "logs/calc_ERF_toa/{vName}_{model}_{experiment}_{freq}.log"

    wildcard_constraints:
        vName = 'ERFt|ERFtcs|ERFtaf|ERFtcsaf'
    
    script:
        "../scripts/compute_ERF_TOA.py"



rule calculate_SW_ERF:
    input:
        exp_downwelling_SW = lambda w: get_paths(w,VARS[w.vName][1],w.experiment),
        exp_upwelling_SW = lambda w: get_paths(w,VARS[w.vName][0],w.experiment), 
        ctrl_downwelling_SW = lambda w:get_control_path(w, VARS[w.vName][1]),
        ctrl_upwelling_SW = lambda w:get_control_path(w, VARS[w.vName][0])
    
    output:
        outpath = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    
    wildcard_constraints:
        vName='ERFtsw|ERFtswcs|ERFsurfsw|ERFsurfswcs|ERFtswcsaf|ERFtswaf'
    log:
        "logs/calc_ERF_SW/{vName}_{model}_{experiment}_{freq}.log"  
    script:
        "../scripts/compute_ERF_SW.py"

rule calculate_ERF_TOA_LW:
    input:
        exp_upwelling_LW = lambda w: get_paths(w,VARS[w.vName][0],w.experiment),
        ctrl_upwelling_LW = lambda w: get_control_path(w,VARS[w.vName][0])
    output:
        outpath = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    log:
        "logs/calc_ERF_LW/{vName}_{model}_{experiment}_{freq}.log"  
    wildcard_constraints:
        vName='ERFtlw|ERFtlwaf|ERFtlwcs|ERFtlwcsaf'

    script:
        "../scripts/compute_ERF_LW_TOA.py"


rule calc_cloud_radiative_effect:
    input:
        ERFaf =  lambda w: outdir + f'{w.experiment}/ERFs/{VARS[w.vName][0]}/{VARS[w.vName][0]}_{w.experiment}_{w.model}_{w.freq}.nc',
        ERFafcs = lambda w: outdir + f'{w.experiment}/ERFs/{VARS[w.vName][1]}/{VARS[w.vName][1]}_{w.experiment}_{w.model}_{w.freq}.nc'
    output:
        outpath= outdir+'{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    wildcard_constraints:
        vName='CloudEff|SWCloudEff|LWCloudEff'

    log:
        "logs/{vName}_{model}_{experiment}_{freq}.log"
    script:
        "../scripts/compute_cloud_effect.py"

rule calc_direct_radiative_effect:
    input:
        ERFt =  lambda w:  outdir+ f'{w.experiment}/ERFs/{VARS[w.vName][0]}/{VARS[w.vName][0]}_{w.experiment}_{w.model}_{w.freq}.nc',
        ERFtaf = lambda w: outdir +f'{w.experiment}/ERFs/{VARS[w.vName][1]}/{VARS[w.vName][1]}_{w.experiment}_{w.model}_{w.freq}.nc'
    output:
        outpath=outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'

    log:
        "logs/rad/{vName}_{model}_{experiment}_{freq}.log"
    wildcard_constraints:
        vName='SWDirectEff|SWDirectEff_cs|LWDirectEff|LWDirectEff_cs|DirectEff'
    script:
        "../scripts/compute_direct_radiative_effect.py"

rule calc_absorption:
    input:
        delta_rad_surf = lambda w: outdir + f'{w.experiment}/ERFs/{VARS[w.vName][1]}/{VARS[w.vName][1]}_{w.experiment}_{w.model}_{w.freq}.nc',
        delta_rad_toa = lambda w: outdir + f'{w.experiment}/ERFs/{VARS[w.vName][0]}/{VARS[w.vName][0]}_{w.experiment}_{w.model}_{w.freq}.nc'
    output:
        outpath= outdir+'{experiment}/ERFs/{vName}/{vName}_{experiment}_{model}_{freq}.nc'
    wildcard_constraints:
        vName='atmabsSW|atmabs'
    log:
        "logs/calc_absorption/{vName}_{model}_{experiment}_{freq}.log"
    run:
        import xarray as xr
        from pyclim_noresm.aerosol_forcing import calc_atm_abs
        import time
        vName_rad_surf = VARS[wildcards.vName][1]
        vName_rad_toa = VARS[wildcards.vName][0]
        delta_rad_toa = xr.open_dataset(input.delta_rad_toa)
        delta_rad_surf = xr.open_dataset(input.delta_rad_surf)
        atm_abs = calc_atm_abs(delta_rad_surf[vName_rad_surf],delta_rad_toa[vName_rad_toa])

        atm_abs = atm_abs.to_dataset(name=wildcards.vName)
        atm_abs.attrs=copy.copy(delta_rad_surf.attrs)
        atm_abs.attrs['title'] = 'Atmospheric absorbtion'
        atm_abs.attrs['history'] = f'@{time.ctime()} Generated by: {__file__} ' + atm_abs.attrs['history']
        atm_abs.attrs['source'] = ', '.join(input) + ', ' + atm_abs.attrs['source']
        atm_abs.to_netcdf(output.outpath)