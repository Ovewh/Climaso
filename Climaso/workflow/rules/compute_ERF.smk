
import os
import glob
VARS = config['variables']


def glob_control_path(w, var='', subdir='Amon'):
    paths = glob.glob(f'{ROOT_PATH}/{CMIP_VER}/RFMIP/**/{w.model}/' +
                        f'piClim-control/**/{subdir}/{var}/**/latest/*.nc')
    if paths:
        return paths
    else:
        paths = glob.glob(f'{ROOT_PATH}/{CMIP_VER}/AerChemMIP/**/{w.model}/' +
                        f'piClim-control/**/{subdir}/{var}/**/latest/*.nc')
        return paths


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
        outpath = 'results/{vName}_{experiment}_{model}_{freq}.nc'
    
    run:
        from pyclim_noresm.aerosol_forcing import calc_total_ERF_TOA, merge_exp_ctrl
        from pyclim_noresm.general_util_funcs import yearly_avg
        import xarray as xr
        vName_dw_SW = VARS[wildcards.vName][1]
        vName_up_SW = VARS[wildcards.vName][0]
        vName_up_LW = VARS[wildcards.vName][2]
        if len(input.exp_downwelling_SW)>1:
            exp_dw_SW = xr.open_mfdataset(input.exp_downwelling_SW, chunks={'time':120}, 
                                        data_vars=[vName_dw_SW])
            exp_up_SW = xr.open_mfdataset(input.exp_upwelling_SW, chunks={'time':120},
                                        data_vars=[vName_up_SW])
            exp_up_LW = xr.open_mfdataset(input.exp_upwelling_LW, chunks={'time':120},
                                        data_vars=[vName_up_LW])
            ctrl_dw_SW = xr.open_mfdataset(input.ctrl_downwelling_SW, chunks={'time':120},
                                        data_vars=[vName_dw_SW])
            ctrl_up_SW = xr.open_mfdataset(input.ctrl_upwelling_SW,chunks={'time':120},
                                        data_vars=[vName_up_SW])
            ctrl_up_LW = xr.open_mfdataset(input.ctrl_upwelling_LW,chunks={'time':120},
                                        data_vars=[vName_up_LW])
        else:
            exp_dw_SW = xr.open_dataset(input.exp_downwelling_SW, chunks={'time':120})
            exp_up_SW = xr.open_dataset(input.exp_upwelling_SW, chunks={'time':120})
            exp_up_LW = xr.open_dataset(input.exp_upwelling_LW, chunks={'time':120})
            ctrl_dw_SW = xr.open_dataset(input.ctrl_downwelling_SW, chunks={'time':120})
            ctrl_up_SW = xr.open_dataset(input.ctrl_upwelling_SW,chunks={'time':120})
            ctrl_up_LW = xr.open_dataset(input.ctrl_upwelling_LW,chunks={'time':120})

        dw_SW = merge_exp_ctrl(exp_dw_SW, ctrl_dw_SW)
        up_LW = merge_exp_ctrl(exp_up_LW, ctrl_up_LW)
        up_SW = merge_exp_ctrl(exp_up_SW, ctrl_up_SW)

        ERF = calc_total_ERF_TOA(dw_SW[vName_dw_SW], up_SW[vName_up_SW],
                        up_LW[vName_up_LW],
                        dw_SW[f'control_{vName_dw_SW}'].rename(vName_dw_SW)
                        , up_SW[f'control_{vName_up_SW}'].rename(vName_up_SW),
                        up_LW[f'control_{vName_up_LW}'].rename(vName_up_LW))

        if wildcards.freq == 'Ayear':
            ERF = yearly_avg(ERF)

        ERF = ERF.dataset(name=wildcards.vName)
        ERF.to_netcdf(output.outpath)



        



