from pyclim_noresm.aerosol_forcing import calc_total_ERF_TOA, merge_exp_ctrl
from utils import load_CMIP_data, copy_meta_data_CMIP
from pyclim_noresm.general_util_funcs import yearly_avg
import time

VARS = snakemake.config['variables']

vName_dw_SW = VARS[snakemake.wildcards.vName][1]
vName_up_SW = VARS[snakemake.wildcards.vName][0]
vName_up_LW = VARS[snakemake.wildcards.vName][2]
exp_dw_SW = load_CMIP_data(snakemake.input.exp_downwelling_SW, data_vars=[vName_dw_SW])
exp_up_SW = load_CMIP_data(snakemake.input.exp_upwelling_SW, data_vars=[vName_up_SW])
exp_up_LW = load_CMIP_data(snakemake.input.exp_upwelling_LW, data_vars=[vName_up_LW])
ctrl_dw_SW = load_CMIP_data(snakemake.input.ctrl_downwelling_SW,data_vars=[vName_dw_SW])
ctrl_up_SW = load_CMIP_data(snakemake.input.ctrl_upwelling_SW, data_vars=[vName_up_SW])
ctrl_up_LW = load_CMIP_data(snakemake.input.ctrl_upwelling_LW,data_vars=[vName_up_LW])

dw_SW = merge_exp_ctrl(exp_dw_SW, ctrl_dw_SW)
up_LW = merge_exp_ctrl(exp_up_LW, ctrl_up_LW)
up_SW = merge_exp_ctrl(exp_up_SW, ctrl_up_SW)

ERF = calc_total_ERF_TOA(dw_SW[vName_dw_SW], up_SW[vName_up_SW],
                up_LW[vName_up_LW],
                dw_SW[f'control_{vName_dw_SW}'].rename(vName_dw_SW)
                , up_SW[f'control_{vName_up_SW}'].rename(vName_up_SW),
                up_LW[f'control_{vName_up_LW}'].rename(vName_up_LW))

if snakemake.wildcards.freq == 'Ayear':
    ERF = yearly_avg(ERF)

ERF = ERF.to_dataset(name=snakemake.wildcards.vName)
attrs = copy_meta_data_CMIP(exp_dw_SW.attrs)
ERF = ERF.assign_attrs(**attrs)
ERF = ERF.assign_attrs(varialbe_id=snakemake.wildcards.vName)
ERF.attrs['source'] = ', '.join(snakemake.input)
ERF.attrs['history'] = f'@{time.ctime()} Generated by: {__file__}'
ERF.attrs['frequency'] = snakemake.wildcards.freq
ERF.attrs['title'] = 'Effective radiative forcing TOA'
ERF.to_netcdf(snakemake.output.outpath)