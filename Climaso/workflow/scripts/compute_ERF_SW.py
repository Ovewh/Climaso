from pyclim_noresm.aerosol_forcing import merge_exp_ctrl,calc_SW_ERF
from pyclim_noresm.general_util_funcs import yearly_avg
from .utils import load_CMIP_data

VARS = snakemake.config['variables']

vName_dw_SW = VARS[snakemake.wildcards.vName][1]
vName_up_SW = VARS[snakemake.wildcards.vName][0]
exp_dw_SW = load_CMIP_data(snakemake.input.exp_downwelling_SW, data_vars=[vName_dw_SW])
exp_up_SW = load_CMIP_data(snakemake.input.exp_upwelling_SW, data_vars=[vName_up_SW])
ctrl_dw_SW = load_CMIP_data(snakemake.input.ctrl_downwelling_SW,data_vars=[vName_dw_SW])
ctrl_up_SW = load_CMIP_data(snakemake.input.ctrl_upwelling_SW, data_vars=[vName_up_SW])
up_SW = merge_exp_ctrl(exp_up_SW, ctrl_up_SW)
dw_SW = merge_exp_ctrl(exp_dw_SW, ctrl_dw_SW)

ERF = calc_SW_ERF(dw_SW[vName_dw_SW], up_SW[vName_up_SW],
                dw_SW[f'control_{vName_dw_SW}'].rename(vName_dw_SW)
                , up_SW[f'control_{vName_up_SW}'].rename(vName_up_SW))
if snakemake.wildcards.freq == 'Ayear':
    ERF = yearly_avg(ERF)

ERF = ERF.dataset(name=snakemake.wildcards.vName)
ERF.to_netcdf(snakemake.output.outpath)
