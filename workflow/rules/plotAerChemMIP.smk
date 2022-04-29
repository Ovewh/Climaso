
rule plot_ERFs:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc', 
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    
    wildcard_constraints:
        vName = 'ERFtsw|ERFtlw|ERFtcs|ERFt'

    notebook:
        "../notebooks/plot_ERFs_AerChemMIP.py.ipynb"

rule plot_atm_abs:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc', 
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    
    wildcard_constraints:
        vName ='atmabsSW|atmabs'

    notebook:
        "../notebooks/plot_atmabs_AerChemMIP.py.ipynb"



rule plot_global_avaraged_ERFs:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc', 
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_globalaverage_AerChemMIP.png'
    wildcard_constraints:
        vName = 'ERFtsw|ERFtlw|ERFtcs|ERFt'

    notebook:
        "../notebooks/plot_ERFs_AerChemMIP_gobal_average.py.ipynb"

rule plot_ERFaci:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc',
            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL'],
            allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    wildcard_constraints:
        vName = 'SWDirectEff|LWDirectEff|SWCloudEff|LWCloudEff|LWDirectEff_cs|SWDirectEff_cs|CloudEff|DirectEff'
    
    notebook:
        "../notebooks/plot_Radiative_effects_AerChemMIP.py.ipynb"
    
rule plot_change_cdnc:
    input:
        path_exp=expand(outdir+'piClim-2xdust/cdnc/cdnc_piClim-2xdust_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem',  'GISS-E2-1-G', 
                        'MIROC6','UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/cdnc/cdnc_piClim-control_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem',  'GISS-E2-1-G', 
                        'MIROC6','UKESM1-0-LL',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/cdnc_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'

    notebook:
        "../notebooks/plot_CDNC_change.py.ipynb"

rule plot_change_lwp:
    input:
        path_exp = expand(outdir+'piClim-2xdust/lwp/lwp_piClim-2xdust_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem',
                        'UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/lwp/lwp_piClim-control_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem', 
                        'UKESM1-0-LL',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/lwp_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'

    notebook:
        "../notebooks/plot_lwp_change.py.ipynb"

rule plot_change_clt:
    input:
        path_exp = expand(outdir+'piClim-2xdust/clt/clt_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/clt/clt_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/clt_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'

    notebook:
        "../notebooks/plot_clt_change.py.ipynb"

rule plot_change_ts:
    input:
        path_exp = expand(outdir+'piClim-2xdust/ts/ts_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 
                         'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/ts/ts_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G',  
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/ts_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'

    notebook:
        "../notebooks/plot_ts_change.py.ipynb"

rule plot_change_tas:
    input:
        path_exp = expand(outdir+'piClim-2xdust/tas/tas_piClim-2xdust_{model}_Ayear.nc',
                 model=['GISS-E2-1-G', 
                        'MIROC6', 'GFDL-ESM4', 
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/tas/tas_piClim-control_{model}_Ayear.nc',
                 model=['GISS-E2-1-G',  
                        'MIROC6', 'GFDL-ESM4', 
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/tas_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'

    notebook:
        "../notebooks/plot_ts_change.py.ipynb"


rule plot_change_prs:
    input:
        path_exp = expand(outdir+'piClim-2xdust/pr/pr_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/pr/pr_piClim-control_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/pr_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs|rel'

    notebook:
        "../notebooks/plot_pr_change.py.ipynb"




rule plot_emidust:
    input:
        paths=expand(outdir+'{experiment}/emidust/emidust_{experiment}_{model}_Amon.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        areacello = expand('input_data/gridarea_{model}.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
            
    output:
        outpath=outdir+'figs/AerChemMIP/emidust_{experiment}_Ayear_map.png' 
    wildcard_constraints:
        experiment="|".join(ALL_ALLOWED)
    
    notebook:
        "../notebooks/plot_emidust.py.ipynb"


rule plot_depdust:
    input:
        paths=expand(outdir+'{experiment}/emidust/depdust_{experiment}_{model}_Amon.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        areacello = expand('input_data/gridarea_{model}.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
    
    output:
        outpath=outdir+'figs/AerChemMIP/depdust_{experiment}_Ayear_map.png' 
    wildcard_constraints:
        experiment="|".join(CONTROL_EXPS)
    
    notebook:
        "../notebooks/plot_emidust.py.ipynb"