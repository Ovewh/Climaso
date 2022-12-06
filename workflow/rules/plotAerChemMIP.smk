
rule plot_ERFs:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc', 
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    
    wildcard_constraints:
        vName = 'ERFtsw|ERFtlw|ERFtcs|ERFt|ERFsurf|ERFsurfsw'
    params:
        draw_error_mask=False,
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
            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL', 'GFDL-ESM4','IPSL-CM6A-LR-INCA'],
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
    params:
        scaling_factor=1e-6,
        units= '# cm-3',
        label='CDNC',
        rel_minmax=[-50,50]
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/cdnc_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

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

    params:
        scaling_factor=1e3,
        units="g m-2",
        rel_minmax=[-30,50],
        label="Liquid water Path"

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

rule plot_change_clivi:
    input:  
        path_exp = expand(outdir+'piClim-2xdust/clivi/clivi_piClim-2xdust_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM']),
        path_ctrl = expand(outdir + 'piClim-control/clivi/clivi_piClim-control_{model}_Ayear.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/clivi_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'
    params:
        scaling_factor=1000,
        units= '[g m-2]',
        label='Ice water path',
        abs_minmax=[-10,20],
        rel_minmax=[-50,50]

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

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
    params:
        label='Cloud fraction',
        abs_minmax=[-8,8],
        rel_minmax=[-15,15],
        draw_error_mask=True
    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

rule plot_change_ts:
    input:
        path_exp = expand(outdir+'piClim-2xdust/ts/ts_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 
                         'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/ts/ts_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G',  
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/ts_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'
    params:
        label= '$\Delta$ T',
        rel_minmax=[-30,50],
        units=['K'],
        draw_error_mask=False,
        abs_minmax=[-0.5,0.5]

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

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

    params:
        label= '$\Delta$ T',
        rel_minmax=[-30,50],
        units=['K']

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"


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

    params:
        label='$\Delta$ Precipitation',
        abs_minmax=[-0.3,0.3],
        rel_minmax=[-60,60],
        scaling_factor=1000,
        units = "[g m-2 s-1]",
        cmap='BrBG',
        draw_error_mask=True


    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"


rule plot_emidust:
    input:
        paths=expand(outdir+'{experiment}/emidust/emidust_{experiment}_{model}_Amon.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        areacello = expand('workflow/input_data/gridarea_{model}.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
            
    output:
        outpath=outdir+'figs/AerChemMIP/emidust_{experiment}_Ayear_map.png' 
    wildcard_constraints:
        experiment="|".join(ALL_ALLOWED)
    
    notebook:
        "../notebooks/plot_emidust.py.ipynb"

rule plot_feedback_decomposed:
    input:
        paths_ERFt=expand(outdir + '{experiment}/Feedback_per_emis/ERFt_{variable}_{experiment}_{model}_Ayear.yaml',
                    model=[ 'GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        paths_CloudEff = expand(outdir + '{experiment}/Feedback_per_emis/CloudEff_{variable}_{experiment}_{model}_Ayear.yaml',
                            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL', 'GFDL-ESM4'],
                            allow_missing=True),
        paths_DirectEff = expand(outdir + '{experiment}/Feedback_per_emis/DirectEff_{variable}_{experiment}_{model}_Ayear.yaml',
                            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL', 'GFDL-ESM4'],
                            allow_missing=True)

    output:
        outpath=outdir+'figs/AerChemMIP/Feedbacks/Feedback_per_emis_decomposed_{variable}_{experiment}.png'

    notebook:
        "../notebooks/plot_feedback.py.ipynb"
rule plot_depdust:
    input:
        paths=expand(outdir+'{experiment}/depdust/depdust_{experiment}_{model}_Amon.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        areacello = expand('workflow/input_data/gridarea_{model}.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
    params:
        regrid=False

    output:
        outpath=outdir+'figs/AerChemMIP/depdust_{experiment}_Ayear_map.png' 
    wildcard_constraints:
        experiment="|".join(CONTROL_EXPS)
    
    notebook:
        "../notebooks/plot_emidust.py.ipynb"

rule plot_level_cloud_changes:
    input:
        paths = expand(outdir+'piClim-2xdust/delta_{variable}/delta_{variable}_piClim-2xdust_{model}.nc',
                        model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6' ,'IPSL-CM6A-LR-INCA',
                        'MPI-ESM-1-2-HAM', 
                        'UKESM1-0-LL',
                            
                        'CNRM-ESM2-1','NorESM2-LM'],allow_missing=True)
    output:
        outpath = outdir + 'figs/AerChemMIP/delta_{variable}_{plevel}_piClim-2xdust.png'

    wildcard_constraints:
        plevel='low|middle|high'


    notebook:
        "../notebooks/plot_level_cloud_change.py.ipynb"

rule plot_single_model_cloud_canges:
    input:
        path = outdir+'piClim-2xdust/delta_{variable}/delta_{variable}_piClim-2xdust_{model}.nc'
    output:
        outpath = outdir+'figs/AerChemMIP/delta_2xdust/single_model_plots/delta_{variable}_{plevel}_piClim-2xdust_{model}.png'

    wildcard_constraints:
        plevel='low|middle|high'
    
    notebook:
        "../notebooks/plot_level_cloud_change_single_model.py.ipynb"



rule plot_level_changes:
    input:
        expand(rules.plot_level_cloud_changes.output, variable=['cli', 'cl'], plevel=['low','middle','high'],
            model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6','CNRM-ESM2-1','NorESM2-LM'])

# rule plot_all_AerChemMIP:
#     input:
#         expand(rules.plot_level_changes.output, 
#                 model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6' ,'IPSL-CM6A-LR-INCA',
#                         'CNRM-ESM2-1','NorESM2-LM'], variable='mmrdust', plevel=['low','middle','high'])
#         expand(rules.plot_level_cloud_changes.output, variable=['cli', 'cl'], plevel=['low','middle','high'],
#             model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6','CNRM-ESM2-1','NorESM2-LM', 'MPI-ESM-1-2-HAM'
#             ,'IPSL-CM6A-LR-INCA','UKESM1-0-LL'])

rule generate_table:
    input:
        clt_exp=rules.plot_change_clt.input.path_exp,
        clt_ctrl=rules.plot_change_clt.input.path_ctrl,
        lwp_ctrl=rules.plot_change_lwp.input.path_ctrl,
        lwp_exp = rules.plot_change_lwp.input.path_exp,
        prs_exp = rules.plot_change_prs.input.path_exp,
        prs_ctrl = rules.plot_change_prs.input.path_ctrl,
        clivi_exp = rules.plot_change_clivi.input.path_exp,
        clivi_ctrl = rules.plot_change_clivi.input.path_ctrl,
        areacello=rules.plot_depdust.input.areacello,
        emidust_ctrl = expand(rules.plot_emidust.input.paths,zip, 
                    experiment=['piClim-control']),
        emidust_exp = expand(rules.plot_emidust.input.paths,zip, 
                    experiment=['piClim-2xdust']),
        ERFaci = expand(rules.plot_ERFaci.input.paths, zip, 
                    vName=['CloudEff','DirectEff','SWCloudEff', 'LWCloudEff','LWDirectEff','SWDirectEff']),
        ERFt = expand(rules.plot_ERFs.input.paths, zip,
                     vName=['ERFt','ERFtsw','ERFtlw']),
        atmabs = expand(rules.plot_atm_abs.input.paths, zip, 
                    vName=['atmabs']),
        feedback_tot = expand(rules.plot_feedback_decomposed.input.paths_ERFt, 
                                experiment=['piClim-2xdust'], variable=['emidust']),
        feedback_Clouds = expand(rules.plot_feedback_decomposed.input.paths_CloudEff, 
                                experiment=['piClim-2xdust'], variable=['emidust']),
        feedback_Direct = expand(rules.plot_feedback_decomposed.input.paths_DirectEff, 
                                experiment=['piClim-2xdust'], variable=['emidust'])
    log:
        "logs/generate_table.log"
    output:
        outpath=outdir+'aerChemMIP_2xdust_table.csv',
        forcing_table=outdir+'forcing_table.png',
        diagnostics_table_path=outdir+'diagnostics_table.png',
        feedback_table = outdir+'feedback_table.png'
    
    notebook:
        "../notebooks/generate_table.py.ipynb"

rule boot_strap_sampling_dust_forcing:
    input:  
        rules.generate_table.output.outpath
    output:
        outpath = outdir +'boot_strapped_forcing_estimates.csv',
        outplot = outdir + 'boot_strapped_forcing_boxplot.png'

    notebook:
        "../notebooks/boot_strap_table.py.ipynb"