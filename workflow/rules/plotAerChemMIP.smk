
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
        add_global_avg=True,
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
        vName ='atmabsLW|atmabsSW|atmabs'

    notebook:
        "../notebooks/plot_atmabs_AerChemMIP.py.ipynb"



rule plot_global_avaraged_ERFs:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc', 
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM',], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_globalaverage_AerChemMIP.png'
    wildcard_constraints:
        vName = 'ERFtsw|ERFtlw|ERFtcs|ERFt'

    notebook:
        "../notebooks/plot_ERFs_AerChemMIP_gobal_average.py.ipynb"

rule plot_albedo_radiative_effect:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc',
            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL', 'GFDL-ESM4','IPSL-CM6A-LR-INCA', 
            ],
            allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    wildcard_constraints:
        vName = 'ERFtswcsaf|ERFtcsaf|ERFtlwcsaf'
    params:
        draw_error_mask=False
    notebook:
        "../notebooks/plot_Radiative_effects_AerChemMIP.py.ipynb"

rule plot_direct_and_cloud_effect:
    input:
        paths=expand(outdir+'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc',
            model=['MPI-ESM-1-2-HAM','EC-Earth3-AerChem','CNRM-ESM2-1','NorESM2-LM','UKESM1-0-LL', 'GFDL-ESM4','IPSL-CM6A-LR-INCA',
            ],
            allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/ERFs/{vName}_piClim-2xdust_AerChemMIP.png'
    wildcard_constraints:
        vName = 'SWDirectEff|LWDirectEff|SWCloudEff|LWCloudEff|LWDirectEff_cs|SWDirectEff_cs|CloudEff|DirectEff'
    params:
        draw_error_mask=False

    notebook:
        "../notebooks/plot_Radiative_effects_AerChemMIP.py.ipynb"
    
rule plot_change_cdncvi:
    input:
        path_exp=expand(outdir+'piClim-2xdust/derived_variables/cdncvi/cdncvi_{model}_piClim-2xdust_Ayear.nc',
                model=['EC-Earth3-AerChem',  'GISS-E2-1-G', 
                        'MIROC6','GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/derived_variables/cdncvi/cdncvi_{model}_piClim-control_Ayear.nc',
                model=['EC-Earth3-AerChem',  'GISS-E2-1-G', 
                        'MIROC6','GFDL-ESM4',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    params:
        units= '# cm-2',
        label='CDNC',
        add_global_avg=True,
        rel_minmax=[-50,50]
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/cdncvi_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='rel|abs'

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

rule plot_change_lwp:
    input:
        path_exp = expand(outdir+'piClim-2xdust/lwp/lwp_piClim-2xdust_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem',
                        'UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'NorESM2.0.6dev-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/lwp/lwp_piClim-control_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem', 
                        'UKESM1-0-LL',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'NorESM2.0.6dev-LM'], allow_missing=True)
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

rule plot_change_clwvi:
    input:
        path_exp = expand(outdir+'piClim-2xdust/clwvi/clwvi_piClim-2xdust_{model}_Ayear.nc',
                model=['GISS-E2-1-G','EC-Earth3-AerChem',#'IPSL-CM6A-LR-INCA', 
                         'GFDL-ESM4', 'UKESM1-0-LL',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM']),
        path_ctrl=expand(outdir+'piClim-control/clwvi/clwvi_piClim-control_{model}_Ayear.nc',
                model=[ 'GISS-E2-1-G','EC-Earth3-AerChem',#'IPSL-CM6A-LR-INCA', 
                        'GFDL-ESM4', 'UKESM1-0-LL',  'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/clwvi_piClim-2xdust_AerChemMIP_{kind}.png'
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
                        'CNRM-ESM2-1','NorESM2-LM']),
        path_ctrl=expand(outdir+'piClim-control/clt/clt_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])
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

rule plot_change_aod:
    input:
        path_exp = expand(outdir+'piClim-2xdust/od550aer/od550aer_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA',
                         'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/od550aer/od550aer_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G',  'IPSL-CM6A-LR-INCA',
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/od550aer_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'
    params:
        label= '$\Delta$ Od550aer',
        rel_minmax=[-0.2,0.45],
        units='AOD',
        draw_error_mask=False,
        add_global_avg=True,
        abs_minmax=[-0.2,0.5]

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"


rule plot_change_aaod:
    input:
        path_exp = expand(outdir+'piClim-2xdust/abs550aer/abs550aer_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA',
                         'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/abs550aer/abs550aer_piClim-control_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G',  'IPSL-CM6A-LR-INCA',
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM', 'UKESM1-0-LL','MIROC6'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/abs550aer_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'
    params:
        label= '$\Delta$ AAOD 550mn',
        units='AOD',
        draw_error_mask=False,
        add_global_avg=True,
        abs_minmax=[-0.1,0.15],

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

rule plot_change_tas:
    input:
        path_exp = expand(outdir+'piClim-2xdust/tas/tas_piClim-2xdust_{model}_Ayear.nc',
                 model=['GISS-E2-1-G','MPI-ESM-1-2-HAM', 'IPSL-CM6A-LR-INCA',
                        'MIROC6', 'GFDL-ESM4', 'UKESM1-0-LL',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/tas/tas_piClim-control_{model}_Ayear.nc',
                 model=['GISS-E2-1-G', 'IPSL-CM6A-LR-INCA','MPI-ESM-1-2-HAM', 
                        'MIROC6', 'GFDL-ESM4', 'UKESM1-0-LL',   
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

rule plot_ccn_change:
    input:
        path_ctrl= expand(outdir+'piClim-control/ccn/ccn_piClim-control_{model}_Ayear.nc',
                 model=['MPI-ESM-1-2-HAM','NorESM2-LM'], allow_missing=True),
        path_exp = expand(outdir+'piClim-2xdust/ccn/ccn_piClim-2xdust_{model}_Ayear.nc',
                 model=['MPI-ESM-1-2-HAM','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/ccn_piClim-2xdust_AerChemMIP_{kind}.png'

    wildcard_constraints:
        kind='abs|rel'

    params:
        label = '$\Delta$ CCN',
        scaling_factor = 1e-6,
        units = 'cm$^{-3}$',


    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"

rule plot_change_prs:
    input:
        path_exp = expand(outdir+'piClim-2xdust/pr/pr_piClim-2xdust_{model}_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM']),
        path_ctrl=expand(outdir+'piClim-control/pr/pr_piClim-control_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 
                        'MIROC6','UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'])

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/pr_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs|rel'

    params:
        label='$\Delta$ Precipitation',
        rel_minmax=[-60,60],
        scaling_factor=1,
        units = "[mm year-1]",
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

rule plot_depdust:
    input:
        paths=expand(outdir+'{experiment}/depdust/depdust_{experiment}_{model}_Amon.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA',
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True),
        areacello = expand('workflow/input_data/gridarea_{model}.nc',
                    model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA',
                        'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
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
                        model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6' ,
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


rule plot_lwp_aerchemmip:
    input:
        expand(outdir+"{experiment}/lwp/lwp_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem',
                    'UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                    'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/lwp_{experiment}_map.png'
    params:
        vmin=0,
        vmax=90,
        cb_extend='max',
        label='LWP [g/m2]',
        nlevels=16,
        scaling_factor=1000,
        add_global_avg=True,
        time_slice=slice(3,-2)

    notebook:
        "../notebooks/plot_absolute_fields.py.ipynb"

rule plot_cli_aerchemmip:
    input:
        expand(outdir+"{experiment}/cli/cli_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem',
                    'UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                    'CNRM-ESM2-1','NorESM2-LM', 'GISS-E2-1-G'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/cli_{experiment}_{plevel}_map.png'
    params:
        vmin=0,
        vmax=0.025,
        cb_extend='max',
        label='cli [g/g]',
        nlevels=16,
        scaling_factor=1000,
        method='sum',
        add_global_avg=True,
        time_slice=slice(3,-2)
    wildcard_constraints:
        plevel='low|middle|high'

    conda:
        "../envs/comp_cat.yaml"
    notebook:
        "../notebooks/plot_absolute_fields_cloudlevels.py.ipynb"

rule plot_cl_aerchemmip:
    input:
        expand(outdir+"{experiment}/cl/cl_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem',
                    'UKESM1-0-LL', 'MPI-ESM-1-2-HAM',
                    'CNRM-ESM2-1','NorESM2-LM', 'GISS-E2-1-G'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/cl_{experiment}_{plevel}_map.png'
    params:
        vmin=0,
        vmax=70,
        cb_extend='max',
        label='cl %',
        nlevels=16,
        scaling_factor=1,
        method='sum',
        add_global_avg=True,
        time_slice=slice(3,-1),
        cmap = 'Blues'
    wildcard_constraints:
        plevel='low|middle|high'

    conda:
        "../envs/comp_cat.yaml"
    notebook:
        "../notebooks/plot_absolute_fields_cloudlevels.py.ipynb"


rule plot_pr_aerchemmip:
    input:
        expand(outdir+"{experiment}/pr/pr_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem', 'GISS-E2-1-G', 'IPSL-CM6A-LR-INCA', 'MIROC6',
                        'UKESM1-0-LL', 'GFDL-ESM4', 'MPI-ESM-1-2-HAM',
                        'CNRM-ESM2-1','NorESM2-LM'], allow_missing=True)
    output:
        outpath=outdir+'figs/AerChemMIP/pr_{experiment}_map.png'
    params:
        vmin=0,
        cb_extend='max',
        label='precipitation [kg/m2]',
        nlevels=16,
        add_global_avg=True,
        time_slice=slice(3,-2)

    notebook:
        "../notebooks/plot_absolute_fields.py.ipynb"

rule plot_level_changes:
    input:
        expand(rules.plot_level_cloud_changes.output, variable=['cli', 'cl'], plevel=['low','middle','high'],
            model=['GISS-E2-1-G','EC-Earth3-AerChem', 'MIROC6','CNRM-ESM2-1','NorESM2-LM','MPI-ESM-1-2-HAM', 'UKESM1-0-LL'])

rule plot_changes_dust_loading:
    input:
        expand(rules.plot_single_model_cloud_canges.output, 
                model=['EC-Earth3-AerChem', 
                    'MIROC6' ,
                    'MPI-ESM-1-2-HAM', 'UKESM1-0-LL',
                        'CNRM-ESM2-1','NorESM2-LM'], variable='loaddust', plevel=['low','middle','high'])

rule plot_change_concdust:
    input:
        path_exp = expand(outdir + "piClim-2xdust/derived_variables/concdust/concdust_{model}_piClim-2xdust_Ayear.nc",
                model=['EC-Earth3-AerChem','GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1', 'GFDL-ESM4',
                        'NorESM2-LM', 'UKESM1-0-LL', 'IPSL-CM6A-LR-INCA']),
        path_ctrl = expand(outdir + "piClim-control/derived_variables/concdust/concdust_{model}_piClim-control_Ayear.nc", 
                    model = ['EC-Earth3-AerChem', 'GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1', 'GFDL-ESM4',
                        'NorESM2-LM', 'UKESM1-0-LL', 'IPSL-CM6A-LR-INCA']),
        areacello = 'workflow/input_data/common_grid.nc'

    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/concdust_piClim-2xdust_AerChemMIP_{kind}.png'

    wildcard_constraints:
        kind='abs|rel'

    params:
        label = '$\Delta$ concdust',
        scaling_factor = 1e-9,
        units = '[Tg]',
        abs_minmax=[1e-5,1]


    notebook:
        "../notebooks/dust_analysis/plot_change_dustburden.py.ipynb"

rule plot_change_concso4:
    input:
        path_exp = expand(outdir + "piClim-2xdust/derived_variables/concso4/concso4_{model}_piClim-2xdust_Ayear.nc",
                model=['GISS-E2-1-G', 'MIROC6', 'MPI-ESM-1-2-HAM','EC-Earth3-AerChem',
                    'CNRM-ESM2-1', 'GFDL-ESM4','IPSL-CM6A-LR-INCA','UKESM1-0-LL',
                        'NorESM2-LM']),
        path_ctrl = expand(outdir + "piClim-control/derived_variables/concso4/concso4_{model}_piClim-control_Ayear.nc", 
                    model = ['GISS-E2-1-G', 'MIROC6', 'MPI-ESM-1-2-HAM','EC-Earth3-AerChem',
                     'CNRM-ESM2-1', 'GFDL-ESM4','IPSL-CM6A-LR-INCA','UKESM1-0-LL',
                        'NorESM2-LM']),
        areacello = 'workflow/input_data/common_grid.nc'

    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/concso4_piClim-2xdust_AerChemMIP_{kind}.png'

    wildcard_constraints:
        kind='abs|rel'

    params:
        label = '$\Delta$ SO4 burden',
        scaling_factor = 1e-9,
        units = '[Tg]',
        abs_minmax=[1e-5,1]


    notebook:
        "../notebooks/dust_analysis/plot_change_so4burden.py.ipynb"

rule plot_change_concpm1:
    input:
        path_exp = expand(outdir + "piClim-2xdust/derived_variables/concpm1/concpm1_{model}_piClim-2xdust_Ayear.nc",
                model=[ 'MPI-ESM-1-2-HAM'
                   , 'GFDL-ESM4',
                        'NorESM2-LM']),
        path_ctrl = expand(outdir + "piClim-control/derived_variables/concpm1/concpm1_{model}_piClim-control_Ayear.nc", 
                    model = [ 'MPI-ESM-1-2-HAM'
                     , 'GFDL-ESM4',
                        'NorESM2-LM']),
        areacello = 'workflow/input_data/common_grid.nc'

    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/concpm1_piClim-2xdust_AerChemMIP_{kind}.png'

    wildcard_constraints:
        kind='abs|rel'

    params:
        label = '$\Delta$ PM1 burden',
        scaling_factor = 1e-9,
        units = '[Tg]',
        abs_minmax=[1e-5,1]


    notebook:
        "../notebooks/dust_analysis/plot_change_pm1burden.py.ipynb"



rule plot_change_conch2oaer:
    input:
        path_exp = expand(outdir+'piClim-2xdust/derived_variables/conch2oaer/conch2oaer_{model}_piClim-2xdust_Ayear.nc',
                 model=['EC-Earth3-AerChem', 'IPSL-CM6A-LR-INCA',
                          'MPI-ESM-1-2-HAM','UKESM1-0-LL',
                       'NorESM2-LM', 'MIROC6'], allow_missing=True),
        path_ctrl=expand(outdir+'piClim-control/derived_variables/conch2oaer/conch2oaer_{model}_piClim-control_Ayear.nc',
                 model=['EC-Earth3-AerChem',  'IPSL-CM6A-LR-INCA',
                         'MPI-ESM-1-2-HAM','UKESM1-0-LL',
                        'NorESM2-LM', 'MIROC6'], allow_missing=True)

    
    output:
        outpath=outdir+'figs/AerChemMIP/delta_2xdust/conch2oaer_piClim-2xdust_AerChemMIP_{kind}.png'
    wildcard_constraints:
        kind='abs'
    params:
        # label= '$\Delta$ AAOD 550mn',
        # units='AOD',
        # draw_error_mask=False,
        add_global_avg=True,
        abs_minmax=[-1,1],
        log_norm=True,

    notebook:
        "../notebooks/plot_change_notebook.py.ipynb"


rule generate_text_table:
    input:
        clt_exp=rules.plot_change_clt.input.path_exp,
        clt_ctrl=rules.plot_change_clt.input.path_ctrl,
        lwp_ctrl=rules.plot_change_clwvi.input.path_ctrl,
        lwp_exp = rules.plot_change_clwvi.input.path_exp,
        prs_exp = rules.plot_change_prs.input.path_exp,
        prs_ctrl = rules.plot_change_prs.input.path_ctrl,
        clivi_exp = rules.plot_change_clivi.input.path_exp,
        clivi_ctrl = rules.plot_change_clivi.input.path_ctrl,
        od550aer_ctrl = rules.plot_change_aod.input.path_ctrl,
        od550aer_exp = rules.plot_change_aod.input.path_exp,
        abs550aer_ctrl = rules.plot_change_aaod.input.path_ctrl,
        abs550aer_exp = rules.plot_change_aaod.input.path_exp,
        concdust_ctrl = expand(outdir + "piClim-control/derived_variables/concdust/concdust_{model}_piClim-control_Ayear.nc", 
                    model = ['EC-Earth3-AerChem', 'GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'UKESM1-0-LL', 'CNRM-ESM2-1',
                        'NorESM2-LM', 'GFDL-ESM4']),
        concdust_exp = expand(outdir + "piClim-2xdust/derived_variables/concdust/concdust_{model}_piClim-2xdust_Ayear.nc",
                model=['EC-Earth3-AerChem','GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'UKESM1-0-LL',
                        'CNRM-ESM2-1','NorESM2-LM', 'GFDL-ESM4']),
        areacello=rules.plot_depdust.input.areacello,
        dulifetime=expand(outdir + "piClim-2xdust/dulifetime/lifetime_piClim-2xdust_dulifetime_{model}_Ayear.yaml",
                model=['EC-Earth3-AerChem','GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'IPSL-CM6A-LR-INCA',
                        'CNRM-ESM2-1','NorESM2-LM', 'GFDL-ESM4']),
        delta_MEEabs = expand(outdir + "piClim-2xdust/delta_dustMEEabs/dustMEEabs_piClim-2xdust_{model}_Ayear.yaml",
                model=['EC-Earth3-AerChem','GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1'
                       ,'NorESM2-LM', 'GFDL-ESM4']),
            
        delta_MEE = expand(outdir + "piClim-2xdust/delta_dustMEE/dustMEE_piClim-2xdust_{model}_Ayear.yaml",
                model=['EC-Earth3-AerChem','GISS-E2-1-G', 'MIROC6',
                    'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1'
                       ,'NorESM2-LM', 'GFDL-ESM4']),

        
        emidust_ctrl = expand(rules.plot_emidust.input.paths,zip, 
                    experiment=['piClim-control']),
        emidust_exp = expand(rules.plot_emidust.input.paths,zip, 
                    experiment=['piClim-2xdust']),
        ERFaci = expand(rules.plot_direct_and_cloud_effect.input.paths, zip, 
                    vName=['CloudEff','DirectEff','SWCloudEff', 'LWCloudEff','LWDirectEff','SWDirectEff']),
        ERFt = expand(rules.plot_ERFs.input.paths, zip,
                     vName=['ERFt','ERFtsw','ERFtlw']),
        atmabs = expand(rules.plot_atm_abs.input.paths, zip, 
                    vName=['atmabs','atmabsSW']),
        albedo = expand(rules.plot_albedo_radiative_effect.input.paths, zip, 
                    vName=['ERFtswcsaf','ERFtcsaf','ERFtlwcsaf']),

    log:
        "logs/generate_table.log"
    output:
        outpath=outdir+'aerChemMIP_2xdust_table.csv'
    
    notebook:
        "../notebooks/generate_text_table.py.ipynb"


rule generate_table:
    input:
        rules.generate_text_table.output
    log:
        "logs/generate_table.log"
    output:
        #outpath=outdir+'aerChemMIP_2xdust_table.csv',
        forcing_table=outdir+'forcing_table.png',
        diagnostics_table_path=outdir+'diagnostics_table.png',
        optics_diag_table = outdir+'optics_diag_table.png',
    
    notebook:
        "../notebooks/generate_table.py.ipynb"

rule boot_strap_sampling_dust_forcing:
    input:  
        rules.generate_text_table.output
    output:
        outpath = outdir +'boot_strapped_forcing_estimates.csv',
        outplot = outdir + 'boot_strapped_forcing_boxplot.png'

    notebook:
        "../notebooks/boot_strap_table.py.ipynb"