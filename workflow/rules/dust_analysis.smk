
rule refractive_index_and_absorption:
    input:
        atmabs = expand(rules.plot_atm_abs.input.paths, zip, 
                    vName=['atmabs']),
        abs550_ctrl = expand(outdir + 'piClim-control/abs550aer/abs550aer_piClim-control_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G','MIROC6','GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                     'CNRM-ESM2-1','MPI-ESM-1-2-HAM','EC-Earth3-AerChem','UKESM1-0-LL']),
        
        abs550_exp = expand(outdir + 'piClim-2xdust/abs550aer/abs550aer_piClim-2xdust_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G', 'MIROC6','GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                    'CNRM-ESM2-1','MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem','UKESM1-0-LL']),

        oddust550_ctrl = expand(outdir + 'piClim-control/od550dust/od550dust_piClim-control_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G','GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                     'CNRM-ESM2-1','MPI-ESM-1-2-HAM','EC-Earth3-AerChem','UKESM1-0-LL']),
        
        oddust550_exp = expand(outdir + 'piClim-2xdust/od550dust/od550dust_piClim-2xdust_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G', 'GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                    'CNRM-ESM2-1','MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem','UKESM1-0-LL']),

        diag_table = expand(outdir+'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G', 'GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                    'CNRM-ESM2-1','MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem','UKESM1-0-LL'])
    
    output: 
        absortion_plot = outdir+'figs/AerChemMIP/SWDirectEff_AAOD_refractive_index.png'

    notebook:
        "../notebooks/dust_analysis/optical_properties_absorption.py.ipynb"


rule vertical_profiles:
    input:
        paths = expand(outdir + "piClim-control/loaddust/loaddust_piClim-control_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem', 
                    'MPI-ESM-1-2-HAM','MIROC6',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1' ], allow_missing=True),

    output:
        path = outdir + "figs/AerChemMIP/dust_vertical_profiles.png",

    conda:
        "../envs/comp_cat.yaml"
    
    notebook:
        "../notebooks/dust_analysis/vertical_profiles.py.ipynb"

def variable_translator_lifetime(w):
    variable = w.variable
    vdict = {
        'dulifetime' : ['concdust','depdust'],
        'so4lifetime' : ['concso4','depso4'],
        'concss' : ['concss', 'depss'],
    }
    return vdict[variable]


rule calculate_lifetime:
    input:
        paths = lambda w: expand(outdir + f"{{experiment}}/derived_variables/{variable_translator_lifetime(w)[0]}/{variable_translator_lifetime(w)[0]}_{{model}}_{{experiment}}_Ayear.nc", 
                   allow_missing=True),
        paths_emissions = lambda w: expand(outdir + f"{{experiment}}/{variable_translator_lifetime(w)[1]}/{variable_translator_lifetime(w)[1]}_{{experiment}}_{{model}}_Ayear.nc",
                allow_missing=True),
        path_area = "workflow/input_data/common_grid.nc"
    output:
        table = outdir + "{experiment}/{variable}/lifetime_{experiment}_{variable}_{model}_Ayear.yaml"
    wildcard_constraints:
        variables = "dulifetime|so4lifetime|pm1lifetime|concss"

    notebook:
        "../notebooks/dust_analysis/lifetime.py.ipynb"

rule aerosol_species_lifetime_changes:
    input:
        du = expand(outdir+"piClim-2xdust/dulifetime/lifetime_piClim-2xdust_dulifetime_{model}_Ayear.yaml", 
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','MIROC6', 'GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1', 'IPSL-CM6A-LR-INCA']),
        du_ctrl = expand(outdir+"piClim-control/dulifetime/lifetime_piClim-control_dulifetime_{model}_Ayear.yaml", 
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','MIROC6', 'GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1', 'IPSL-CM6A-LR-INCA']),
        so = expand(outdir+"piClim-2xdust/so4lifetime/lifetime_piClim-2xdust_so4lifetime_{model}_Ayear.yaml", 
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','MIROC6', 'GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1', 'IPSL-CM6A-LR-INCA']),
        so_ctrl = expand(outdir+"piClim-control/so4lifetime/lifetime_piClim-control_so4lifetime_{model}_Ayear.yaml", 
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','MIROC6', 'GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1', 'IPSL-CM6A-LR-INCA']),
    output: 
        outpath = outdir + "figs/AerChemMIP/lifetimechanges_piClim-2xdust.png"
    notebook:
        "../notebooks/dust_analysis/lifetime_changes.py.ipynb"

def variable_translator_MEE(w):
    variable = w.variable
    vdict = {
        'dustMEE' : ['concdust','od550aer'],
        'ssMEE' : ['concss','od550aer'],
        'bcMEE' : ['concbc','od550aer'],
        'so4MEE' : ['concso4','od550aer'],
        'dustMEEabs': ['concdust','abs550aer'],
        'bcMEEabs': ['concbc','abs550aer']
    }
    return vdict[variable]

rule calculate_MEE_change:
    input:
        burden = lambda w: outdir + f"{{experiment}}/derived_variables/{variable_translator_MEE(w)[0]}/{variable_translator_MEE(w)[0]}_{{model}}_{{experiment}}_Ayear.nc",
        od550aer = lambda w: outdir + f"{{experiment}}/{variable_translator_MEE(w)[1]}/{variable_translator_MEE(w)[1]}_{{experiment}}_{{model}}_Ayear.nc",
        burden_ctrl = lambda w: outdir + f"piClim-control/derived_variables/{variable_translator_MEE(w)[0]}/{variable_translator_MEE(w)[0]}_{{model}}_piClim-control_Ayear.nc",
        od550aer_ctrl = lambda w: outdir + f"piClim-control/{variable_translator_MEE(w)[1]}/{variable_translator_MEE(w)[1]}_piClim-control_{{model}}_Ayear.nc"
    output:
        outpath = outdir + "{experiment}/delta_{variable}/{variable}_{experiment}_{model}_Ayear.yaml"
    wildcard_constraints:
        variables = "dustMEE|ssMEE|bcMEE|so4MEE|dustMEEabs|bcMEEabs"
    notebook:
        "../notebooks/dust_analysis/MME.py.ipynb"

rule make_dust_diag_file:
    input:
        catalog = ancient(rules.build_catalogues.output.json),
        mask = outdir + 'masks/dust_regions.nc',
        burden_exp = outdir + '{experiment}/derived_variables/concdust/concdust_{model}_{experiment}_Ayear.nc',
        universial_area_mask = 'workflow/input_data/common_grid.nc', 
        model_area_mask = 'workflow/input_data/gridarea_{model}.nc'
    output:
        dust_diag_exp = outdir + 'dust_diag_files/dust_diag_{model}_{experiment}.nc',
    
    wildcard_constraints:
        experiment = 'piClim-2xdust|piClim-control'
    
    notebook:
        "../notebooks/dust_analysis/make_dust_diag_file.py.ipynb"

rule calc_dust_regional_erf_table:
    input:
        catalog = ancient(rules.build_catalogues.output.json),
        data_tracker = ancient('config/.data_trackers/piClim-2xdust_{model}_CMIP6.yaml'),
        mask = outdir + 'masks/dust_regions.nc'
    output:
        outpath_masked = outdir + 'piClim-2xdust/ERFs/ERF_tables/dusty/piClim-2xdust_{model}.csv',
        outpath_unmasked = outdir + 'piClim-2xdust/ERFs/ERF_tables/nodusty/piClim-2xdust_{model}.csv',
        outpath_all = outdir + 'piClim-2xdust/ERFs/ERF_tables/all/piClim-2xdust_{model}.csv'
    
    log:
        "logs/erf_tables/{model}_piClim-2xdust_regional.log"
    notebook:
        "../notebooks/forcing_calculations/calc_global_regional_erf.py.ipynb"

rule plot_diagnostic_table:
    input:
        ctrl_data = expand(outdir + 'dust_diag_files/dust_diag_{model}_piClim-control.nc',
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        exp_data = expand(outdir + 'dust_diag_files/dust_diag_{model}_piClim-2xdust.nc',
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        mask = outdir + 'masks/dust_regions.nc',
    output:    
        outpath = outdir + 'figs/AerChemMIP/dust_diagnostic_table.pdf'

    notebook:
        "../notebooks/dust_analysis/dust_diagnostic_table.py.ipynb"

rule plot_forcing_decomposition_nodust:
    input:
        expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/nodusty/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        directory(outdir+'figs/AerChemMIP/ERFfigures/nodust/')
    notebook:
        "../notebooks/dust_analysis/forcing_plot.py.ipynb"



rule plot_forcing_decomposition_dusty:
    input:
        expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/dusty/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        directory(outdir+'figs/AerChemMIP/ERFfigures/dust/')
    notebook:
        "../notebooks/dust_analysis/forcing_plot.py.ipynb"


rule plot_dusty_vs_no_dusty_changes:
    input:
        dusty_region=expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/dusty/piClim-2xdust_{model}.csv',
                model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4', 'CNRM-ESM2-1' ]),
        nodusty_region=expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/nodusty/piClim-2xdust_{model}.csv',
                model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4','CNRM-ESM2-1']),
        all_region = expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/nodusty/piClim-2xdust_{model}.csv',
                model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4', 'CNRM-ESM2-1' ]),
    output:
        outpath = outdir+'figs/AerChemMIP/dusty_vs_non_dusty.png'
    
    notebook:
        "../notebooks/dust_analysis/dusty_vs_non_dusty.py.ipynb"

rule create_diagnotics_table:
    input:
        catalog = ancient(rules.build_catalogues.output.json),
        data_tracker = ancient('config/.data_trackers/piClim-2xdust_{model}_CMIP6.yaml'),
        mask = outdir + 'masks/dust_regions.nc',
        burden = expand(outdir+"piClim-2xdust/derived_variables/concdust/concdust_{model}_piClim-2xdust_Ayear.nc", 
                                        model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        burden_ctrl = expand( outdir + "piClim-control/derived_variables/concdust/concdust_{model}_piClim-control_Ayear.nc",
                            model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        
    


rule what_does_2xdust_mean:
    input:
        oddust550_ctrl = expand(outdir + 'piClim-control/od550dust/od550dust_piClim-control_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G','GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                     'CNRM-ESM2-1','MPI-ESM-1-2-HAM','EC-Earth3-AerChem','UKESM1-0-LL']),
        
        oddust550_exp = expand(outdir + 'piClim-2xdust/od550dust/od550dust_piClim-2xdust_{model}_Ayear.nc', 
                    model=['GFDL-ESM4', 'NorESM2-LM', 'GISS-E2-1-G', 'GISS-E2-1-G','IPSL-CM6A-LR-INCA',
                    'CNRM-ESM2-1','MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem','UKESM1-0-LL'])
    output:
        path = outdir + "figs/AerChemMIP/what_does_2xdust_mean.png",
    
    notebook:
        "../notebooks/dust_analysis/what_does_2xdust_mean.py.ipynb"    


rule plot_dust_radiation_interactions:
    input:
        catalog = ancient(rules.build_catalogues.output.json),
        mask = outdir + 'masks/dust_regions.nc',
        burden = expand(outdir+"piClim-2xdust/derived_variables/concdust/concdust_{model}_piClim-2xdust_Ayear.nc", 
                                        model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        od550aer = expand(outdir + "piClim-2xdust/od550aer/od550aer_piClim-2xdust_{model}_Ayear.nc", 
                        model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        
        burden_ctrl = expand( outdir + "piClim-control/derived_variables/concdust/concdust_{model}_piClim-control_Ayear.nc",
                            model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        od550aer_ctrl = expand(outdir + "piClim-control/od550aer/od550aer_piClim-control_{model}_Ayear.nc", 
                          model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        outdir+'figs/AerChemMIP/dust_radiation_interactions.pdf'

    params:
        models = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']
    notebook:
        "../notebooks/dust_analysis/dust_radiation_interactions.py.ipynb"    
