
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

        diag_table = outdir+'aerChemMIP_2xdust_table.csv'
    
    output: 
        absortion_plot = outdir+'figs/AerChemMIP/absorption_refractive_index.png'

    notebook:
        "../notebooks/dust_analysis/optical_properties_absorption.py.ipynb"


rule vertical_profiles:
    input:
        paths = expand(outdir + "{experiment}/loaddust/loaddust_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem', 
                    'MPI-ESM-1-2-HAM','MIROC6','UKESM1-0-LL',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1' ], allow_missing=True),

    output:
        path = outdir + "figs/AerChemMIP/vertical_profiles_{experiment}.png",

    conda:
        "../envs/comp_cat.yaml"
    
    notebook:
        "../notebooks/dust_analysis/vertical_profiles.py.ipynb"

rule calculate_lifetime:
    input:
        paths = expand(outdir + "{experiment}/concdust/concdust_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem', 
                    'MPI-ESM-1-2-HAM','MIROC6','GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1' ], allow_missing=True),
        paths_emissions = expand(outdir + "{experiment}/emidust/emidust_{experiment}_{model}_Ayear.nc",
                model=['EC-Earth3-AerChem', 
                    'MPI-ESM-1-2-HAM','MIROC6','GISS-E2-1-G',
                        'NorESM2-LM','GFDL-ESM4','CNRM-ESM2-1' ], allow_missing=True),
        path_area = "workflow/input_data/common_grid.nc"
    output:
        path = outdir + "figs/AerChemMIP/lifetime_{experiment}.png",
        table = outdir + "figs/AerChemMIP/lifetime_{experiment}.yaml"

    notebook:
        "../notebooks/dust_analysis/lifetime.py.ipynb"
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
        