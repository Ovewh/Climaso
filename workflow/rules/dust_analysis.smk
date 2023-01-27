
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
        absortion_plot = outdir+'figs/AerChemMIP/SWDirectEff_AAOD_refractive_index.png'

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
        variables = "dustMME|ssMME|bcMME|so4MME|dustMMEabs|bcMMEabs"
    notebook:
        "../notebooks/dust_analysis/MME.py.ipynb"

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
        