
def get_variables(path):
    import os
    if os.path.exists(path):
        
        with open(path) as f:
            varName = yaml.safe_load(f)

    else: 
        with open('config/reference_data_request.yaml') as f:

            varName = yaml.safe_load(f)
    return varName
    
rule generate_text_table:
    input:
        'config/.data_trackers/piClim-2xdust_{model}_CMIP6.yaml',
        'config/.data_trackers/piClim-control_{model}_CMIP6.yaml',
        experiments_diags = lambda w: expand(outdir+f"piClim-2xdust/{{variable}}/{{variable}}_piClim-2xdust_{w.model}_Ayear.nc",
                        variable=get_variables(f"config/.data_trackers/piClim-2xdust_{w.model}_CMIP6.yaml")['vars'], 
                                allow_missing=True),
        control_diags = lambda w: expand(outdir+f'piClim-control/{{variable}}/{{variable}}_piClim-control_{w.model}_Ayear.nc',
                        variable=get_variables(f"config/.data_trackers/piClim-control_{w.model}_CMIP6.yaml")['vars'], 
                        allow_missing=True),

        exp_derived_diag = lambda w: expand(outdir+f'piClim-2xdust/derived_variables/{{variable}}/{{variable}}_{w.model}_piClim-2xdust_Ayear.nc',
                        variable=get_variables(f"config/.data_trackers/piClim-2xdust_{w.model}_CMIP6.yaml")['derived_vars'],
                         allow_missing=True),
        
        control_derived_diag = lambda w: expand(outdir+f'piClim-control/derived_variables/{{variable}}/{{variable}}_{w.model}_piClim-control_Ayear.nc',
                                variable= get_variables(f"config/.data_trackers/piClim-control_{w.model}_CMIP6.yaml")['derived_vars'], allow_missing=True),
        erfs = lambda w:  expand(outdir+f'piClim-2xdust/ERFs/{{variable}}/{{variable}}_piClim-2xdust_{w.model}_Ayear.nc',
                        variable=get_variables(f"config/.data_trackers/piClim-2xdust_{w.model}_CMIP6.yaml")['ERFs'], allow_missing=True),
        areacello='workflow/input_data/common_grid.nc',

    log:
        "logs/generate_text_table/aerChemMIP_piClim_2xdust_table_{model}.log"
    params:
        time_slice = slice(3,None)

    output:
        metadata=outdir+'diagnostics_piClim_2xdust/{model}/metadata_table_{model}.csv',
        diff=outdir+'diagnostics_piClim_2xdust/{model}/diff_{model}.csv',
        exp=outdir+'diagnostics_piClim_2xdust/{model}/exp_{model}.csv',
        ctrl=outdir+'diagnostics_piClim_2xdust/{model}/ctrl_{model}.csv',
        info_yaml=outdir+'diagnostics_piClim_2xdust/{model}/analysis_info_{model}.yaml',
    
    notebook:
        "../notebooks/generate_text_table.py.ipynb"


rule plot_diagnostics_table:
    input:
        diff=expand(outdir+'diagnostics_piClim_2xdust/{model}/diff_{model}.csv', model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA','CNRM-ESM2-1', 'GFDL-ESM4']),
        exp=expand(outdir+'diagnostics_piClim_2xdust/{model}/metadata_table_{model}.csv',
                 model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA','CNRM-ESM2-1', 'GFDL-ESM4']),
    log:
        "logs/generate_table.log"
    output:
        #outpath=outdir+'aerChemMIP_2xdust_table.csv',
        diagnostics_table_path=outdir+'cloud_diagnostics.pdf',
        optics_diag_table = outdir+'optics_diag.pdf',
    notebook:
        "../notebooks/generate_table.py.ipynb"

rule plot_forcing_decomposition_paper:
    input:
        expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        outdir+'figs/AerChemMIP/dust_radiative_effects.pdf'
    notebook:
        "../notebooks/dust_analysis/forcing_plot_paper.py.ipynb"    

rule plot_forcing_decomposition_cacti:
    input:
        expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        directory(outdir+'figs/AerChemMIP/ERFfigures/')
    notebook:
        "../notebooks/dust_analysis/forcing_plot.py.ipynb"