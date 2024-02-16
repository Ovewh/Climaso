
    
rule text_diagnostic_table:
    input: 
        ctrl_data = expand(outdir + 'dust_diag_files/dust_diag_{model}_piClim-control.nc',
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4','CNRM-ESM2-1']),
        exp_data = expand(outdir + 'dust_diag_files/dust_diag_{model}_piClim-2xdust.nc',
                model=['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4', 'CNRM-ESM2-1']),
        mask = outdir + 'masks/dust_regions.nc',

    output:
        outpath=outdir+'tables/AerChemMIP/dust_abs_diagnostic_table.csv',
        relpath=outdir+'tables/AerChemMIP/dust_rel_diagnostic_table.csv'
    notebook:
        "../notebooks/dust_analysis/dust_diagnostic_table.py.ipynb"        

rule plot_forcing_decomposition_paper:
    input:
        forcing_tables = expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        diag_tables = rules.text_diagnostic_table.output

    output:
        outdir+'figs/AerChemMIP/dust_radiative_effects.pdf'
    notebook:
        "../notebooks/dust_analysis/forcing_plot_paper.py.ipynb"    


rule plot_forcing_effciency:
    input:
        forcing_tables = expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4']),
        diag_tables = rules.text_diagnostic_table.output

    output:
        outdir+'figs/AerChemMIP/dust_radiative_efficiency.pdf'
    notebook:
        "../notebooks/dust_analysis/forcing_efficiency_plot.py.ipynb"

rule plot_forcing_decomposition_cacti:
    input:
        expand(outdir + 'piClim-2xdust/ERFs/ERF_tables/piClim-2xdust_{model}.csv',
        model = ['NorESM2-LM', 'MPI-ESM-1-2-HAM', 'CNRM-ESM2-1','EC-Earth3-AerChem', 'GISS-E2-1-G',
                        'UKESM1-0-LL', 'MIROC6', 'IPSL-CM6A-LR-INCA', 'GFDL-ESM4'])

    output:
        directory(outdir+'figs/AerChemMIP/ERFfigures/')
    notebook:
        "../notebooks/dust_analysis/forcing_plot.py.ipynb"