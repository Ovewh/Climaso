rule dust_variability:
    input:
        noresm = expand(rules.change_historical_ts.output.outpath, 
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                ext='csv', model=['NorESM2-LM']),
        ecEarth = expand(rules.change_historical_ts.output.outpath, 
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                ext='csv', model=['EC-Earth3-AerChem']),
        mpiEsm = expand(rules.change_historical_ts.output.outpath, 
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550aer','od550lt1aer'],
                ext='csv', model=['MPI-ESM-1-2-HAM'])
    output:
        forcing_variability_ts = outdir +'figs/antropogenic_forcing_variability/variability_ts.png',
        fraction_fine_mode = outdir + 'figs/antropogenic_forcing_variability/antropogenic_fine_mode_fraction.png',
        

    notebook:
        "../notebooks/antropogenic_forcing_variability_analysis.py.ipynb"
rule aod_variability_maps_forces_piClim:
    input:
        noresm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550dust','od550aer','od550lt1aer'],
                model='NorESM2-LM'),
        ecEarth = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550aer'],
                model='EC-Earth3-AerChem'),
        mpiEsm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550aer'],
                model='MPI-ESM-1-2-HAM')
    output:
        aod_varibility_map = outdir + 'figs/antropogenic_forcing_variability/variability_map_piClim.png'
    params:
        fig_title = 'Variability in Antropogenic \n AOD piClim-aer - piClim-control'

    notebook:
        "../notebooks/antropogenic_aod_variability_maps_piClim.py.ipynb"

# rule aod_variability_maps_forces_piClim4xCO2:
#     input:
#         noresm = expand(rules.calc_change_historical.output.outpath,
#                 hist_pert='abrupt-4xCO2', hist_base='piControl',
#                 vName = ['od550aer'],
#                 model='NorESM2-LM'),
#         ecEarth = expand(rules.calc_change_historical.output.outpath,
#                 hist_pert='abrupt-4xCO2', hist_base='piControl',
#                 vName = ['od550aer'],
#                 model='EC-Earth3-AerChem'),
#         mpiEsm = expand(rules.calc_change_historical.output.outpath,
#                 hist_pert='abrupt-4xCO2', hist_base='piControl',
#                 vName = ['od550aer'],
#                 model='MPI-ESM-1-2-HAM')
#     output:
#         aod_varibility_map = outdir + 'figs/antropogenic_forcing_variability/variability_map_4xCO2.png'
#     params:
#         fig_title = 'Variability in Antropogenic \n AOD abrupt-4xCO2 - piControl'

#     notebook:
#         "../notebooks/antropogenic_aod_variability_maps_piClim.py.ipynb"


rule aod_variability_maps_forces:
    input:
        noresm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                model='NorESM2-LM'),
        ecEarth = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                model='EC-Earth3-AerChem'),
        mpiEsm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550aer','od550lt1aer'],
                model='MPI-ESM-1-2-HAM')
    output:
        aod_varibility_map = outdir + 'figs/antropogenic_forcing_variability/variability_map.png'
    
    params:

    notebook:
        "../notebooks/antropogenic_aod_variability_maps.py.ipynb"


rule aod_variability_maps_forces_detrend:
    input:
        noresm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                model='NorESM2-LM'),
        ecEarth = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer'],
                model='EC-Earth3-AerChem'),
        mpiEsm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550aer','od550lt1aer'],
                model='MPI-ESM-1-2-HAM')
    output:
        aod_varibility_map = outdir + 'figs/antropogenic_forcing_variability/variability_map_detrended.png'
    
    params:
        detrend = True
    notebook:
        "../notebooks/antropogenic_aod_variability_maps.py.ipynb"

rule dust_wind_interactions_NorESM:
    input:
        expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                vName = ['od550dust','od550aer','od550lt1aer','sfcWind','emidust'],
                model='NorESM2-LM')
    output:
        wind_aod = outdir + 'figs/antropogenic_forcing_variability/aod_wind_dust.png'
    
    notebook:
        "../notebooks/wind_antropogenic_aod.py.ipynb"

rule dust_aod_variability:
    input:
        noresm_histSST = expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                model='NorESM2-LM', allow_missing=True),
        noresm_piClim = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                model='NorESM2-LM', allow_missing=True),
        ecEarth_histSST =  expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                model='EC-Earth3-AerChem', allow_missing=True),
        ecEarth_piClim =  expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                model='EC-Earth3-AerChem', allow_missing=True)
    output:
        dod_comp = outdir + 'figs/antropogenic_forcing_variability/variability_comparison_{vName}.png'
    notebook:
        "../compare_dod.py.ipynb"
        