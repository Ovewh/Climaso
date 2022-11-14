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
        "../notebooks/compare_dod.py.ipynb"


rule regional_forcing_trends_and_variability:
    input:
        histSST=expand(outdir+'histSST/ERFs/{vName}/{vName}_histSST_{model}_Ayear.nc',
        model=['EC-Earth3-AerChem','NorESM2-LM','MPI-ESM-1-2-HAM'],vName=['ERFt','DirectEff','CloudEff']),
        piClim_aer = expand(outdir+'piClim-aer/ERFs/{vName}/{vName}_piClim-aer_{model}_Ayear.nc',
        model=['EC-Earth3-AerChem','NorESM2-LM','MPI-ESM-1-2-HAM'],vName=['ERFt'])

    output:
        outpath=outdir + 'figs/antropogenic_forcing_variability/forcing_variability_comparison_trends.png'

    notebook:
        "../notebooks/regional_trends_and_variability.py.ipynb" 


rule regional_trends_and_variability:
    input:
        expand(rules.calc_change_historical.output.outpath,
                hist_pert='histSST', hist_base='histSST-piAer',
                model=['NorESM2-LM', 'EC-Earth3-AerChem', 'MPI-ESM-1-2-HAM'], allow_missing=True)
    output:
        outpath=outdir + 'figs/antropogenic_forcing_variability/{vName}_variability_comparison_trends.png'

    notebook:
        "../notebooks/regional_trends_and_variability.py.ipynb"        
      
       
rule compare_forcing_variabiliy:
    input:
        piCli2xdust=expand(outdir + 'piClim-2xdust/ERFs/{vName}/{vName}_piClim-2xdust_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','NorESM2-LM'], allow_missing=True),
        piCliaer = expand(outdir + 'piClim-aer/ERFs/{vName}/{vName}_piClim-aer_{model}_Ayear.nc',
                model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','NorESM2-LM'], allow_missing=True),
    output:
        outpath = outdir + 'figs/antropogenic_forcing_variability/piClimAer-piClim2dust_{vName}_variability_comp.png'

    notebook:
        "../notebooks/forcing_variability_comparison_piCliaer_piClidust.py.ipynb"

rule piClim_interannual:
    input:
        picli_vars = expand(outdir + 'piClim-control/{vName}/{vName}_piClim-control_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut'], 
                    model=['NorESM2-LM','MPI-ESM-1-2-HAM','EC-Earth3-AerChem']),
        piaer_vars = expand(outdir + 'piClim-aer/{vName}/{vName}_piClim-aer_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut'], 
                    model=['NorESM2-LM','MPI-ESM-1-2-HAM','EC-Earth3-AerChem'])
    output: 
        outpath= outdir + 'figs/antropogenic_forcing_variability/toa_variability_picli_piaer.png'
    notebook:
        "../notebooks/inter_annual_variability/TOA_variability_comparison_piCliaer_piCli.py.ipynb"

rule toa_imbalance_dust_correlations:
    input:
        picli_rad = expand(outdir + 'piClim-control/{vName}/{vName}_piClim-control_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','rlutcs','rsutcs'], 
                    model=['NorESM2-LM','EC-Earth3-AerChem','MPI-ESM-1-2-HAM']),
        piaer_rad = expand(outdir + 'piClim-aer/{vName}/{vName}_piClim-aer_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','rlutcs','rsutcs'], 
                    model=['NorESM2-LM','EC-Earth3-AerChem','MPI-ESM-1-2-HAM']),
        piaer_load = expand(outdir + 'piClim-aer/{vName}/{vName}_piClim-aer_{model}_Ayear.nc', 
                        vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM', 'NorESM2-LM']),
        picli_load = expand(outdir + 'piClim-control/{vName}/{vName}_piClim-control_{model}_Ayear.nc',
                        vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM','NorESM2-LM']),
    output: 
        picli_allsky= outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_picli_allsky.png',
        picli_clearsky = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_picli_clearsky.png',
        piaer_allsky  = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_piaer_allsky.png',
        piaer_clearsky = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_piaer_clearsky.png'
    notebook:
        "../notebooks/inter_annual_variability/TOA_variability_and_dust_piCliaer_piCli.py.ipynb"

rule toa_imbalance_dust_correlations_histSST:
    input:
        histSST_rad = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','rlutcs','rsutcs'], 
                    model=['NorESM2-LM','EC-Earth3-AerChem','MPI-ESM-1-2-HAM']),
        histSST_load = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                        vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM', 'NorESM2-LM']),

    output: 
        histSST_allsky_end= outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_histSST_allsky_1985-2015.png',
        histSST_clearsky_end = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_histSST_clearsky_1985-2015.png',
        histSST_allsky_beginning  = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_histSST_allsky_1860-1890.png',
        histSST_clearsky_beginning = outdir + 'figs/antropogenic_forcing_variability/toa_variability_and_dust_histSST_clearsky_1860-1890.png'
    notebook:
        "../notebooks/inter_annual_variability/TOA_variability_and_dust_histSST.py.ipynb"

rule toa_dustload_relationship:
    input:
        noresm_rad = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc',
                            vName=['rsdt','rlut','rlutcs','rsutcs'],
                            experiment=['piClim-control','piClim-aer', 'piClim-2xdust','piClim-2xss','piClim-lu','piClim-SO2', 'piClim-BC']),
        noresm_load = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc',
                        vName=['mmrdust','airmass'], experiment=['piClim-control','piClim-aer', 'piClim-2xdust', 
                                'piClim-2xss','piClim-lu','piClim-SO2', 'piClim-BC'] ),     
        
        ec_rad = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_EC-Earth3-AerChem_Ayear.nc',
                            vName=['rsdt','rlutcs','rsutcs'],
                            experiment=['piClim-control','piClim-aer', 'piClim-2xdust','piClim-CH4', 'piClim-NTCF']),    
        ec_load = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_EC-Earth3-AerChem_Ayear.nc',
                        vName=['mmrdust','airmass'], experiment=['piClim-control','piClim-aer', 'piClim-2xdust','piClim-CH4', 'piClim-NTCF']),        
        
        
        mpi_load = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_MPI-ESM-1-2-HAM_Ayear.nc',
                        vName=['mmrdust','airmass'], experiment=['piClim-control','piClim-aer', 'piClim-2xdust','piClim-2xss', 
                                                                'piClim-OC','piClim-SO2']),  
        
        mpi_rad = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_MPI-ESM-1-2-HAM_Ayear.nc',
                            vName=['rsdt','rlutcs','rsutcs'],
                            experiment=['piClim-control','piClim-aer', 'piClim-2xdust', 'piClim-2xss', 
                                'piClim-OC','piClim-SO2']),    

        histSST_rad = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                    vName=['rsdt','rlutcs','rsutcs'], 
                    model=['NorESM2-LM','EC-Earth3-AerChem','MPI-ESM-1-2-HAM']),
        histSST_load = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                        vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM', 'NorESM2-LM']),
        
    output:
        toa_dustload_txt = outdir + 'figs/antropogenic_forcing_variability/dustload_toa_imbalance.csv',
        toa_dustload_png = outdir + 'figs/antropogenic_forcing_variability/dustload_toa_imbalance.png',
        toa_dustload_sw_png = outdir + 'figs/antropogenic_forcing_variability/dustload_toa_sw_imbalance.png'

    notebook:
        "../notebooks/inter_annual_variability/dustload_toa_imbalance_relationship.py.ipynb"

rule toa_load_relationships_NorESM:
    input:
        noresm_rad = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc',
                            vName=['rsdt','rlut','rlutcs','rsutcs'],
                            experiment=['piClim-control','piClim-aer', 'piClim-2xdust','piClim-2xss','piClim-lu','piClim-SO2', 'piClim-BC']),
        noresm_load = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc',
                        vName=['mmrdust','mmrss','mmrsoa','mmrbc','mmrso4','airmass',], experiment=['piClim-control','piClim-aer', 'piClim-2xdust', 
                                'piClim-2xss','piClim-lu','piClim-SO2', 'piClim-BC'] ),     
        noresm_clt = expand(outdir + '{experiment}/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc',
                        vName=['clt'], experiment=['piClim-control','piClim-aer', 'piClim-2xdust', 
                                'piClim-2xss','piClim-lu','piClim-SO2', 'piClim-BC'] ),     

        
    output:
        toa_load_sw_png = outdir + 'figs/antropogenic_forcing_variability/load_toa_imbalance_sw_NorESM.png',
        toa_load_png = outdir + 'figs/antropogenic_forcing_variability/load_toa_imbalance_NorESM.png'

    notebook:
        "../notebooks/inter_annual_variability/load_toa_imbalance_relationship_NorESM.py.ipynb"


rule remove_dust_variability:
    input:
        histSST_ERFcs = expand(outdir + 'histSST/ERFs/ERFtcs/ERFtcs_histSST_{model}_Ayear.nc',
            model=['EC-Earth3-AerChem', 'NorESM2-LM', 'MPI-ESM-1-2-HAM']),
        histSST_load = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                    vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM', 'NorESM2-LM']),   
        histSST_piaer_load = expand(outdir + 'histSST-piAer/{vName}/{vName}_histSST-piAer_{model}_Ayear.nc',
                    vName=['mmrdust','airmass'], model=['EC-Earth3-AerChem','MPI-ESM-1-2-HAM', 'NorESM2-LM']),
        toa_dustload_txt = outdir + 'figs/antropogenic_forcing_variability/dustload_toa_imbalance.csv'   
    output:
        outpath = outdir + 'figs/antropogenic_forcing_variability/clearsky_ERF_without_dustload.png',
        outpath1950 = outdir + 'figs/antropogenic_forcing_variability/clearsky_ERF_without_dustload_1950-1980.png',
    notebook:
         "../notebooks/inter_annual_variability/remove_dust_loadcontrib_from_ERFcs.py.ipynb"

rule noresm_dustload_forcing_variability:
    input:
        picli_vars = expand(outdir + 'piClim-control/{vName}/{vName}_piClim-control_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut', 'emidust', 'loaddust'], 
                    model=['NorESM2-LM']),
        piaer_vars = expand(outdir + 'piClim-aer/{vName}/{vName}_piClim-aer_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','emidust', 'loaddust'], 
                    model=['NorESM2-LM']),

        histSST = expand(outdir + 'histSST/{vName}/{vName}_histSST_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','emidust', 'loaddust'], 
                    model=['NorESM2-LM']),
        hist = expand(outdir + 'historical/{vName}/{vName}_historical_{model}_Ayear.nc',
                    vName=['rsut','rsdt','rlut','emidust', 'loaddust'], 
                    model=['NorESM2-LM'])
    
    output: 
        outpath = outdir + 'figs/antropogenic_forcing_variability/toa_variability_dust_NorESM2-LM.png'
    
    notebook:
        "../notebooks/inter_annual_variability/TOA_variability_and_dustload_noresm.py.ipynb" 

    
rule identify_forcing_regions:
    input:
        forcing_ec = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_EC-Earth3-AerChem_Ayear.nc',
        forcing_mpi = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_MPI-ESM-1-2-HAM_Ayear.nc',
        forcing_noresm = outdir + '{experiment}/ERFs/{vName}/{vName}_{experiment}_NorESM2-LM_Ayear.nc'
    output:
        mask = outdir + 'FORCeS_models_forcing_regions_{vName}_{experiment}.nc',
        plot = outdir + 'FORCeS_models_forcing_regions_plot_{vName}_{experiment}.png'
    notebook:
        "../examination_of_aerosol_forcing_regions.py.ipynb"

rule forcing_variability_aod:
    input:
        forcing_eff_noresm = outdir + 'forcing_efficiency/piClim-aer/ERFt_NorESM2-LM_delta_piClim-aer_piClim-control.nc',
        forcing_eff_mpi = outdir+'forcing_efficiency/piClim-aer/ERFt_MPI-ESM-1-2-HAM_delta_piClim-aer_piClim-control.nc',
        forcing_eff_ec = outdir+'forcing_efficiency/piClim-aer/ERFt_EC-Earth3-AerChem_delta_piClim-aer_piClim-control.nc',
        aod_noresm = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550dust','od550aer'],
                model='NorESM2-LM'),
        aod_mpi = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550aer'],
                model='MPI-ESM-1-2-HAM'),
        aod_ec = expand(rules.calc_change_historical.output.outpath,
                hist_pert='piClim-aer', hist_base='piClim-control',
                vName = ['od550aer', 'od550dust'],
                model='EC-Earth3-AerChem')
    output:
        outpath = outdir + 'figs/antropogenic_forcing_variability/forcing_variability_aod.png'
    notebook:
        "../notebooks/examine_variability_impact_on_forcing.py.ipynb"
