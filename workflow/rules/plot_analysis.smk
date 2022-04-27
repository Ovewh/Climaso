
rule plot_delta_ERF_control_experiment:
    input: 
        expand('results/ERFt_{experiment}_{models}_Ayear.nc', models=MODELS, allow_missing=True)
    output:
        outpath='figures/ERFt_{experiment}.png'
        
    