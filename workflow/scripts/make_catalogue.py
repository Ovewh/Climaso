import ecgtools
from ecgtools import Builder
from ecgtools.parsers.cmip import parse_cmip6_using_directories
lustre = snakemake.config.get('seperate_Noresm', True)
if snakemake.wildcards.source == 'noresm' or snakemake.wildcards.source == 'noresmnirdtoolkit':
    exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*']
elif lustre == False:
    exclude_patterns=['*/files/*']

else:
    exclude_patterns=['*/files/*', '*/latest','.cmorout/*', '*/NorCPM1/*', '*/NorESM1-F/*','*/NorESM2-LM/*','*/NorESM2-MM/*']    
builder = Builder(paths=[snakemake.params.root_path], depth=snakemake.params.depth,
                joblib_parallel_kwargs={'n_jobs': snakemake.threads, 'verbose':13},
                exclude_patterns=exclude_patterns)
builder.build(parsing_func=parse_cmip6_using_directories)

builder.clean_dataframe()
df = builder.df
df.to_csv(snakemake.output.outpath, compression='gzip', index=False)
