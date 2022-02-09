import pathlib as pl

def generate_experiment_lookup_table(root_dir, cmip_ver='CMIP6'):
    root = pl.Path(root_dir+'/'+cmip_ver)
    activity_paths = [elem for elem in root.iterdir() if elem.is_dir()]
    lookup_dir = {}
    for path in activity_paths:
        
        institu = [elem for elem in path.iterdir() if elem.is_dir()]
        activity = path.parts[-1]
        for inst in institu:
            for model in inst.iterdir():
                if model.is_dir():
                    for experiment in model.iterdir():
                        lookup_dir[experiment.parts[-1]]=activity
    return lookup_dir
    
def generate_institu_lookup_table(root_dir, cmip_ver='CMIP6'):
    root = pl.Path(root_dir+'/'+cmip_ver)
    activity_paths = [elem for elem in root.iterdir() if elem.is_dir()]
    lookup_dir = {}
    for path in activity_paths:
        
        institu = [elem for elem in path.iterdir() if elem.is_dir()]
        for inst in institu:
            for model in inst.iterdir():
                if model.is_dir():
                    lookup_dir[model.parts[-1]] = inst.parts[-1]
    return lookup_dir


def generate_file_ending_lookup_table(root_dir, cmip_ver='CMIP6'):
    root = pl.Path(root_dir+'/'+cmip_ver)
    activity_paths = [elem for elem in root.iterdir() if elem.is_dir()]
    lookup_dir = {}
    for path in activity_paths:
        if path.is_dir():
            activity = path.parts[-1]
            lookup_dir[activity] = {}
            institu = [elem for elem in path.iterdir() if elem.is_dir()]
            for inst in institu:
                for model in inst.iterdir():
                    if model.is_dir():
                        model_name = model.parts[-1] 
                        lookup_dir[activity][model_name] = {}
                        for experiment in model.iterdir():
                            if experiment.is_dir():
                                experiment_name = experiment.parts[-1]
                                lookup_dir[activity][model_name][experiment_name] = {}
                                varient = next(experiment.iterdir())
                                for table_id in varient.iterdir():
                                    if table_id.is_dir():
                                        variable = next(table_id.iterdir())

                                        fnames = []
                                        for gl in variable.iterdir():
                                            gl_iter = gl.iterdir()
                                            version = next(gl_iter)
                                            if version.is_symlink():
                                                version = next(gl_iter)
                                            if version.is_dir():    

                                                for fname in version.iterdir():
                                                    fnames.append(fname.parts[-1].split('_')[-1])
                                        lookup_dir[activity][model_name][experiment_name][table_id.parts[-1]] = fnames
                
    return lookup_dir