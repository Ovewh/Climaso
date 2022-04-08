from attr import attr
import xarray as xr


def load_CMIP_data(path, **dataset_kwargs):

    if len(path) >1:
        return xr.open_mfdataset(path, chunks={'time':120}, **dataset_kwargs)
    else:
        print(dataset_kwargs)
        dataset_kwargs.pop('data_vars')
        return xr.open_dataset(path[0], chunks={'time':120}, **dataset_kwargs)


def copy_meta_data_CMIP(attrs):
    remove_keys = ['source','table_info', 'variable_id','tracking_id',
                    'history','creation_date', 'branch_time_in_child',
                    'branch_time_in_parent','Conventions']
    new_attrs = {k: attrs[k] for k in set(list(attrs.keys())) - set(remove_keys)}
    return new_attrs

def transelate_aerocom_helper(wildcards):
    if wildcards.freq=='9999':
        freq='clim'
    elif wildcards.freq=='monthly':
        freq='Amon'
    elif wildcards.freq=='yearly':
        freq='Ayear'
    else:
        freq=wildcards.freq
        
    return freq