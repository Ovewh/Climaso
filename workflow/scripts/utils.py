from attr import attr
import xarray as xr
import xesmf as xe

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

def regrid_global(ds: xr.DataArray, base_ds: xr.Dataset, lon: int =3, lat: int=2):
    """
    Return regridded model output to get all the 
    models on a common grid.

    params
    --------
        ds: Dataarray containting data field that should be regriddedd 
    """
    ds_out = xe.util.grid_global(lon,lat, cf=True)
    regridder = xe.Regridder(ds, ds_out, 'conservative')
    ds = regridder(ds, keep_attrs=True)
    return ds


def transelate_aerocom_helper(wildcards):
    print(wildcards)
    if wildcards.freq=='2010':
        freq='clim'
    elif wildcards.freq=='monthly':
        freq='Amon'
    elif wildcards.freq=='yearly':
        freq='Ayear'
    else:
        freq=wildcards.freq
        
    return freq