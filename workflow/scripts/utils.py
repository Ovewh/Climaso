import xarray as xr
import xesmf as xe
from pyclim_noresm.aerosol_forcing import merge_exp_ctrl

def load_CMIP_data(path, **dataset_kwargs):

    if len(path) >1:
        return xr.open_mfdataset(path, chunks={'time':120}, **dataset_kwargs)
    else:
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


def calc_relative_change(ds_ctrl: xr.Dataset, ds_exp: xr.Dataset):
    """
    Calculate the relative change of a diagnostic between
    the control and experiment. The data is on model level 
    it is integrated vertically.

    params
    ------
        ds_ctrl: dataset containing the control experiment
        ds_exp: dataset containing the perturbed experiment
    
    returns
    --------
        rel_diff: dataset containing the relative difference
    """
    return _calc_change(ds_ctrl, ds_exp,rel_change=True) 

def calc_abs_change(ds_ctrl:xr.Dataset, ds_exp: xr.Dataset):
    """
    Calculated absolute change of diagnostic between
    the control and experiment. Data that is on model
    level is integrated vertically.

    params:
    -------
        ds_ctrl: dataset containing the control experiment
        ds_exp: dataset containing the perturbed experiment
    
    """
    return _calc_change(ds_ctrl, ds_exp, rel_change=False)

def _calc_change(ds_ctrl:xr.Dataset, ds_exp:xr.Dataset, rel_change: bool =False):

    vName = ds_ctrl.variable_id
    dvars = set(list(ds_ctrl.data_vars))
    keep_vars = set(['lev_bounds','time_bounds', 'lon_bnds','lat_bnds',vName])
    
    with xr.set_options(keep_attrs=True):
        ds_ctrl = ds_ctrl.drop_vars(dvars-keep_vars).squeeze()
        ds_exp = ds_exp.drop_vars(dvars-keep_vars).squeeze()
        ds = merge_exp_ctrl(ds_ctrl, ds_exp)

        if 'lev' in ds_ctrl.dims:
            da_ctrl = ds[f'control_{vName}'].sum(dim='lev')
            da_exp = ds[vName].sum(dim='lev')
        else:
            da_ctrl = ds[f'control_{vName}']
            da_exp = ds[vName]
        if 'year' in da_ctrl.dims:
            da_ctrl = da_ctrl.mean(dim='year')
            da_exp = da_exp.mean(dim='year')
        if rel_change:
            diff = ((da_exp-da_ctrl)/da_ctrl)*100
            diff.attrs['units'] = '%'
            diff.attrs['long_name'] = 'Relative change of {}'.format(da_exp.attrs['long_name']) 
            diff = diff.to_dataset(name=vName)
            diff.attrs = {**diff.attrs, **ds_exp.attrs}
        else:
            diff = da_exp-da_ctrl
            diff.attrs['long_name'] = 'Difference between experiment and control of {}'.format(da_exp.attrs['long_name'])
            diff = diff.to_dataset(name=vName)
            diff.attrs = {**diff.attrs, **ds_exp.attrs}
    return diff

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