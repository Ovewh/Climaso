import xarray as xr
import xesmf as xe
from pyclim_noresm.aerosol_forcing import merge_exp_ctrl
from pyclim_noresm.general_util_funcs import global_avg
import numpy as np




def make_consistent(dsets):
    template_ds = dsets[0]
    fixed_dsets = []
    for ds in dsets:
        ds = ds.reindex(
            {"lon": template_ds.lon, "lat": template_ds.lat}, method="nearest"
        )
        ds = ds.assign(
            {"lon_bnds": template_ds.lon_bnds, "lat_bnds": template_ds.lat_bnds}
        )
        fixed_dsets.append(ds)
    return fixed_dsets


def load_CMIP_data(path, **dataset_kwargs):

    if len(path) > 1:
        try:
            ds = xr.open_mfdataset(path, chunks={"time": 120}, **dataset_kwargs)
        except xr.MergeError:
            dvar = dataset_kwargs.pop("data_vars")
            if not isinstance(dvar, list):
                dvar = list(dvar)
            dsets = [xr.open_dataset(p, **dataset_kwargs) for p in path]
            dsets = make_consistent(dsets)
            ds = xr.concat(dsets, dim="time", data_vars=dvar)
        return ds
    else:
        dataset_kwargs.pop("data_vars")
        return xr.open_dataset(path[0], chunks={"time": 120}, **dataset_kwargs)


def copy_meta_data_CMIP(attrs):
    remove_keys = [
        "source",
        "table_info",
        "variable_id",
        "tracking_id",
        "history",
        "creation_date",
        "branch_time_in_child",
        "branch_time_in_parent",
        "Conventions",
    ]
    new_attrs = {k: attrs[k] for k in set(list(attrs.keys())) - set(remove_keys)}
    return new_attrs


def compute_annual_emission_budget(ds: xr.Dataset, grid_area: xr.Dataset):
    """
    Computes annual emission budget from monthly data

    Params
    ------
        ds :  monthly emission flux,
        grid_area: area of each grid cell
    """
    yearly_budget = ds[ds.variable_id].resample(time="Y").mean() * 365 * 24 * 60 * 60

    yearly_budget = yearly_budget * grid_area["cell_area"]
    yearly_budget = yearly_budget.sum(dim=["lon", "lat"]) * 1e-9
    yearly_budget.attrs["units"] = "Tg year-1"
    mean_yearly_budget = yearly_budget.mean(dim="time")
    return mean_yearly_budget


def regrid_global(ds: xr.Dataset, ds_out: xr.Dataset=None, 
                        lon: int = 3, 
                        lat: int = 2,
                        method: str = "conservative",):
    """
    Return regridded model output to get all the
    models on a common grid.

    params
    --------
        ds: Dataarray containting data field that should be regriddedd
    """
    if not ds_out:
        ds_out = xe.util.grid_global(lon, lat, cf=True)
    regridder = xe.Regridder(ds, ds_out, method=method)
    ds = regridder(ds, keep_attrs=True)
    return ds


def calc_error_gridded(da: xr.DataArray, time_dim: str = "year", kind="std"):
    if isinstance(da, xr.Dataset):
        da = da[da.variable_id]
    with xr.set_options(keep_attrs=True):
        std = da.std(dim=time_dim)
        if kind == "sem" or kind == "SEM":
            st_error = std / np.sqrt(len(da[time_dim]))
        else:
            st_error = std
        try:
            st_error.attrs["long_name"] = "{} of {}".format(
                kind, st_error.attrs["long_name"]
            )
        except KeyError:
            st_error.attrs["long name"] = "{} of {}".format(
                kind, st_error.attrs["long name"]
            )

    return st_error


def calc_error(da: xr.DataArray, time_dim: str = "year", kind="std"):
    """
    Calculates the standard error of the mean.

    params:
        da: DataArray|Dataset containing timeseries data:

    """
    if isinstance(da, xr.Dataset):
        da = da[da.variable_id]
    with xr.set_options(keep_attrs=True):
        if "lat" in da.dims:
            da = global_avg(da)

        std = da.std(dim=time_dim)
        if kind == "sem" or kind == "SEM":
            st_error = std / np.sqrt(len(da[time_dim]))
        else:
            st_error = std
        try:
            st_error.attrs["long_name"] = "{} of {}".format(
                kind, st_error.attrs["long_name"]
            )
        except KeyError:
            st_error.attrs["long name"] = "{} of {}".format(
                kind, st_error.attrs["long name"]
            )

    return st_error


def calc_relative_change(
    ds_ctrl: xr.Dataset, ds_exp: xr.Dataset, time_slice: slice = None, time_average: bool=False
):
    """
    Calculate the relative change of a diagnostic between
    the control and experiment. The data is on model level
    it is integrated vertically.

    params
    ------
        ds_ctrl: dataset containing the control experiment
        ds_exp: dataset containing the perturbed experiment
        time_slice: slice time based on index

    returns
    --------
        rel_diff: dataset containing the relative difference
    """
    return _calc_change(ds_ctrl, ds_exp, rel_change=True, time_slice=time_slice,time_average=time_average)


def calc_abs_change(
    ds_ctrl: xr.Dataset, ds_exp: xr.Dataset, time_slice=None, time_average: bool=False, vertical_sum: bool=True
):
    """
    Calculated absolute change of diagnostic between
    the control and experiment. Data that is on model
    level is integrated vertically.

    params:
    -------
        ds_ctrl: dataset containing the control experiment
        ds_exp: dataset containing the perturbed experiment
        time_slice : slice time based on index

    """
    return _calc_change(
        ds_ctrl,
        ds_exp,
        rel_change=False,
        time_slice=time_slice,
        time_average=time_average,
        vertical_sum =vertical_sum,
    )


def _calc_change(
    ds_ctrl: xr.Dataset,
    ds_exp: xr.Dataset,
    rel_change: bool = False,
    time_slice=None,
    time_average=False,
    vertical_sum=True,
):

    vName = ds_ctrl.variable_id
    dvars = set(list(ds_ctrl.data_vars))
    keep_vars = set(["lev_bounds", "time_bounds", "lon_bnds", "lat_bnds", vName])

    with xr.set_options(keep_attrs=True):
        ds_ctrl = ds_ctrl.drop_vars(dvars - keep_vars).squeeze()
        ds_exp = ds_exp.drop_vars(dvars - keep_vars).squeeze()
        ds = merge_exp_ctrl(ds_ctrl, ds_exp)

        if "lev" in ds_ctrl.dims and vertical_sum==True:
            da_ctrl = ds[f"control_{vName}"].sum(dim="lev")
            da_exp = ds[vName].sum(dim="lev")
        elif "lev" in ds_ctrl.dims and vertical_sum==False:
            da_ctrl = ds[f"control_{vName}"]
            da_exp = ds[vName]
        else:
            da_ctrl = ds[f"control_{vName}"]
            da_exp = ds[vName]
        
        if "year" in da_ctrl.dims:
            t_dim='year'
        else:
            t_dim='time'
        
        if time_slice:
            da_ctrl = da_ctrl.isel({f'{t_dim}':time_slice})
            da_exp = da_exp.isel({f'{t_dim}':time_slice})
        if time_average:
            da_ctrl = da_ctrl.mean({f'{t_dim}':time_slice})
            da_exp = da_exp.mean({f'{t_dim}':time_slice})

        if rel_change:
            diff = ((da_exp - da_ctrl) / da_ctrl) * 100
            diff.attrs["units"] = "%"
            diff.attrs["long_name"] = "Relative change of {}".format(
                da_exp.attrs["long_name"]
            )
            diff = diff.to_dataset(name=vName)
            diff.attrs = {**diff.attrs, **ds_exp.attrs}
        else:
            diff = da_exp - da_ctrl
            diff.attrs[
                "long_name"
            ] = "Difference between experiment and control of {}".format(
                da_exp.attrs["long_name"]
            )
            diff = diff.to_dataset(name=vName)
            diff.attrs = {**diff.attrs, **ds_exp.attrs}
    return diff


def transelate_aerocom_helper(wildcards):
    if wildcards.freq == "2010":
        freq = "clim"
    elif wildcards.freq == "monthly":
        freq = "Amon"
    elif wildcards.freq == "yearly":
        freq = "Ayear"
    else:
        freq = wildcards.freq

    return freq
