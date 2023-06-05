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
    Returns
    -------
        yearly_budget: annual emission budget

    """
    if ds[ds.variable_id].units == "kg m-2 s-1":
        yearly_budget = ds[ds.variable_id] * 365 * 24 * 60 * 60
    else:
        yearly_budget = ds[ds.variable_id]
    yearly_budget = yearly_budget * grid_area["cell_area"].load()
    yearly_budget = yearly_budget.sum(dim=["lon", "lat"]) * 1e-9
    yearly_budget.attrs["units"] = "Tg year-1"
    std_yearly_budget = yearly_budget.std(dim="time")
    mean_yearly_budget = yearly_budget.mean(dim="time")
    return mean_yearly_budget, std_yearly_budget




def regrid_global(
    ds: xr.Dataset,
    ds_out: xr.Dataset = None,
    lon: int = 3,
    lat: int = 2,
    method: str = "conservative",
    **regridder_kwargs,
):
    """
    Return regridded model output to get all the
    models on a common grid.

    params
    --------
        ds: Dataarray containting data field that should be regriddedd
    """
    if ds_out is None:
        ds_out = xe.util.grid_global(lon, lat, cf=True)
    regridder = xe.Regridder(ds, ds_out, method=method, **regridder_kwargs)
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
    ds_ctrl: xr.Dataset,
    ds_exp: xr.Dataset,
    time_slice: slice = None,
    time_average: bool = False,
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
    return _calc_change(
        ds_ctrl,
        ds_exp,
        rel_change=True,
        time_slice=time_slice,
        time_average=time_average,
    )


def calc_abs_change(
    ds_ctrl: xr.Dataset,
    ds_exp: xr.Dataset,
    time_slice=None,
    time_average: bool = False,
    vertical_sum: bool = True,
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
        vertical_sum=vertical_sum,
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

        if "lev" in ds_ctrl.dims and vertical_sum == True:
            da_ctrl = ds[f"control_{vName}"].sum(dim="lev")
            da_exp = ds[vName].sum(dim="lev")
        elif "lev" in ds_ctrl.dims and vertical_sum == False:
            da_ctrl = ds[f"control_{vName}"]
            da_exp = ds[vName]
        else:
            da_ctrl = ds[f"control_{vName}"]
            da_exp = ds[vName]

        if "year" in da_ctrl.dims:
            t_dim = "year"
        else:
            t_dim = "time"

        if time_slice:
            da_ctrl = da_ctrl.isel({f"{t_dim}": time_slice})
            da_exp = da_exp.isel({f"{t_dim}": time_slice})
        if time_average:
            da_ctrl = da_ctrl.mean({f"{t_dim}": time_slice})
            da_exp = da_exp.mean({f"{t_dim}": time_slice})

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


def model_levels_to_pressure_levels(ds: xr.Dataset | xr.DataArray):
    """
    Convert model levels to pressure levels.
    """
    import time
    import pathlib as pl

    try:
        import geocat.comp as geocomp
    except ImportError:
        raise ImportError(
            "geocat-comp is required for this function. Run snakemake with --use-conda to enable it."
        )
    ds = ds.copy()
    if "ap" in ds.data_vars:
        ds = ds.rename_vars({"ap": "a"})

    if "lev_bounds" in ds.data_vars:
        ds = ds.rename_vars({"lev_bounds": "lev_bnds"})

    if isinstance(ds, xr.Dataset):
        da = ds[ds.variable_id].copy()
    else:
        da = ds.copy()
    if ds.cf["Z"].formula == "p = a*p0 + b*ps":
        da = geocomp.interpolation.interp_hybrid_to_pressure(
            data=da,
            ps=ds["ps"],
            hyam=ds["a"],
            hybm=ds["b"],
            p0=ds.get("p0", 100000.0),
        )
    elif ds.cf["Z"].formula in ["p = ap + b*ps", "p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)"]:
        da = geocomp.interpolation.interp_hybrid_to_pressure(
            data=da, ps=ds["ps"], hyam=ds["a"], hybm=ds["b"], p0=ds.get("p0", 1)
        )

    else:
        print("formula not reconized")

    outds = da.to_dataset(name=ds.variable_id)
    outds.attrs = ds.attrs

    outds.attrs["history"] = (
        outds.attrs["history"]
        + f"@{time.ctime()} converted to standard pressure levels {pl.Path.cwd().parts[-1]}"
    )
    return outds


def pressure_to_height(pressure, p0=100000, t0=288.15):
    """
    Convert atmospheric pressure to altitude using the hypsometric equation.

    Parameters:
        pressure (float): Atmospheric pressure in Pascals.

    Returns:
        float: Altitude in meters.

    """
    # Define constants
    R = 287.058  # gas constant for dry air (J/(kg*K))
    g = 9.81  # acceleration due to gravity (m/s^2)
    L = 0.0065  # temperature lapse rate (°C/m)


    # Calculate temperature and height
    temperature = t0 - L * pressure / p0
    height = ((R / g * L) * np.log(p0 / pressure)
              + (t0 / L) * (1 - (pressure / p0)**(R * g / (L * R))))

    return height


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


def read_list_input_paths(path_list: list, models_pos: int = -2):
    """
    Read a list of input paths and return a dictionary
    of models and corresponding datasets and the variable name.

    Parameters
    ----------
        path_list : list

    Returns
    -------
        out_dict : dict
        vname : str

    """
    models = {pst.split("_")[models_pos].split(".")[0]: pst for pst in path_list}
    out_dict = {}
    for model, path in models.items():
        try:
            ds = xr.open_dataset(path)
        except ValueError:
            ds = xr.open_dataset(path, decode_times=False)
            ds = ds.drop(["time_bnds"])
            ds = xr.decode_cf(ds)

        out_dict[model] = ds.copy()
    vname = ds.variable_id
    return out_dict, vname


def calculate_pooled_variance(da_ctrl, da_exp):
    """
    Calculate the pooled variance of two samples.

    Parameters
    ----------
        da_ctrl : xr.DataArray or np.ndarray
        da_exp : xr.DataArray or np.ndarray

    Returns
    -------
        pooled_var : float
    """

    n1 = len(da_ctrl)
    n2 = len(da_exp)
    var1 = np.var(da_ctrl)
    var2 = np.var(da_exp)
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    return pooled_var


def calculate_CI(
    da_ctrl: xr.DataArray or np.ndarray,
    da_exp: xr.DataArray or np.ndarray,
    alpha: float = 0.05,
    global_mean: bool = False,
    mask: xr.DataArray = None,
    weights: xr.DataArray = None,
):
    """
    Calculate the confidence interval of the difference between two samples.

    Parameters
    ----------
        da_ctrl : xr.DataArray or np.ndarray
        da_exp : xr.DataArray or np.ndarray
        alpha : float, optional
            Significance level, by default 0.05
        global_mean : bool, optional
            If True, calculate the global mean of the two samples, by default False
        mask : xr.DataArray, optional
            Mask to apply to the two samples, by default None
        weights : xr.DataArray, optional
            Weights to apply to the two samples, by default None

    Returns
    -------
        CI : tuple
            Confidence interval
        diff : float
            Difference between the two samples
    """
    from scipy.stats import t

    if global_mean:
        ctrlM = global_avg(da_ctrl)
        expM = global_avg(da_exp)
        diff = np.mean(da_exp) - np.mean(da_ctrl)
    elif mask is not None or weights is not None:
        ctrlM = masked_average(da_ctrl, mask=mask, weights=weights, dim=["lat", "lon"])
        expM = masked_average(da_exp, mask=mask, weights=weights, dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    else:
        ctrlM = da_ctrl.mean(dim=["lat", "lon"])
        expM = da_exp.mean(dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    pooled_var = calculate_pooled_variance(ctrlM, expM)
    std_err = np.sqrt(pooled_var) * np.sqrt(1 / len(ctrlM) + 1 / len(expM))
    t_value = t.ppf(1 - alpha / 2, len(ctrlM) + len(expM) - 2)
    CI = (diff - t_value * std_err, diff + t_value * std_err)
    return CI[0], CI[1], diff


def t_test_diff_sample_means(
    da_ctrl: xr.DataArray,
    da_exp: xr.DataArray,
    global_mean: bool = False,
    mask: xr.DataArray = None,
    weights: xr.DataArray = None,
):
    """
    Perform a t-test for the difference between two samples.

    Parameters
    ----------
        da_ctrl : xr.DataArray or np.ndarray
        da_exp : xr.DataArray or np.ndarray
        global_mean : bool, optional
            If True, calculate the global mean, by default False
        mask : xr.DataArray, optional
            Mask to apply to the data, by default None
        weights : xr.DataArray, optional
            Weights to apply to the data, by default None

    Returns
    -------
        t_value : float
            t-value of the test
        p_value : float
            p-value of the test
        diff : float
            Difference between the two samples
    """
    from scipy.stats import ttest_ind

    if global_mean:
        ctrlM = global_avg(da_ctrl)
        expM = global_avg(da_exp)
        diff = np.mean(da_exp) - np.mean(da_ctrl)
    elif mask is not None or weights is not None:
        ctrlM = masked_average(da_ctrl, mask=mask, weights=weights, dim=["lat", "lon"])
        expM = masked_average(da_exp, mask=mask, weights=weights, dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    else:
        ctrlM = da_ctrl.mean(dim=["lat", "lon"])
        expM = da_exp.mean(dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    t_value, p_value = ttest_ind(expM, ctrlM, equal_var=True)
    return t_value, p_value, diff


def diff_means_greater_than_varability(
    da_ctrl: xr.DataArray,
    da_exp: xr.DataArray,
    global_mean: bool = False,
    mask: xr.DataArray = None,
    weights: xr.DataArray = None,
):
    """
    Check if the difference between two samples is greater than the variability of the samples.

    Parameters
    ----------
        da_ctrl : xr.DataArray or np.ndarray
        da_exp : xr.DataArray or np.ndarray
        global_mean : bool, optional
            If True, calculate the global mean of the samples, by default False
        mask : xr.DataArray, optional
            Mask to apply to the samples, by default None
        weights : xr.DataArray, optional
            Weights to apply to the samples, by default None

    Returns
    -------
        bool
    """
    if global_mean:
        ctrlM = global_avg(da_ctrl)
        expM = global_avg(da_exp)
        diff = np.mean(da_exp) - np.mean(da_ctrl)
    elif mask is not None or weights is not None:
        ctrlM = masked_average(da_ctrl, mask=mask, weights=weights, dim=["lat", "lon"])
        expM = masked_average(da_exp, mask=mask, weights=weights, dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    else:
        ctrlM = da_ctrl.mean(dim=["lat", "lon"])
        expM = da_exp.mean(dim=["lat", "lon"])
        diff = expM.mean() - ctrlM.mean()

    pooled_var = calculate_pooled_variance(da_ctrl, da_exp)
    return np.abs(diff) > np.sqrt(pooled_var)


def masked_average(
    da: xr.DataArray, mask: xr.DataArray = None, weights: xr.DataArray = None, dim=None
):
    """
    Calculate the average of a DataArray over a mask.

    Parameters
    ----------
        da : xr.DataArray
        mask : xr.DataArray
        weights : xr.DataArray
        dim : str, optional

    Returns
    -------
        avg : float
    """
    da_copy = da.copy()
    if mask is not None:
        if mask.dtype == "bool":
            da_copy = da_copy.where(mask)
        else:
            da_copy = da_copy.where(mask.isnull() == False)
        if weights is not None:
            weights = weights.where(mask)
    if weights is not None:
        avg = da_copy.weighted(weights).mean(dim=dim)
    else:
        avg = da_copy.mean(dim=dim)

    return avg

def make_equal_dimmension(da: xr.DataArray, refrence_da: xr.DataArray, dim: str = "time"):
    """
    Make a DataArray have the same dimensions as a refrence DataArray.

    Parameters
    ----------
        da : xr.DataArray
        refrence_da : xr.DataArray
        dim : str, optional

    Returns
    -------
        da : xr.DataArray
    """
    da = da.copy()
    if dim not in da.dims:
        da = da.expand_dims(dim=dim)
    da = da.reindex_like(refrence_da)
    return da