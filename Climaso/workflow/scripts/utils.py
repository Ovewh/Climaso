import xarray as xr

def load_CMIP_data(path, **dataset_kwargs):
    if len(path) >1:
        return xr.open_mfdataset(path, chunks={'time':120}, **dataset_kwargs)
    else:
        return xr.open_dataset(path, chunks={'time':120}, **dataset_kwargs)