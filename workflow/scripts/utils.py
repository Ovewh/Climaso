import xarray as xr


def load_CMIP_data(path, **dataset_kwargs):

    if len(path) >1:
        print(path)
        return xr.open_mfdataset(path, chunks={'time':120}, **dataset_kwargs)
    else:
        dataset_kwargs.pop('data_vars')
        return xr.open_dataset(path[0], chunks={'time':120}, **dataset_kwargs)

