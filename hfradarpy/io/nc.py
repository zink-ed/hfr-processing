import datetime as dt
import numpy as np
import pandas as pd
import pprint
import xarray as xr
from collections import OrderedDict
import os

datetime_format = "%Y%m%dT%H%M%SZ"


def make_encoding(ds, time_start="seconds since 1970-01-01 00:00:00", comp_level=4, chunksize=10000, fillvalue=-999.00):
    encoding = {}

    for k in ds.data_vars:
        values = ds[k].values
        shape = values.shape

        encoding[k] = {"zlib": True, "complevel": comp_level, "_FillValue": np.float32(fillvalue)}

        if 0 not in shape:
            if values.dtype.kind == "O":
                values = values.astype("str")

            if values.dtype.kind == "S":
                size = values.dtype.itemsize
                if size > 1:
                    shape = shape + (size,)

            dim0 = min(shape[0], chunksize)
            shape = (dim0,) + shape[1:]
            encoding[k]["chunksizes"] = shape

    # add the encoding for time so xarray exports the proper time
    encoding["time"] = dict(units=time_start, calendar="gregorian", zlib=False, _FillValue=None, dtype=np.double)
    # encoding['site_code_flags'] = dict(zlib=True, _FillValue=int(0))

    return encoding


def required_global_attributes(required_attributes, time_start, time_end):
    """
    Call this function to make your netcdf have CF compliant attributes

    Args:
        required_attributes (dict): any additional attributes you want to add.
        time_start (_type_): start time of dataset
        time_end (_type_): end time of dataset

    Returns:
        _type_: _description_
    """
    datetime_format = "%Y%m%dT%H%M%SZ"
    created = pd.Timestamp(dt.datetime.utcnow()).strftime(datetime_format)  # creation time Timestamp

    time_start = pd.Timestamp(str(time_start)).strftime(datetime_format)
    time_end = pd.Timestamp(str(time_end)).strftime(datetime_format)

    # Required global attributes
    global_attrs = [
        ("publisher_name", required_attributes["publisher_name"]),
        ("publisher_email", required_attributes["publisher_email"]),
        ("publisher_url", required_attributes["publisher_url"]),
        ("title", required_attributes["title"]),
        ("summary", required_attributes["summary"]),
        ("keywords", required_attributes["keywords"]),
        ("Conventions", "CF-1.6, ACDD-1.3"),
        ("naming_authority", required_attributes["naming_authority"]),
        ("history", required_attributes["history"]),
        ("source", required_attributes["source"]),
        ("processing_level", required_attributes["processing_level"]),
        ("comment", required_attributes["comment"]),
        ("acknowledgment", required_attributes["acknowledgment"]),
        ("standard_name_vocabulary", "CF Standard Name Table v41"),
        ("date_created", created),
        ("creator_name", required_attributes["creator_name"]),
        ("creator_email", required_attributes["creator_email"]),
        ("creator_url", required_attributes["creator_url"]),
        ("institution", required_attributes["institution"]),
        ("project", required_attributes["project"]),
        ("geospatial_lat_min", -90),
        ("geospatial_lat_max", 90),
        ("geospatial_lon_min", -180),
        ("geospatial_lon_max", 180),
        ("geospatial_vertical_min", 0.0),
        ("geospatial_vertical_max", 0.0),
        ("geospatial_vertical_positive", "down"),
        ("time_coverage_start", time_start),
        ("time_coverage_end", time_end),
        ("sea_name", required_attributes["sea_name"]),
        ("creator_type", "person"),
        ("creator_institution", required_attributes["creator_institution"]),
        ("contributor_name", required_attributes["contributor_name"]),
        ("contributor_role", required_attributes["contributor_role"]),
        ("geospatial_lat_units", "degrees_north"),
        ("geospatial_lon_units", "degrees_east"),
        ("date_modified", created),
        ("date_issued", created),
        ("date_metadata_modified", created),
        ("keywords_vocabulary", "GCMD Science Keywords"),
        ("platform", required_attributes["platform"]),
        ("instrument", required_attributes["instrument"]),
        ("cdm_data_type", required_attributes["cdm_data_type"]),
        ("references", required_attributes["references"]),
    ]

    if "ncei_template_version" in required_attributes.keys():
        global_attrs = [("ncei_template_version", required_attributes["ncei_template_version"])] + global_attrs

    global_attrs = OrderedDict(global_attrs)
    return global_attrs


def aggregate(files, save_dir, save_filename=None):
    """
    This function allows you to aggregate multiple netcdf files into a single file.
    It will concatenate on the coordinates of the netcdf files

    Args:
        files (list): list of files or regular expression to location of files you want to open
        save_dir (str): directory to save aggregated netcdf file
        save_filename (str, optional): filename of aggregated netcdf file. Defaults to None.

    Returns:
        str: filename of aggregated netcdf file.
    """
    # Create save directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)

    print("Aggregating the following datasets:")
    pprint.pprint(files)

    # Opening files lazily (not into memory) using xarray
    ds = xr.open_mfdataset(files, combine="by_coords")
    ds.attrs["time_coverage_start"] = pd.Timestamp(str(ds["time"].min().data)).strftime(datetime_format)
    ds.attrs["time_coverage_end"] = pd.Timestamp(str(ds["time"].max().data)).strftime(datetime_format)

    # Encode variables for efficiency reasons
    encoding = make_encoding(ds)
    for v in list(ds.coords):
        encoding[v] = dict(zlib=False, _FillValue=False)

    print("Saving aggregated datasets as netCDF4 file.")
    ds.load()  # Load lazy arrays of open files into memory. Performance is better once loaded

    if save_filename:
        save_file = "{}/{}-{}_{}".format(save_dir, ds.time_coverage_start, ds.time_coverage_end, save_filename)
        ds.to_netcdf(save_file, encoding=encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])
    else:
        save_file = "{}/{}-{}_totals_aggregated.nc".format(save_dir, ds.time_coverage_start, ds.time_coverage_end)
        ds.to_netcdf(save_file, encoding=encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])
    ds.close()
    return save_file
