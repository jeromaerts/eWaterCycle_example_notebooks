"""
Authors: Jerom Aerts, Pieter Hazenberg
Contact: j.p.m.aerts@tudelft.nl

This script runs HydroMT for multiple basins defined by shapefiles. The scripts iterates
multiple strahler order numbers. The best strahler order is selected based on the match
between derived parameter set catchment area and the shapefile defined catchment area.
The resulting parameter set is derived using the optimal strahler order value to catchment
area ratio.
"""

import glob
import os
import shutil
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapefile

# Define the shapefile input directory and get shapefiles.
SHAPEDIR = Path("/gpfs/home6/jaerts/hydromt/catchment_shapefiles/")
shapefiles = glob.glob(f"{SHAPEDIR}/*.shp")

# Define parameter set output directory.
PARAMDIR = Path("/gpfs/home6/jaerts/hydromt/hydromt_wflow_parameters/")

# Set hydromt sources .yml and .ini file
datasource = f"{PARAMDIR}/data_sources_deltares.yml"
hydromt_config = f"{PARAMDIR}/config_hydromt_wflow_sbm.ini"

# Derive wflow_sbm at 1km pixel resolution.
resolution = 0.00833333
resolution_name = "1km"

# Check if resolution information exists.
RESDIR = Path(f"{PARAMDIR}/{resolution_name}/")

# Create directory structure.
if not os.path.exists(RESDIR):
    os.makedirs(RESDIR)

# Run hydromt for multiple basins defined as shapefiles.
for shapefile in shapefiles:

    # Set output directory
    basin_name = shapefile.split("/")[-1].split(".")[0]

    OUTDIR = Path(f"{RESDIR}/wflow_sbm_{resolution_name}_{basin_name}/")

    # Hardcoded lat lon station location, this should be retrieved from a table or shapefile
    station_lat = 47.596523
    station_lon = 8.6819796

    # Read shapefile.
    gdf_shape = gpd.read_file(shapefile)

    # Get basin geometry
    bbox = gdf_shape.loc[0].geometry.bounds
    left_bound = str(np.round(1e6 * (bbox[0] - 0.25)) / 1e6)
    bottom_bound = str(np.round(1e6 * (bbox[1] - 0.25)) / 1e6)
    right_bound = str(np.round(1e6 * (bbox[2] + 0.25)) / 1e6)
    top_bound = str(np.round(1e6 * (bbox[3] + 0.25)) / 1e6)

    # Set minimum strahler order.
    min_strahler_order = 6
    min_ratio = 10
    area_ratio = min_ratio + 1

    # Create lists for output dataframe
    df_info = pd.DataFrame()
    basin_names = []
    strahler_orders = []
    shape_areas = []
    wflow_areas = []
    area_ratios = []

    # Run over a number of stream order values to ensure that optimal catchment area is obtained
    # Optimal area is defined wrt area as defined within camels dataset
    while area_ratio > min_ratio:
        print("yes")
        for strahler_order in range(min_strahler_order, 9):
            print(f"strahler_order: {strahler_order}")

            if os.path.exists(OUTDIR):
                shutil.rmtree(OUTDIR)

            runbasin = (
                "hydromt build wflow -vv "
                + str(OUTDIR)
                + " -r "
                + str(resolution)
                + " \"{'subbasin': ["
                + str(station_lon)
                + ", "
                + str(station_lat)
                + "], 'strord': "
                + str(strahler_order)
                + ", 'bbox': ["
                + left_bound
                + ", "
                + bottom_bound
                + ", "
                + right_bound
                + ", "
                + top_bound
                + ']}" -i '
                + str(hydromt_config)
                + " -d "
                + datasource
            )

            print(runbasin)
            os.system(runbasin)

            catchment_json_file = f"{OUTDIR}/staticgeoms/basins.geojson"

            # Read catchment geojson file as geopandas object.
            gdf_json = gpd.read_file(catchment_json_file)

            # Obtain catchment area (km2)
            shape_area = (
                gdf_shape["geometry"]
                .to_crs({"init": "epsg:3395"})
                .map(lambda p: p.area / 10 ** 6)
                .values
            )
            wflow_area = (
                gdf_json["geometry"]
                .to_crs({"init": "epsg:3395"})
                .map(lambda p: p.area / 10 ** 6)
                .values
            )
            area_ratio = shape_area / wflow_area

            print(f"shape_area: {shape_area}")
            print(f"wflow_area: {wflow_area}")
            print(f"area ratio: {area_ratio}")

            basin_names.append(basin_name)
            strahler_orders.append(strahler_order)
            shape_areas.append(shape_area)
            wflow_areas.append(wflow_area)
            area_ratios.append(area_ratio)

    df_info["basin_name"] = basin_names
    df_info["strahler_order"] = strahler_orders
    df_info["shape_area"] = shape_areas
    df_info["wflow_area"] = wflow_areas
    df_info["area_ratio"] = area_ratios

    # Select best strahler order value based on area ratio
    df_info = df_info.set_index("area_ratio")
    df_info = df_info.sort_index()
    strahler_order = df_info.strahler_order[0]

# Final run HydroMT with best strahler order
runbasin = (
    "hydromt build wflow -vv "
    + str(OUTDIR)
    + " -r "
    + str(resolution)
    + " \"{'subbasin': ["
    + str(station_lon)
    + ", "
    + str(station_lat)
    + "], 'strord': "
    + str(strahler_order)
    + ", 'bbox': ["
    + left_bound
    + ", "
    + bottom_bound
    + ", "
    + right_bound
    + ", "
    + top_bound
    + ']}" -i '
    + str(hydromt_config)
    + " -d "
    + datasource
)

print(runbasin)
os.system(runbasin)