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

# Set min, max strahler order.
min_strahler_order = 4
max_strahler_order = 9

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
    
    # Obtain epsg number of coordinate system in UTM, for given latitude.
    # Check if northern or southern hemisphere
    if station_lat > 0:
        epsg_c = "EPSG:" + str(32601 + int(np.floor((station_lon - -180.0) / 6)))
    else:
        epsg_c = "EPSG:" + str(32701 + int(np.floor((station_lon - -180.0) / 6)))

    # Create lists for output dataframe
    df_info = pd.DataFrame()
    basin_names = []
    strahler_orders = []
    camels_areas = []
    wflow_areas = []
    area_ratios = []

    # Run over a number of stream order values to ensure that optimal basin area is obtained
    # Optimal area is defined wrt area as defined within camels dataset
    for strahler_order in range(min_strahler_order, max_strahler_order):
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

        # Obtain CAMELS basin area (km2)
        camels_area = row.area_gages2

        # Read basin geojson file as geopandas object
        basin_json_file = f"{OUTDIR}/staticgeoms/basins.geojson"
        gdf_json = gpd.read_file(basin_json_file)

        # Convert coordinates to metric system using epsg code and obtain catchment area (km2)
        wflow_area = (gdf_json.to_crs(epsg_c).area * 1e-6)[0]
        area_ratio = np.abs(np.log(camels_area / wflow_area))

        print(f"camels_area: {camels_area}")
        print(f"wflow_area: {wflow_area}")
        print(f"area ratio: {area_ratio}")

        basin_names.append(basin_name)
        strahler_orders.append(strahler_order)
        camels_areas.append(camels_area)
        wflow_areas.append(wflow_area)
        area_ratios.append(area_ratio)

    # Construct output dataframe
    df_info["basin_name"] = basin_names
    df_info["strahler_order"] = strahler_orders
    df_info["shape_area"] = camels_areas
    df_info["wflow_area"] = wflow_areas
    df_info["area_ratio"] = area_ratios

    # Select best strahler order value based on area ratio
    df_info = df_info.sort_values(
        ["strahler_order", "area_ratio"], ascending=(False, True)
    )
    strahler_order = df_info.loc[0].strahler_order
    print(f"Selected Strahler Order: {strahler_order}")

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