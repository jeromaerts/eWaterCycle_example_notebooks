import glob
import os
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapefile

# Define directory which will contain the wflow-sbm models for each basin in the camels dataset
PARAMDIR = Path("/gpfs/home6/jaerts/eWaterCycle_example_data/wflow_camels_parameters/")
CONFIGDIR = Path("/gpfs/home6/jaerts/eWaterCycle_example_notebooks/hydromt/")

# Set hydromt sources .yml and .ini file
datasource = f"{CONFIGDIR}/data_sources_hydromt.yml"
hydromt_config = f"{CONFIGDIR}/hydromt_wflow_sbm.ini"

# Perform wflow at 3 different resolutions, corresponding to ~ 200m, 1000m, and 3000m pixel resolution
resolution = [0.00166667, 0.00833333, 0.025]
resolution_names = ["200m", "1km", "3km"]

# Set min, max strahler order.
min_strahler_order = 4
max_strahler_order = 9

# Obtain information about each basin within the camels dataset
camels_shapefile = "/gpfs/home6/jaerts/eWaterCycle_example_data/camels_input_files/camels_basin_shapefiles/HCDN_nhru_final_671.shp"
camels_topo_file = (
    "/gpfs/home6/jaerts/eWaterCycle_example_data/camels_input_files/camels_topo.txt"
)
df_camels_topo = pd.read_table(camels_topo_file, delimiter=";")

# Fix missing leading zero Gauge ID
gauge_ids = df_camels_topo.gauge_id.to_list()
gauge_ids = [str(x) for x in gauge_ids]
gauge_ids = ["0" + x if len(x) == 7 else x for x in gauge_ids]
df_camels_topo["gauge_id"] = gauge_ids

# Read CAMELS shapefile.
gdf_shape = gpd.read_file(camels_shapefile)

# Fix missing leading zero Gauge ID
gdf_shape["hru_id"] = gauge_ids
gdf_shape = gdf_shape.set_index("hru_id")


# Start HydroMT loop
for resolution_name in resolution_names:

    # Check if resolution information exists.
    RESDIR = Path(f"{PARAMDIR}/{resolution_name}/")
    
    # Create directory structure.
    if not os.path.exists(RESDIR):
        os.makedirs(RESDIR)
    
    # Iterate gauge IDs
    for index, row in df_camels_topo.iterrows():

        # Set output directory
        basin_name = f"{index}_camels_{row.gauge_id}_{resolution_name}"
        OUTDIR = Path(f"{RESDIR}/{basin_name}/")

        # Set streamflow observation station coordinates
        station_lat = row.gauge_lat
        station_lon = row.gauge_lon

        # Get basin geometry
        bbox = gdf_shape.loc[row.gauge_id].geometry.bounds
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