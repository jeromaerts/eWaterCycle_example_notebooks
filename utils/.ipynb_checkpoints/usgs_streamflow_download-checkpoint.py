# -*- coding: utf-8 -*-
"""
Author Jerom Aerts
j.p.m.aerts@tudelft.nl

This script downloads USGS streamflow data.
"""

import os
import sys
import urllib.request as urllib2

import numpy as np
import pandas as pd


def download_usgs_data(
    usgs_info_file,
    outputfolder,
    output_format,
    startDT,
    endDT,
    parameterCd,
    basin_start,
    basin_end,
    convert_unit_timestep=True
):
    # More information: https://waterservices.usgs.gov/rest/IV-Test-Tool.html
    # output_format e.g. ('json', 'rdb')
    # startDT e.g. ('1980-01-01')
    # endDT e.g.('2018-12-31')
    # parameterCd e.g. ('00060') discharge, cubic feet per second
    # 30208	m3/s
    # excel_format : [USGS_ID, lat, lon]!

    # Load USGS gauge ids
    stations = pd.read_table(usgs_info_file, delimiter=";")
    # Drop n number of rows, Remove after testing!
    stations = stations.iloc[basin_start:basin_end]
    # Create output folder
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    for index, station in stations.iterrows():
        # Fix import error -> adds 0 value leading gauge id
        usgsid = str(np.array(station["gauge_id"], dtype=np.int))
        if len(usgsid) == 7:
            usgsid = "0" + usgsid
        # Create download link
        url = (
            "https://waterservices.usgs.gov/nwis/iv/?format="
            + output_format
            + ",1.0&sites="
            + usgsid
            + "&startDT="
            + startDT
            + "&endDT="
            + endDT
            + "&parameterCd="
            + parameterCd
            + "&siteStatus=all"
        )
        out = outputfolder + "/" + usgsid + "." + output_format
        urllib2.urlretrieve(url, out)

        if convert_unit_timestep is True:

            # Changes format table to [datetime, discharge]
            df = pd.read_table(
                out,
                skiprows=55,
                usecols=[2, 4],
                header=None,
                names=["datetime", "discharge"],
            )  # Read rdb table and set column headers

            # Set Ice and Dis values to 0
            df["discharge"] = df["discharge"].replace("Ice", 0)
            df["discharge"] = df["discharge"].replace("Dis", 0)
            
            # Get equipment malfunction values
            indexvals = df.index[df["discharge"] == "Eqp"].tolist()
            for i in indexvals:
                # Overwrite with previous measurement
                df["discharge"].loc[i] = df["discharge"].loc[i - 1]

            # Convert to datetime and set index
            df.index = pd.to_datetime(
                df["datetime"], infer_datetime_format=True
            )
            # Drop obsolete column
            df = df.drop(columns="datetime")

            # Resample to hourly values
            df = df.resample("H").mean()

            # Print USGS ID for log
            print(" USGSID")
            print(usgsid)
            
            # Set table id and coordinates
            df["USGS_ID"] = usgsid
            df["lat"] = station["gauge_lat"]
            df["lon"] = station["gauge_lon"]
            
            # Convert series to float
            df["discharge"] = df["discharge"].apply(lambda x: float(x))
            # Convert to cubic meters per second
            df["discharge"] = df["discharge"].apply(lambda x: x / 35.315)

            # Write daily UTC files
            df.index = df.index.tz_localize("UTC")
            df = df.resample("D").mean()
            
            df["USGS_ID"] = usgsid
            df.to_csv(outputfolder + "/" + usgsid + "_UTC_daily.csv")
            log = print(usgsid + " Downloaded")
        else:
            exit()
    return log


# Set system variables from bash
usgs_info_file = sys.argv[1]
outputfolder = sys.argv[2]
start_date = sys.argv[3] 
end_date = sys.argv[4]
basin_start = int(sys.argv[5])
basin_end = int(sys.argv[6])

download_usgs_data(
    usgs_info_file, outputfolder, "rdb", start_date, end_date, "00060", basin_start, basin_end
)
