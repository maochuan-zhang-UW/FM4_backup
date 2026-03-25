#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 08:11:46 2025

@author: mczhang
"""

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Initialize the client
client = Client("IRIS", user="mczhang8@uw.edu", password="RxL42LH6XTNFybHb")

# Define parameters
network = "2F"
station = "*"
location = "*"
channel = "*"
starttime = UTCDateTime("2022-10-25T00:00:00")
endtime = UTCDateTime("2022-10-25T01:00:00")

# Download the data
st = client.get_waveforms(network=network, station=station, location=location,
                          channel=channel, starttime=starttime, endtime=endtime)

# Save to a MiniSEED file
st.write('/Users/mczhang/Documents/MATLAB/output_test.mseed', format='MSEED')

# Print summary
print(st)
