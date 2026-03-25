import os
import logging
import numpy as np
import scipy.io as sio
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def download_seismic_data(
    start_date="2022-09-05",
    end_date="2022-09-05T09:00:00",
    network="2F",
    station_pattern="AX*",
    location_pattern="*",
    channel_patterns=["HH?", "EL?"],
    output_folder="data",
    username="mczhang8@uw.edu",
    token="RxL42LH6XTNFybHb",
    increment_hours=1,
    overlap_seconds=15
):
    """
    Download seismic data from IRIS hourly with overlap and save as mat files.

    Args:
        start_date (str): Start date in YYYY-MM-DD format
        end_date (str): End date in YYYY-MM-DD format
        network (str): Network code
        station_pattern (str): Station pattern (e.g., "AX*")
        location_pattern (str): Location pattern
        channel_patterns (list): List of channel patterns
        output_folder (str): Base output directory (data/yyyy/mm/)
        username (str): IRIS username
        token (str): IRIS authentication token
        increment_hours (float): Time increment in hours
        overlap_seconds (float): Overlap in seconds between intervals

    Returns:
        dict: Trace metadata and data for the last interval
    """
    # Initialize IRIS client
    try:
        client = Client("IRIS", user=username, password=token)
        logger.info("Initialized IRIS client")
    except Exception as e:
        logger.error(f"Failed to initialize IRIS client: {str(e)}")
        return {}

    # Convert dates
    try:
        start = UTCDateTime(start_date)
        end = UTCDateTime(end_date)
        logger.info(f"Processing date range: {start_date} to {end_date}")
    except Exception as e:
        logger.error(f"Invalid date format: {str(e)}")
        return {}

    # Get station inventory
    try:
        inventory = client.get_stations(
            network=network,
            station=station_pattern,
            location=location_pattern,
            channel=",".join(channel_patterns),
            starttime=start,
            endtime=end,
            level="channel"
        )
        logger.info(f"Found {len(inventory.networks)} networks, {sum(len(sta) for net in inventory.networks for sta in net.stations)} stations")
        for net in inventory.networks:
            for sta in net.stations:
                logger.info(f"Station {sta.code}: {len(sta.channels)} channels")
    except Exception as e:
        logger.error(f"Failed to retrieve inventory: {str(e)}")
        return {}

    # Convert increment and overlap to seconds
    increment_seconds = increment_hours * 3600
    overlap = overlap_seconds / 86400  # Convert seconds to days for UTCDateTime

    # Initialize last trace dictionary
    last_trace_dict = {}

    # Loop through each hour
    current_date = start
    while current_date < end:
        interval_start = current_date - overlap
        interval_end = current_date + increment_seconds + overlap
        date_str = current_date.strftime('%Y-%m-%d-%H-%M-%S')
        logger.info(f"Processing interval: {interval_start.strftime('%Y-%m-%d %H:%M:%S.%f')} to {interval_end.strftime('%Y-%m-%d %H:%M:%S.%f')}")

        # Initialize data containers for this interval
        trace_dict = {
            'network': [], 'station': [], 'location': [], 'channel': [],
            'sensitivity': [], 'sensitivityFrequency': [], 'data': [],
            'sampleCount': [], 'sampleRate': [], 'startTime': [], 'endTime': []
        }

        # Create output directory (data/yyyy/mm/)
        year = current_date.strftime('%Y')
        month = current_date.strftime('%m')
        output_dir = os.path.join(output_folder, year, month)
        os.makedirs(output_dir, exist_ok=True)

        for net in inventory.networks:
            for sta in net.stations:
                for chan in sta.channels:
                    try:
                        logger.info(f"Downloading {net.code}.{sta.code}.{chan.location_code}.{chan.code}")
                        stream = client.get_waveforms(
                            network=net.code,
                            station=sta.code,
                            location=chan.location_code,
                            channel=chan.code,
                            starttime=interval_start,
                            endtime=interval_end,
                            attach_response=True
                        )

                        if not stream:
                            logger.warning(f"No data for {net.code}.{sta.code}.{chan.location_code}.{chan.code}")
                            continue

                        # Process each trace
                        for trace in stream:
                            try:
                                resp = trace.stats.response._get_overall_sensitivity_and_gain()
                                trace_dict['network'].append(trace.stats.network)
                                trace_dict['station'].append(trace.stats.station)
                                trace_dict['location'].append(trace.stats.location)
                                trace_dict['channel'].append(trace.stats.channel)
                                trace_dict['sensitivity'].append(float(resp[1]))
                                trace_dict['sensitivityFrequency'].append(float(resp[0]))
                                trace_dict['data'].append(trace.data.astype(np.float64))
                                trace_dict['sampleCount'].append(int(trace.stats.npts))
                                trace_dict['sampleRate'].append(float(trace.stats.sampling_rate))
                                trace_dict['startTime'].append(trace.stats.starttime.strftime("%Y-%m-%d:%H:%M:%S.%f"))
                                trace_dict['endTime'].append(trace.stats.endtime.strftime("%Y-%m-%d:%H:%M:%S.%f"))
                                logger.info(f"Processed {trace.stats.network}.{trace.stats.station}.{trace.stats.channel} ({trace.stats.npts} samples)")
                            except Exception as e:
                                logger.error(f"Error processing trace {sta.code}.{chan.code}: {str(e)}")
                                continue

                    except Exception as e:
                        logger.warning(f"Error downloading {net.code}.{sta.code}.{chan.location_code}.{chan.code}: {str(e)}")
                        continue

        # Save trace dictionary to mat file
        mat_file = os.path.join(output_dir, f"{date_str}.mat")
        try:
            sio.savemat(mat_file, {'trace': trace_dict}, format='5', do_compression=True)
            file_size = os.path.getsize(mat_file) / 1024  # Size in KB
            logger.info(f"Saved {mat_file} ({file_size:.2f} KB), {len(trace_dict['data'])} traces")
            last_trace_dict = trace_dict  # Store last interval's data
        except Exception as e:
            logger.error(f"Failed to save {mat_file}: {str(e)}")

        current_date += increment_seconds

    return last_trace_dict

if __name__ == "__main__":
    trace_data = download_seismic_data()
    if trace_data.get('network'):
        logger.info("Data download complete.")
    else:
        logger.warning("No data downloaded.")