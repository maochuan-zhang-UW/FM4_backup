import obspy
from obspy.clients.fdsn.client import Client

# Load the .mseed file
file_path = '/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/output_test.mseed'
st = obspy.read(file_path)

# Filter for stations matching 2F.AX**A with channels either "HH?" or "ELZ"
filtered_st = st.select(station="AX???").select(channel="HH?") + st.select(station="AX???", channel="EL?")

# Define the start and end times for trimming
start_time = obspy.UTCDateTime("2022-10-25T10:00:00")
end_time = obspy.UTCDateTime("2022-10-25T11:00:00")

# Trim the stream to the specified time range
trimmed_st = filtered_st.trim(starttime=start_time, endtime=end_time)

# Normalize the trimmed stream
trimmed_st.normalize()

# Plot the normalized and trimmed waveform
trimmed_st.plot()
# Print the structure of trimmed_st
print("Number of traces in trimmed_st:", len(trimmed_st))
print("\nDetails of each trace in trimmed_st:")
for i, trace in enumerate(trimmed_st):
    print(f"  Station: {trace.stats.station}")
    print(f"  Channel: {trace.stats.channel}")
