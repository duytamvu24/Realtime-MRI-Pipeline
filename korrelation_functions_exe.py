# clean ECG
import numpy as np
from scipy import signal
#import ipywidgets as widgets
import re
import math as mt
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


def clean_data(line):
    c = 1
    data = []
    skip = False
    for l in line:
        if int(l) == 5002:
            skip = True
        #print(c, l, skip)
        if not skip and l != 5000 and l!=6000:
            data.append(l)
        c = c + 1
        if skip:
            if int(l) == 6002:
                skip = False
    # Convert the list to a NumPy array
    data = np.array(data)
    return data
    
from datetime import timedelta

def timeConverter(time_ms: int) -> str:
    # Total seconds
    seconds = time_ms // 1000
    # Remaining milliseconds
    ms = time_ms % 1000

    # Convert into hours, minutes and seconds
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    remaining_seconds = seconds % 60

    # timedelta used to generate hh:mm:ss format
    base_time = timedelta(hours=hours, minutes=minutes, seconds=remaining_seconds)

    # Append milliseconds to the formatted string
    return f"{base_time}.{ms:03d}"


def read_ecg(filename):
    with open(filename) as f:
        lines=f.readlines()
        line = lines[0]
        line = [int(s) for s in line.split() if s.isdigit()]
        data = clean_data(line)
        for line in lines[1:]:
            if "LogStartMDHTime:" in line:
                log_start_time = int(line.split(":")[1].strip())  # Extract numeric value
                print("LogStartMDHTime:", log_start_time)
                real_time_bellow = timeConverter(log_start_time)
                print(f"ECG recording start time: {real_time_bellow}")
                break  # Stop after first occurrence
        return data, log_start_time


def read_spiro_data(name_spiro):
    name_spiro = name_spiro
    
    file = open(name_spiro, 'r')
    lines = file.read().splitlines()   
    file.close()

    # Split columns using tab separator and convert commas to decimal points
    spiro = []
    for line in lines:
        if not line:
            continue
        lineparts = line.split('\t')
        spiro.append([float(i.replace(',','.')) for i in lineparts])
    
    spiro_resorted = list(map(list, zip(*spiro)))

    # Read spirometry data again (structure preserved from original implementation)
    file = open(name_spiro, 'r')
    lines = file.read().splitlines()
    file.close()
    return spiro_resorted




from scipy.signal import find_peaks

def calc_correction(data, start, end):
    inspiration = (data[:] >= 0)
    expiration = np.logical_not(inspiration)

    # Integrate positive (inspiration) and negative (expiration) parts separately
    insp_int = np.sum(data[start:end][inspiration[start:end]])
    exp_int = np.sum(-1*data[start:end][expiration[start:end]])

    # Calculate correction factor
    faktor = insp_int/exp_int
    data_corr = data.copy()

    # Apply correction factor only to inspiration values
    data_corr[:][inspiration] = data[:][inspiration]/faktor

    return (data_corr, insp_int, exp_int, faktor)


def get_volume_from_flow(flow):
    data_corr, insp_int, exp_int, faktor = calc_correction(np.array(flow), 0,-1)

    # Compute cumulative sum to obtain volume from flow
    vol_corr = np.cumsum(data_corr[:])
    
    # Automatic baseline correction according to Halima
    vol_automatic_BC = BC_vol(vol_corr,0, -1)

    # Convert to liters using sampling rate
    new_vol_automatic_BC = vol_automatic_BC * 8 / 1000
    return new_vol_automatic_BC


# Automatic baseline correction if volume has already been calculated
def BC_vol(vol,begin,end):

    peaks_min_idx0, peaks_min = calc_all_minima(vol[begin:end],0,-1)

    # Combine peak values and indices into nested list
    peaks_and_idx = np.vstack((peaks_min, peaks_min_idx0)).T

    # Detect outliers among detected minima
    list_min0=detect_outliers(peaks_and_idx) 
    print(len(list_min0))

    new_peaks_min = np.delete(peaks_min,list_min0)

    # Indices of minima that were not removed
    new_peaks_min_idx = np.delete(peaks_min_idx0,list_min0)

    peak_inter_med = median_window(2,vol[begin:end], new_peaks_min_idx)

    vol_intervall_new = vol[begin:end]-peak_inter_med

    return vol_intervall_new


def calc_all_minima(vol,result_begin,result_end):

    vol_med = vol.copy()[result_begin:result_end]

    # Detect minima by applying peak detection to the inverted signal
    peaks_min_idx0= find_peaks(-vol[result_begin:result_end],distance = 100)[0]

    # Extract corresponding volume values at the detected minima
    peaks_min = vol[peaks_min_idx0 +result_begin]

    print(len(peaks_min))

    return (peaks_min_idx0,peaks_min)


def detect_outliers(peaks_and_idx):
    
    list_min0 =[]
    k=0
    j = 0

    for idx,peaks in enumerate(peaks_and_idx):

        if idx==0:
            if peaks[0]>600:
                # If first peak value is larger than threshold remove it
                list_min0.append(idx)

            if peaks[0]-peaks[1]>100:  # baseline correction towards lower values
                list_min0.append(idx)

        if idx>0:

            diff_time = peaks[1]-k
            diff_peaks = peaks[0]-j

            if diff_time<50 or diff_peaks>50:
                list_min0.append(idx)

        k = peaks[1]
        j = peaks[0]

    return list_min0


def median_window(k,vol_intervall,new_peaks_min_idx):

    # Copy of the volume signal
    vol_med = vol_intervall.copy()

    n = len(new_peaks_min_idx)

    xnew = np.linspace(0,vol_med.size+1 , num=vol_med.size, endpoint=True)

    # Replace each peak value with the median of its neighborhood
    for i in range(k, n-k-1):

        neighborhood = new_peaks_min_idx[(i-k):(i+k+1)] 

        vol_med[new_peaks_min_idx[i]] = np.median(vol_med[neighborhood])

    # Linear interpolation between corrected peak values
    peak_inter_med = np.interp(xnew, new_peaks_min_idx,vol_med[new_peaks_min_idx])

    return (peak_inter_med)
