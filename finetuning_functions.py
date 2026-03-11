# in this notebook: binning scans using spirometry and ecg
import glob
import pydicom
import pydicom.filewriter
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.interpolate as interp
from datetime import datetime, timedelta, date
from math import floor
from collections import OrderedDict 

# in this notebook are functions used in various differente notebooks where Dicom images are handled 
import numpy as np
import os
import math as mt
import pydicom as dcm
import scipy as sc
# Das Notebook MRI_Functions muss im gleichen Verzeichnis wie dieses Notebook liegen.
#%run MRI_Functions_ver2.ipynb



def read_dicomDir(input_dir):
    """
    function to import dicom images, their given aquisition time(sorted and not sorted) and their amount as arrays
    Input: Directory where the transcripted Images are located
    """
    listPathDicom = []
    input_slices = os.listdir(input_dir)
    listPathDicom = [os.path.join(input_dir, path) for path in input_slices]
    listPathDicom.sort()
    # create Dir to save data and pictures
    amountFrames = len(listPathDicom)
    RefDs = [{}]*amountFrames
    AT = np.zeros([amountFrames,2])
    AT[:,0] = range(0,amountFrames)
    print(listPathDicom[0])
    RefDs[0] = dcm.dcmread(listPathDicom[0])
    ArrayDicom = np.zeros([200,200,amountFrames], dtype=RefDs[0].pixel_array.dtype)
    #np.zeros([][200,200,amountFrames]],[amountFrames])
    for filenameDCM in listPathDicom:
        ds = dcm.dcmread(filenameDCM)
        #Das ist hhmmss
        AT[listPathDicom.index(filenameDCM),1]= float(ds.AcquisitionTime)
        ArrayDicom[:, :, listPathDicom.index(filenameDCM)] = ds.pixel_array # store the raw image data
    sortedArr = AT[AT[:,1].argsort()]
    sortedDicom = ArrayDicom[:,:,AT[:,1].argsort()]
    return (sortedDicom,sortedArr,AT,amountFrames)

# 
def timeConverter(time_ms):
    """
    function to convert the aquisition time, aquired through dicom tags, into miliseconds after midnight time format
    Input: Acquisition time from MR-Images
    """
    seconds = time_ms // 1000
    milliseconds = time_ms % 1000

    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    remaining_seconds = seconds % 60

    time = timedelta(hours=hours, minutes=minutes, seconds=remaining_seconds, milliseconds=milliseconds)
    return str(time)



def read_bellow(filename):
    # read out bellow file and defining the first value that is part of the respiratory bellows curve
    with open(filename) as f:
        lines=f.readlines()
        for line in lines[1:]:
            if "LogStartMDHTime:" in line:
                log_start_time = int(line.split(":")[1].strip())  # Zahl extrahieren
                print("LogStartMDHTime:", log_start_time)
                real_time_bellow = timeConverter(log_start_time)
                break  # Falls du nur den ersten Treffer brauchst
            
        line = lines[0][9:]
        line = [int(s) for s in line.split() if s.isdigit()]
        line = np.array(line)
        line = line[np.logical_and(line!=6000, line!=5000)]
        to_delete = list(zip(np.where(line==5002)[0], np.where(line==6002)[0]))
        to_delete_idxs = np.concatenate([np.arange(*tup) for tup in to_delete])
        line = np.delete(line, to_delete_idxs)
        data = np.delete(line, np.where(line==6002)[0])
        #exporting the data for the respiratory bellows curve as "cleaned" data
        np.savetxt(f'{filename}_cleaned.resp',data)
        return data, real_time_bellow

def cut_bellow_data(sig_rb, AT, real_time_bellow):
    """
    Cuts the respiratory bellow signal so that it matches
    the time range of the MRI slices.
    """
    # Get the timestamp of the first MRI slice
    time_first_slice = str(AT[0][1])
    time_obj_first_slice = datetime.strptime(time_first_slice, "%H%M%S.%f")

    # Format the time for display
    formatted_time_first_slice = time_obj_first_slice.strftime("%H:%M:%S.%f")[:-3]
    print(formatted_time_first_slice)

    # Get the timestamp of the last MRI slice
    time_last_slice = str(AT[-1][1])
    last_obj_last_slice = datetime.strptime(time_last_slice, "%H%M%S.%f")

    # Format the time for display
    formatted_last_time_last_slice = last_obj_last_slice.strftime("%H:%M:%S.%f")[:-3]
    # Convert the start time of the bellow signal to datetime
    real_time_bellow_formatted = datetime.strptime(real_time_bellow, '%H:%M:%S.%f')
    # Sampling interval of the bellow signal (seconds per sample)
    bellow_sr = 0.0025
    # Calculate the time difference between bellow start and first MRI slice
    delta_start_to_first_slice = time_obj_first_slice - real_time_bellow_formatted
    print("delta: " + str(delta_start_to_first_slice))
    # Convert time difference to seconds
    seconds_first_slice_bellow = delta_start_to_first_slice.total_seconds()
    # Calculate time difference between bellow start and last MRI slice
    seconds_last_slice_bellow = (last_obj_last_slice - real_time_bellow_formatted).total_seconds()
    # Convert the time differences into sample indices
    bellow_time_stamps_to_first_slice = int(seconds_first_slice_bellow / bellow_sr)
    bellow_time_stamps_to_last_slice = int(seconds_last_slice_bellow / bellow_sr)
    # Cut the bellow signal so it matches the MRI time window
    sig_rb = sig_rb[bellow_time_stamps_to_first_slice:bellow_time_stamps_to_last_slice]
    print(f"Difference between first slice and bellow start: {seconds_first_slice_bellow:.3f} seconds")
    return sig_rb

def cut_spiro(indizes_filename, spiro_npy_filename, AT, zeit_start_zeitstempel_spiro, manuelle_verschiebung_s):
    """
    Cuts the spirometry volume data so that it matches the time range
    of the MRI slice acquisition.
    """
    # Load data
    indizes = np.load(f'{indizes_filename}')
    index_spiro = indizes[1]
    spiro_vol = np.load(f'{spiro_npy_filename}')
    # Cut the data starting from the predefined index
    spiro_vol = spiro_vol[index_spiro:]
    # Get the timestamp of the first MRI slice
    time_first_slice = str(AT[0][1])
    time_obj_first_slice = datetime.strptime(time_first_slice, "%H%M%S.%f")
    # Format for printing
    formatted_time_first_slice = time_obj_first_slice.strftime("%H:%M:%S.%f")[:-3]
    print(formatted_time_first_slice)
    # Get timestamp of the last MRI slice
    time_last_slice = str(AT[-1][1])
    last_obj_last_slice = datetime.strptime(time_last_slice, "%H%M%S.%f")
    formatted_last_time_last_slice = last_obj_last_slice.strftime("%H:%M:%S.%f")[:-3]
    print(f"Cutting off spiro data from beginning up to timestamp {index_spiro}")
    # Convert timestamps to datetime objects
    t1 = datetime.strptime(formatted_time_first_slice, "%H:%M:%S.%f")
    print(f"The first slice begins at {t1}")
    t2 = datetime.strptime(zeit_start_zeitstempel_spiro, "%H:%M:%S")
    t3 = datetime.strptime(formatted_last_time_last_slice, "%H:%M:%S.%f")
    # Calculate time difference between spirometry start and first MRI slice
    diff = t1 - t2
    diff_seconds = diff.total_seconds()
    print(f"Difference: {diff_seconds} seconds")
    # Calculate time duration between first and last MRI slice
    diff_zeit_seconds_between_first_last = (t3 - t1).total_seconds()
    print(f"Time span of slice volume acquisition: {diff_zeit_seconds_between_first_last}")
    # Sampling interval of spirometry data (8 ms)
    spiro_sr = 0.008
    # Convert time difference to number of samples
    anzahl_samples = diff_seconds / spiro_sr
    print(int(anzahl_samples))
    # Manual time shift to correct synchronization
    manuelle_verschiebung_timesteps = manuelle_verschiebung_s / spiro_sr
    # Number of samples covering the MRI acquisition period
    spiro_slice_timesteps = diff_zeit_seconds_between_first_last / spiro_sr
    # Number of samples between spiro start and first slice
    sig_curve_vol_timestamps = diff_seconds / spiro_sr
    print(f"Total cut position: {sig_curve_vol_timestamps + index_spiro}")
    # Cut the spirometry volume data to match MRI acquisition time
    volume_data_spiro = spiro_vol[
        int(sig_curve_vol_timestamps + manuelle_verschiebung_timesteps):
        int(sig_curve_vol_timestamps + manuelle_verschiebung_timesteps + spiro_slice_timesteps)
    ]
    print(f"Length of spiro data: {len(volume_data_spiro)} samples")
    print(f"Duration in seconds: {len(volume_data_spiro) * 0.008}")
    return volume_data_spiro

def berechne_zeitinformationen(acquisition_times, startzeit_spiro):
    """
    Computes timing information between the start of the spirometry recording
    and the MRI slice acquisition.

    Parameters
    ----------
    acquisition_times : list
        List containing [slice_index, timestamp] where the timestamp is in the format HHMMSS.fff

    startzeit_spiro : str
        Start time of the spirometry recording in the format HH:MM:SS

    Returns
    -------
    float
        Time difference in seconds between the first and the last MRI slice
    """
    # Extract raw timestamps of the first and last MRI slice
    zeit_erster_slice_raw = str(acquisition_times[0][1])
    zeit_letzter_slice_raw = str(acquisition_times[-1][1])
    zeit_obj_erster_slice = datetime.strptime(zeit_erster_slice_raw, "%H%M%S.%f")
    zeit_obj_letzter_slice = datetime.strptime(zeit_letzter_slice_raw, "%H%M%S.%f")
    zeit_erster_slice_fmt = zeit_obj_erster_slice.strftime("%H:%M:%S.%f")[:-3]
    zeit_letzter_slice_fmt = zeit_obj_letzter_slice.strftime("%H:%M:%S.%f")[:-3]
    zeit_start_spiro_obj = datetime.strptime(startzeit_spiro, "%H:%M:%S")
    # Calculate time difference between spirometry start and first MRI slice
    differenz_start_spiro_zu_erstem_slice = (
        zeit_obj_erster_slice - zeit_start_spiro_obj
    ).total_seconds()

    # Calculate time span between first and last MRI slice
    differenz_erster_zu_letzter_slice = (
        zeit_obj_letzter_slice - zeit_obj_erster_slice
    ).total_seconds()
    print(f"First slice starts at {zeit_erster_slice_fmt}")
    print(f"Time difference between spirometry start and first slice: {differenz_start_spiro_zu_erstem_slice:.3f} s")
    print(f"Time span between first and last slice: {differenz_erster_zu_letzter_slice:.3f} s")
    return differenz_erster_zu_letzter_slice

# 
def pick_three_closest_timestamps(timestamps, time):
    """helper function to pick the three closest elements to time in list timestamps"""
    timestamps = np.asarray(timestamps) 
    idx = np.argpartition(np.abs(timestamps - time), 3)
    return idx[:3] 




def convert_dicom_time(dicomtime):
    """
    convert dicom timestamps to format used by the spirometer
    """
    dicomtime = dicomtime.strip('b\'')
    dicomhours = int(str(dicomtime)[0:2])
    dicomminutes = int(str(dicomtime)[2:4])
    dicomseconds = int(str(dicomtime)[4:6])
    dicommilliseconds = int(float(str(dicomtime)[6:])*1000)
    return datetime(2000, 1, 1, dicomhours, dicomminutes, dicomseconds, dicommilliseconds*1000)
