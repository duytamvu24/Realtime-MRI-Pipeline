# in this notebook: binning scans using spirometry and ecg
import glob
import pydicom.filewriter
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.interpolate as interp

from math import floor
import datetime 
from datetime import datetime, date
from collections import OrderedDict 

# in this notebook are functions used in various differente notebooks where Dicom images are handled 
import numpy as np
import math as mt
import pydicom as dcm
import scipy as sc
# Das Notebook MRI_Functions muss im gleichen Verzeichnis wie dieses Notebook liegen.
#%run MRI_Functions_ver2.ipynb

# convert dicom timestamps to format used by the spirometer
def convert_dicom_time(dicomtime):
    """
    Convert timestamp from dicom
    """
    dicomtime = dicomtime.strip('b\'')
    dicomhours = int(str(dicomtime)[0:2])
    dicomminutes = int(str(dicomtime)[2:4])
    dicomseconds = int(str(dicomtime)[4:6])
    dicommilliseconds = int(float(str(dicomtime)[6:])*1000)
    return datetime(2000, 1, 1, dicomhours, dicomminutes, dicomseconds, dicommilliseconds*1000)





# function to import dicom images, their given aquisition time(sorted and not sorted) and their amount as arrays
# function to import dicom images, their given aquisition time(sorted and not sorted) and their amount as arrays
def read_dicomDir(input_dir):
    listPathDicom = []
    input_slices = os.listdir(input_dir)
    listPathDicom = [os.path.join(input_dir, path) for path in input_slices]
    listPathDicom.sort()
    # create Dir to save data and pictures
    # Anzahl DicomDateien im Slice = Frames
    amountFrames = len(listPathDicom)
    RefDs = [{}]*amountFrames
    AT = np.zeros([amountFrames,2])
    AT[:,0] = range(0,amountFrames)
    RefDs[0] = dcm.dcmread(listPathDicom[0])
    ArrayDicom = np.zeros([200,200,amountFrames], dtype=RefDs[0].pixel_array.dtype)
    #np.zeros([][200,200,amountFrames]],[amountFrames])
    for filenameDCM in listPathDicom:
        ds = dcm.dcmread(filenameDCM)
        #Das ist hhmmss
        AT[listPathDicom.index(filenameDCM),1]= float(ds.AcquisitionTime)
        ArrayDicom[:, :, listPathDicom.index(filenameDCM)] = ds.pixel_array # store the raw image data
    #Ordnen des gesamten Datensatzes nach der AT
    sortedArr = AT[AT[:,1].argsort()]
    sortedDicom = ArrayDicom[:,:,AT[:,1].argsort()]
    return (sortedDicom,sortedArr,AT,amountFrames)

# function to convert the aquisition time, aquired through dicom tags, into miliseconds after midnight time format

def timeConverter(time_AT):
    hours = mt.floor(time_AT/10000)
    hours_ms = hours*60*60*1000
    min =  mt.floor((time_AT-10000*hours)/100)
    min_ms = min*60*1000
    sec = mt.floor(time_AT-(10000*hours)-(100*min))
    sec_ms = sec*1000
    ms = mt.floor((time_AT-(10000*hours)-(100*min)-sec)*1000)
    com_ms = ms+sec_ms+min_ms+hours_ms
    return com_ms

# helper function to pick the three closest elements to time in list timestamps
def pick_three_closest_timestamps(timestamps, time):
    """return indices of three closest elements to time in timestamps"""
    timestamps = np.asarray(timestamps) 
    idx = np.argpartition(np.abs(timestamps - time), 3)
    return idx[:3] 


def read_spiro_data(name_spiro):
    name_spiro = name_spiro
    
    file = open(name_spiro, 'r')
    lines = file.read().splitlines()
    file.close()
    # trennt die /t voneinander, wechselt , in ein punkt
    spiro = []
    for line in lines:
        if not line:
            continue
        lineparts = line.split('\t')
        spiro.append([float(i.replace(',','.')) for i in lineparts])
    
    spiro_resorted = list(map(list, zip(*spiro)))
    # Liesst die Daten aus den drei Spiro-Dateien aus
    file = open(name_spiro, 'r')
    lines = file.read().splitlines()
    file.close()
    return spiro_resorted

# helper function that creates folders
def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

# fill empty bins by picking scans from same ecg group but different volume group for exspiration
def find_closest_from_other_bin(mid, scans_sort_by_ecg,vol_group, slc, ecg_group):
    vol_groups = list(scans_sort_by_ecg[slc][ecg_group][1].keys())
    distances = None
    for vol_group in vol_groups:
        group = scans_sort_by_ecg[slc][ecg_group][1][vol_group]
        volumes = [float(scan.get_item((0x08, 0x04)).value) for scan in group]
        if len(volumes)>0:
            final_group = group.copy()
            distances = np.array(volumes)-mid
    print(vol_group)
    return final_group[np.argmin(distances)] if distances is not None else None

def prepare_scans(path):
    scans = {} # dictionary containing a dictonary of dicom files: slice, scan 
    i = 0
    dirs = os.listdir(path)
    for folder in dirs:
        # open each folder and load all dicom files
        # place these dicom files in a dictionary
        # add dictionary to dictionary of all scans
        current = path+folder 
        scan = {}
        if folder=='.DS_Store':
            # omit system files 
            continue
        
        print(folder)
        for file in os.listdir(current):
            if file=='.DS_Store':
                continue
            i += 1
            current_file = current+'/'+file
            scan[file] = pydicom.dcmread(current_file, force=True)
        print("Read files from", current)
        scans[folder] = scan
    print("Read", i, "files")

def berechne_total_cutoff(hh, minuten, sek, ms):
    h_timestamp_spiro = (hh * 3600000)
    min_timestamp_spiro = (minuten * 60000)
    s_timestamp_spiro = (sek * 1000)
    ms_timestamp_spiro = ms
    start_spiro_time = h_timestamp_spiro + min_timestamp_spiro + s_timestamp_spiro + ms_timestamp_spiro
    
    #start_spiro_time = time_stamp - (stamp_spiro_no * SR_spiro)
    # Wie viel muss von der Spirodatei weggeschnitten werden, damit wir die differenz zwischen erstem bild
    # und spirostart haben (ohne korrelation berücksichtigt)
    cutoff = int((msnmn-start_spiro_time)/(SR_spiro))
    
    # Schneide so viele Punkte weg, bis zur Ruhephase aus der Spirometriedatei VOM ZEITSTEMPEL
    print(f"Cutoff Zeitstempel in 8ms bis vom Zeitstempel bis zum Start der Bilder: {cutoff}")
    # Auf Cutoff kommt noch die Verschiebung aus Korrelation, da diese hier noch nicht berücksichtigt wurde, also von Anfang bis Zeitstempel
    indizes = np.load('indizes.npy')
    korrelation_verschiebung = indizes[1]
    print(indizes)
    print(f"Abschneiden der Vol Daten vom Anfang bis zum Zeitstempel {korrelation_verschiebung}")
    total_cutoff = cutoff + korrelation_verschiebung
    return total_cutoff


def get_cutted_spiro_data(data_vol, hh, minuten, sek, ms):
    from datetime import datetime, timedelta
    # respiratory flow was measured every 10ms 
    data_vol = np.load("spiro.npy")
    # Erstellt erstmal nur die Timestamps, hat bisher noch nichts selber mit den Tags zu tun
    
    # Gebe hier die Uhrzeit an 
    total_cutoff = berechne_total_cutoff(hh, minuten, sek, ms)
    
    # Theoretisch reicht es auch, einfach die Anfangszeit der Spirodaten einzugeben, die Funktion später filtert nach genauer Uhrzeit
    # und nimmt die drei nahesten Punkte
    # convert flow to ml/sec
    data_vol_cutoff = data_vol[total_cutoff:]
    return data_vol_cutoff


def flow_volume_into_scans(scans_sort, timestamps, data_flow, data_vol_cutoff):
    # pick the three timestamps from the spirometry that are clostest to the acquisition time of the scan
    # calculate the median of the flow and volume measurements belonging to these three time stamps
    # save flow and respiratory volume median into dicom tags
    # flow to (0x08, 0x03), volume to (0x08, 0x04), für jede dicom, l = [trigger time, volume, flow]
    # ecg_v_info ist liste aus slice einträgen, in jedem eintrag sind l drin
    k = []
    for i in scans_sort:
        l = []
        print(i)
        for j in scans_sort[i]:
            current_scan = scans_sort[i][j]
            #print(current_scan)
            # time of current_scan
            time = convert_dicom_time(str(current_scan.get_item((0x08, 0x32)).value).strip('b\''))
            threenn = pick_three_closest_timestamps(timestamps, time)
            #flow 
            f = np.median(data_flow_cutoff[threenn])
            scans_sort[i][j].add_new((0x08, 0x03),'LO', str(f))
            # volume
            v = np.median(data_vol_cutoff[threenn])
            scans_sort[i][j].add_new((0x08, 0x04),'LO', str(v))
            # infos for plotting
            l += [[float(str(current_scan.get_item((0x18, 0x1060)).value).strip('b\'')), float(v), float(f)]]
        k += l
    ecg_v_f = np.array(k)
    #l = [[time, v1, f1], [time2, v2, l2,], [time3, v3, l3]]
    return scans_sort, ecg_v_f

from datetime import datetime

def dicomtime_first_image(dicom_time_str):
    # Beispiel: '2000-01-01 18:18:20.805000'
    dt = datetime.strptime(dicom_time_str, "%Y-%m-%d %H:%M:%S.%f")
    # Millisekunden seit Mitternacht:
    ms_since_midnight = (
        dt.hour * 3600000 +
        dt.minute * 60000 +
        dt.second * 1000 +
        dt.microsecond // 1000
    )
    return ms_since_midnight