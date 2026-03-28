#Importieren der Ateminformation aus der ASC Datei, die von JScope generiert wird
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import scipy as sc
from datetime import date, datetime, timedelta


#Datum:
day = date.today()
dateiname_day = day.isoformat()
#rcParams['figure.figsize'] = 10, 5

# Bestimmung der lokalen Minima zu Beginn und zum Ende der Spirometriekurve.
# dafür wird ein Wert, ab einem bestimmmten Startzeitpunkt, für den folgendes gilt:
# der Flow ist zu diesem Zeitpunkt=0
# die letzten 10 Werte (10 ist hier definiert mit der Variable k) vor dem Zeitpunkt entsprechen einer Exspiration,
# sie sind aufsummiert also kleiner als 0
# die ersten 10 Werte (10=k) nach dem Minimum sind aufsummiert größer als 0, sie entsprechen also einer Inspiration
# das bestimmte Minimum zu Beginn ergibt den Wert 5208, dieser wird als Kontrolle in der Grafik angezeigt.
# die Grafik zeigt den Flow in l/s (y-Achse) zwischen 5000 und 5300 (x-Achse, entspricht 0-300)
# result_begin ist der Beginn des Intervalls in dem der Drift bearbeitet wird
def find_local_minimum(data, begin, k=10):
    data = data[begin:]
    idx = 10
    while data[(idx-k):idx,0].sum() >= 0 or not np.isclose(data[idx, 0],0) or data[idx:(idx+k),0].sum() <= 0:
        idx +=1
    return idx + begin


# Klassifikation der Werte für den Flow in Inspiration falls >=0 und in Exspiration falls <=0
# Bestimmung des Integrals für alle inspiratorischen Flowwerte= Atemvolumen der Inspiration
# Bestimmung des Integrals für alle exspiratorischen Flowwerte= Atemvolumen der Exspiration
# da diese nicht gleich sind, muss eine Korrektur durchgeführt werden, dafür wird das inspiratorische Volumen
# durch das exspiratorische Volumen geteilt, das Ergebnis ist der ermittelte Korrekturfaktor für die Inspiration
# Erstellung eines neuen Datensatzes für den Flow:
# die inspiratorischen Atemflusswerte werden dafür durch den Korrekturfaktor geteilt
# die exspiratorischen Atemflusswerte werden dafür unverändert übernommen
#data_corr alle korrigierten Daten, nur insp, nur exp, faktor

# Schema des Eintragens: data_corr[anfang_BC:ende_BC,Spaltenvariation]

def calc_correction(data, start, end):
    inspiration = (data[:] >= 0)
    expiration = np.logical_not(inspiration)
    insp_int = np.sum(data[start:end][inspiration[start:end]])
    exp_int = np.sum(-1*data[start:end][expiration[start:end]])
    #Integral nur der Werte über und unter null
    faktor = insp_int/exp_int
    data_corr = data.copy()
    #Korrekturfaktor nur auf Inspiration
    data_corr[:][inspiration] = data[:][inspiration]/faktor
    return (data_corr, insp_int, exp_int, faktor)

#1.Minima berechen:
#Die Minima berechnen sich aus den Maxima der um X gespiegelten Volumendatei in dem vorgegebenen Intervall
#In Peaks_min befinden sich dann alle Minima
#Alle Peaks entfernen, die nicht mehr benötigt werden
#distance = 100, mind. 400 ist ein atemzug
def calc_all_minima(vol,result_begin,result_end):
    vol_med = vol.copy()[result_begin:result_end]
    peaks_min_idx0= find_peaks(-vol[result_begin:result_end],distance = 100)[0] # Indizes berechnen, [0 ist der array selber]
    peaks_min = vol[peaks_min_idx0 +result_begin] # Werte aus dem Volumenarray in peaks_min abspeichern 
    print(len(peaks_min))
    return (peaks_min_idx0,peaks_min)

def calc_all_maxima(vol,result_begin,result_end):
    vol_med = vol.copy()[result_begin:result_end]
    peaks_max_idx0,_= np.ravel(find_peaks(vol[result_begin:result_end])) # Indizes berechnen
    peaks_max = vol[peaks_max_idx0 +result_begin] # 
    return (peaks_max_idx0,peaks_max) 


#2.Median Windowing + lineare Interpolation 
#Minima glätten mittels gleitender Medianwert Berechnung:
#lokale Minima mittles Median Glätten:
#Für jedes lokale Minimum Min(k) wird folgender Wert berechnet:
#von z.B. den k lokalen Minima vor und k lokalen Minima nach und eben diesem lokalen Minimum wird das Mediumberechnet und als Wert übernommen.

#lineare Interpolation: Peaks auf die Menge an Werten die Baselinekorrigiert werden interpolieren
def median_window(k,vol_intervall,new_peaks_min_idx):
    #Kopie des Volumens
    vol_med = vol_intervall.copy()

    n = len(new_peaks_min_idx)
    xnew = np.linspace(0,vol_med.size+1 , num=vol_med.size, endpoint=True)
    #Jeder Peakwert wird ersetzt, und zwar durch einen Wert, der benachbarten Werte um die Peaks der median davon
    for i in range(k, n-k-1):
        neighborhood = new_peaks_min_idx[(i-k):(i+k+1)] 
        vol_med[new_peaks_min_idx[i]] = np.median(vol_med[neighborhood])
    peak_inter_med = np.interp(xnew, new_peaks_min_idx,vol_med[new_peaks_min_idx])# linear interpoliertes Array der Medianwerte

    return (peak_inter_med)


#4.Letzte Outlier mit einem Tiefpassfilter(butterworth) korrigieren in den lin interpolierten Minimas
#Abtastrate = 100Hz
#Grenzfrequenz 0.025 Hz
#nyqist Frequenz: 1/2* Abtastrate (Abtasttherorem einhalten)
#normale Grzfrequenz: Grnezfrequenz/ nyqistfrequenz
def lowpass_Filter(peak_inter_med):
    samplingRate = 100
    cutoff = 0.05
    nyq = 0.5*samplingRate
    order = 2
    #data_med = peak_inter_med
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = sc.signal.butter(order, normal_cutoff, btype='low', analog=False) # Filteraufstellen
    data_fil_med = sc.signal.filtfilt(b,a,peak_inter_med)# Filter anwenden
    return data_fil_med


# diese Funktion entdeckt die Outlier die fälschlicherweise bei der Bestimmung der Minima gefunden wurden.
def detect_outliers(peaks_and_idx):#
    
    list_min0 =[]
    k=0
    j = 0
    for idx,peaks in enumerate(peaks_and_idx):
        if idx==0:
            if peaks[0]>600:
                #Wenn größer als 600 wird der erste wert hinzugefügt
                list_min0.append(idx)
            if peaks[0]-peaks[1]>100: # bc nach unten 
                list_min0.append(idx)
        if idx>0:
            diff_time = peaks[1]-k
            diff_peaks = peaks[0]-j
            if diff_time<50 or diff_peaks>50:
                list_min0.append(idx)
        k = peaks[1]
        j = peaks[0]
    return list_min0


# automtatische BC bei unbekannten Volumen und bekannten Flow
def BC(flow,begin,end):
    
    vol_com = np.cumsum(flow[:])*10
    peaks_min_idx0, peaks_min = calc_all_minima(vol_com[begin:end],0,-1)
    peaks_and_idx = np.vstack((peaks_min, peaks_min_idx0)).T
    list_min0=detect_outliers(peaks_and_idx)    
    new_peaks_min = np.delete(peaks_min,list_min0)
    new_peaks_min_idx = np.delete(peaks_min_idx0,list_min0)
    peak_inter_med = median_window(2,vol_com[begin:end], new_peaks_min_idx)
    vol_intervall_new = vol_com[begin:end]-peak_inter_med
    return vol_intervall_new


# automatische BC wenn das VOlumen bereits berechnet ist 
def BC_vol(vol,begin,end):
    peaks_min_idx0, peaks_min = calc_all_minima(vol[begin:end],0,-1)
    #Ordnet in nestled list die peak werte mit dem index zu 
    peaks_and_idx = np.vstack((peaks_min, peaks_min_idx0)).T
    #Detektiert outliers, hier einmal nachfragen, 
    list_min0=detect_outliers(peaks_and_idx) 
    print(len(list_min0))
    new_peaks_min = np.delete(peaks_min,list_min0)
    #index aller nicht gelöschten minima
    new_peaks_min_idx = np.delete(peaks_min_idx0,list_min0)
    peak_inter_med = median_window(2,vol[begin:end], new_peaks_min_idx)
    vol_intervall_new = vol[begin:end]-peak_inter_med
    return vol_intervall_new
    