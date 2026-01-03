#clean ECG
import numpy as np
from scipy import signal
#import ipywidgets as widgets
from IPython.display import display
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
    #line = [int(s) for s in line.split() if s.isdigit()]
    data = np.array(data)
    return data
    
from datetime import timedelta

def timeConverter(time_ms: int) -> str:
    # Gesamtsekunden
    seconds = time_ms // 1000
    # Rest-Millisekunden
    ms = time_ms % 1000

    # Stunden, Minuten, Sekunden
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    remaining_seconds = seconds % 60

    # timedelta nur für hh:mm:ss
    base_time = timedelta(hours=hours, minutes=minutes, seconds=remaining_seconds)

    # String mit Millisekunden hinten anhängen
    return f"{base_time}.{ms:03d}"


def read_ecg(filename):
    with open(filename) as f:
        lines=f.readlines()
        line = lines[0]
        line = [int(s) for s in line.split() if s.isdigit()]
        data = clean_data(line)
        for line in lines[1:]:
            if "LogStartMDHTime:" in line:
                log_start_time = int(line.split(":")[1].strip())  # Zahl extrahieren
                print("LogStartMDHTime:", log_start_time)
                real_time_bellow = timeConverter(log_start_time)
                print(f"Uhrzeit EKG Beginn: {real_time_bellow}")
                break  # Falls du nur den ersten Treffer brauchst
        return data, log_start_time

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


def identify_ecg(data):
    #%matplotlib widget
    import numpy as np
    import matplotlib.pyplot as plt
    from IPython.display import display
    # Interaktive Darstellung in Jupyter
    
    #Wie viel Prozent vom Datensatz willst du sehen? Optimal, um feinere Punkte zu sehen
    perc = 0
    perc_end = 0.65
    # Beispiel-Daten für das Signal  
    y_data = data[6::4]
    perc_datenpunkt = int(perc * len(y_data))
    perc_datenpunkt_end = int(perc_end * len(y_data))
    y_data = y_data[perc_datenpunkt:perc_datenpunkt_end]
    x_data = np.arange(len(y_data)) * 0.0025
    
    # Liste für gespeicherte Punkte
    time_stamp_ecg = []
    
    # Default X-Range für den Slider
    x_min, x_max = 0, 3000
    
    def find_nearest_point(x_click):
        """Findet den nächsten X-Wert in den Daten und gibt X und Y zurück."""
        idx = (np.abs(x_data - x_click)).argmin()  # Index des nächsten Punktes
        return idx, y_data[idx]
    
    def onclick(event):
        """Setzt einen Punkt an der nächstgelegenen X-Position der Daten."""
        if event.xdata is not None and event.ydata is not None:
            nearest_x, nearest_y = find_nearest_point(event.xdata)
            time_stamp_ecg.append((nearest_x, nearest_y))
            print(f"Punkt gesetzt bei: ({nearest_x:.2f}, {nearest_y:.2f})")
            update_plot()
    
    def update_plot(change=None):
        """Zeichnet das Signal mit allen gesetzten Punkten neu."""
        ax.clear()  # Löscht den aktuellen Plot
        ax.plot(x_data, y_data, label="ECG Signal")
        
        # Alle gespeicherten Punkte einzeichnen
        if time_stamp_ecg:
            x_pts, y_pts = zip(*time_stamp_ecg)  # Listen mit X- und Y-Werten
            ax.scatter(x_data[x_pts], y_pts, color="red", s=100, label="Manuelle Punkte")
    
        ax.set_xlim(x_range_slider.value)  # X-Range basierend auf Slider setzen
        ax.legend()
        ax.set_xlabel("Zeit in s")
        ax.set_ylabel("Amplitude")
        ax.set_title("Signal mit interaktivem Bereich & manuell gewählten Punkten")
        ax.grid(True)
        fig.canvas.draw()  # Aktualisiert den Plot
    
    def initialize_plot():
        """Zeichnet das Signal mit allen gesetzten Punkten neu."""
        ax.clear()  # Löscht den aktuellen Plot
        ax.plot(x_data, y_data, label="ECG Signal")
        
        # Alle gespeicherten Punkte einzeichnen
        if time_stamp_ecg:
            x_pts, y_pts = zip(*time_stamp_ecg)  # Listen mit X- und Y-Werten
            ax.scatter(x_pts, y_pts, color="red", s=100, label="Manuelle Punkte")
    
        #ax.set_xlim(x_range_slider.value)  # X-Range basierend auf Slider setzen
        ax.legend()
        ax.set_xlabel("Zeit in s")
        ax.set_ylabel("Amplitude")
        ax.set_title("Signal mit interaktivem Bereich & manuell gewählten Punkten")
        ax.grid(True)
        fig.canvas.draw()  # Aktualisiert den Plot
    
    # Erstellt die Figur mit figsize (16, 4)
    fig, ax = plt.subplots(figsize=(10, 4))
    initialize_plot()
    #update_plot()
    
    # Event-Listener für Mausklicks hinzufügen
    fig.canvas.mpl_connect("button_press_event", onclick)
    
    
    
    # Erstelle interaktiven Slider für X-Range
    x_range_slider = widgets.FloatRangeSlider(
        value=[x_min, x_max], 
        min=0, 
        # Bis zu welchem Datenpunkt wird gezeigt -> es soll alles gezeigt werden
        max=x_data[-1], 
        step=2, 
        description='X-R',
        layout=widgets.Layout(width='100%')
    )
    
    # Aktualisiert den Plot, wenn sich der Slider bewegt
    x_range_slider.observe(update_plot, names='value')
    
    # Zeigt den Slider an
    display(x_range_slider)



from scipy.signal import find_peaks

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

def get_volume_from_flow(flow):
    data_corr, insp_int, exp_int, faktor = calc_correction(np.array(flow), 0,-1)
    vol_corr = np.cumsum(data_corr[:]) # komplett#
    
    #Baseline Korrektur nach Halima
    vol_automatic_BC = BC_vol(vol_corr,0, -1)
    new_vol_automatic_BC = vol_automatic_BC * 8 / 1000
    return new_vol_automatic_BC

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

def calc_all_minima(vol,result_begin,result_end):
    vol_med = vol.copy()[result_begin:result_end]
    peaks_min_idx0= find_peaks(-vol[result_begin:result_end],distance = 100)[0] # Indizes berechnen, [0 ist der array selber]
    peaks_min = vol[peaks_min_idx0 +result_begin] # Werte aus dem Volumenarray in peaks_min abspeichern 
    print(len(peaks_min))
    return (peaks_min_idx0,peaks_min)

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