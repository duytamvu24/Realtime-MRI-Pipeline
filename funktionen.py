#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.signal import find_peaks
import tkinter as tk
import json
import BC_Spiro_Functions_flow_npy_general as sm

#plt.style.use(['science', 'notebook', 'grid'])


kanal_infos = {
    1: {"y_label": "Flow [L/s]"},
    2: {"y_label": "O2-Volume fraction in %"},
    3: {"y_label": "CO2-Volume fraction in %"},
    4: {"y_label": "Volumen [L]"},
}

find_peaks_params = {
    1: {"height": 0.1, "distance": 150, "prominence": 0.3},
    2: {"height": 20, "distance": 200, "prominence": 0.8},
    3: {"height": 3.5, "distance": 200, "prominence": 0.8},
    4: {"distance" : 100, "prominence" : 0.3, "height" : 0.2}
}

find_valleys_params = {
    1: {"height": -0.1, "distance": 150, "prominence": 0.3},
    2: {"height": -17, "distance": 200, "prominence": 0.8},
    3: {"height": -0.1, "distance": 200, "prominence": 0.8},
    4: {"distance": 100, "prominence" : 0.3, "height" : 0.2}
}

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

def get_volume_from_flow(flow):
    data_corr, insp_int, exp_int, faktor = sm.calc_correction(np.array(flow), 0,-1)
    vol_corr = np.cumsum(data_corr[:]) # komplett#
    
    #Baseline Korrektur nach Halima
    vol_automatic_BC = sm.BC_vol(vol_corr,0, -1)
    new_vol_automatic_BC = vol_automatic_BC * 8 / 1000
    return new_vol_automatic_BC

def flag_belastung(daten, beginn_belastung):
    """
    Unterteilt die Daten in Belastung und Ruhephase ein, mit Hilfe eines Flags. 
    Erstellt eine Liste mit gleich vielen Elementen, mit 0 Ruhephase, und einer 1 Belastungsphase.
    Input: Anfangszeit/Index von Beginn der Belastungsphase in min/sekunden
    """
    neue_liste = [0] * len(daten)
    # Index definieren, ab dem die Werte zu 1 werden
    timesteps_to_belastung = beginn_belastung / 0.008
    neue_liste[timesteps_to_belastung:] = [1] * (len(neue_liste) - start_index)
    return neue_liste


def plot_graph(x,y,label,x_label, y_label, title,):
    plt.plot(x,y,label = label)
    plt.xlabel(x_label)
    plt.title(title)
    plt.xticks(rotation=90)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()

def time_to_real_time(time_data, str_start_time):
    # Startzeit angeben (z. B. 14:00 Uhr)
    start_time = datetime.strptime(str_start_time, "%H:%M:%S")
    minutes_list_spiro = time_data
    # Update der Liste: Minuten in Uhrzeiten umrechnen
    time_list = [start_time + timedelta(minutes=minutes) for minutes in minutes_list_spiro]
    return time_list

def filtered_data(time_list, data, str_start_time, str_end_time):
    start_time = datetime.strptime(str_start_time, "%H:%M:%S")
    end_time = datetime.strptime(str_end_time, "%H:%M:%S")
    
    #Filter die Daten auf die Anfangs und Endzeit
    filtered_x = [x_time for x_time in time_list if start_time <= x_time <= end_time]
    filtered_y = [data for i in range(len(time_list)) if start_time <= time_list[i] <= end_time]
    return filtered_x, filtered_y



def bestimme_atempausen(liste, grenze):
    #Returned eine Liste mit der selben Dimension wie die Daten zurück, wo Daten wo Atempausen sind, genullt werden, sonst 1
    # Grenze gibt dabei an, welcher Zeitraum zwischen peaks eine grenze detektiert wird
    n = len(liste)
    ergebnis = [1] * n
    start = -1  # Startindex eines gültigen Abschnitts
    
    for i in range(n):
        if -0.05 <= liste[i] <= 0.05:
            if start == -1:
                start = i  # Start eines neuen gültigen Bereichs
        else:
            if start != -1 and (i - start) > grenze:  # Prüfen, ob mehr als `grenze` Elemente passen
                for j in range(start, i):
                    ergebnis[j] = 0
            start = -1  # Zurücksetzen des Startwerts

    # Falls am Ende der Liste noch ein gültiger Bereich existiert
    if start != -1 and (n - start) > grenze:
        for j in range(start, n):
            ergebnis[j] = 0
            
    return ergebnis


def bestimme_atemfrequenz_mit_angepassten_atempausen(zeit, volume):

    # Umrechnung der Indizes der Peaks in Sekunden
    peaks, _ = find_peaks(volume, distance=200, prominence=0.3)
    atempausen = bestimme_atempausen(volume, 200)

    

    # Beispiel:
    n = 1  # Peaks links entfernen
    m = 1  # Peaks rechts entfernen
    filtered_peaks = list(peaks)

    # Finde die Indizes, an denen eine Pause vorliegt
    pause_indices = np.where(np.array(atempausen) == 0)[0]


    l = find_last_indices_of_zeros(atempausen)


    to_delete = []
    if len(l) == 0:
        print("Im angegebenen Bereich gibt es keine Atempausen!")
    else:
        #Erster Peak bestimmen, falls die Daten mit einer Pause beginnen
        if atempausen[0] == 0:
            delete_peak = [peak for peak in peaks if l[0] < peak < l[1]]
            delete_peak = delete_peak[:m] + delete_peak[-n:]
            atempausen = to_zeros(atempausen, l[0], delete_peak[m-1])
            atempausen = to_zeros(atempausen, delete_peak[-n], l[1])
            to_delete.extend(delete_peak) 
        else: 
            # Wenn nicht mit atempause anfängt, dann lösche n erste peaks linka vor der ersten Atempause
            delete_peak = [peak for peak in peaks if peak < l[0]][-n:]
            atempausen = to_zeros(atempausen, delete_peak[-n], l[0])
            to_delete.extend(delete_peak)
        for i in range(len(l[:-1])):
            delete_peak = [peak for peak in peaks if l[i] < peak < l[i+1]]
            if (m + n) >= len(delete_peak):
                delete_peak = delete_peak
                to_delete.extend(delete_peak) 
                atempausen = to_zeros(atempausen, l[i], l[i+1])
            else:
                delete_peak = delete_peak[:m] + delete_peak[-n:]
                atempausen = to_zeros(atempausen, l[i], delete_peak[m-1])
                atempausen = to_zeros(atempausen, delete_peak[-n], l[i+1])
                to_delete.extend(delete_peak)
        # Ende des Datensatzes, überprüfe, ob am Ende eine Atempause ist, wenn nicht, dann lösche m peaks rechts davon
        if atempausen[-1] == 1:
            delete_peak = [peak for peak in peaks if peak > l[-1]][:m]
            atempausen = to_zeros(atempausen, l[-1], delete_peak[m-1])
            to_delete.extend(delete_peak)
    peaks = [element for element in peaks if element not in to_delete]

    whole_time = [x_time for x_time in time_list if start_time <= x_time <= end_time]
    minuten_liste = [time for time in whole_time[::7500]]
    anzahl_pro_minute = []
    for i, minute in zip(range(len(minuten_liste)), minuten_liste):    
        # Zähle die Einträge aus der zweiten Liste, die in diese Minute fallen
        minute_atempause = atempausen[i*7500:i*7500 + 7499]
        nur_pause = [flag for flag in minute_atempause if flag == 1]
        l = len(nur_pause)
        anzahl_pro_minute.append(l)


    # Zeitpunkte der Peaks
    peak_zeiten = np.array(zeit)[peaks]  # Zeitpunkte der Peaks in Sekunden

    minuten_intervall = 60  # 60 Sekunden
    anzahl_minuten = int(np.ceil(spiro_resorted[0][-1] / minuten_intervall))  # Anzahl der Minuten im Datensatz
    atemfrequenz = []

    #Für Jede Minute wird die Atemfrequenz bestimmt
    for minute in range(anzahl_minuten):
        start = minute * minuten_intervall
        ende = (minute + 1) * minuten_intervall
        peaks_in_minute = np.sum((peak_zeiten >= start) & (peak_zeiten < ende))
        frequenz = peaks_in_minute  # Peaks pro Minute (Ein- und Ausatmung sind ein Peak)
        atemfrequenz.append(frequenz)



    #atemfrequenz = funktionen.n_amtungen_pro_minute(flussdaten, zeit)
    #Anpassung der Frequenz mit Berücksichtigung auf Atempausen

    atemfrequenz = [(a / b) * 7500 for a, b in zip(atemfrequenz, anzahl_pro_minute)]
    atemfrequenz[0] = atemfrequenz[0]
    # Zeitpunkte für die Minuten
    minuten_zeiten = np.arange(anzahl_minuten)
    # Startzeit angeben
    startzeit = datetime(2025, 1, 7, 18, 10, 40)  # 7. Januar 2025, 12:00 Uhr
    # Minuten in Uhrzeiten umwandeln
    uhrzeiten = [startzeit + timedelta(minutes=int(m)) for m in minuten_zeiten]
    uhrzeiten_formatiert = [zeit.strftime("%H:%M:%S") for zeit in uhrzeiten]
    return atemfrequenz2, uhrzeiten_formatiert


# in welchem Format sollen die anfangszeiten und endezeiten sein? Uhrzeit, Zeit nach start oder Anzahl Punkte 
def get_data_phase(data, anfang_ruhe, ende_ruhe, startzeit):
    # Falls die Eingaben Strings sind, in datetime umwandeln
    if isinstance(anfang_ruhe, str):
        anfang_ruhe = datetime.strptime(anfang_ruhe, "%H:%M:%S")
        ende_ruhe = datetime.strptime(ende_ruhe, "%H:%M:%S")
        startzeit = datetime.strptime(startzeit, "%H:%M:%S")

    # Dauer berechnen (Differenz zur Startzeit)
    if isinstance(anfang_ruhe, datetime) and isinstance(ende_ruhe, datetime):
        anfang_diff = (anfang_ruhe - startzeit).total_seconds()
        ende_diff = (ende_ruhe - startzeit).total_seconds()
        
        # In Einheiten umrechnen (1 Einheit = 8ms = 0.008s)
        anfang_index = int(anfang_diff / 0.008)
        ende_index = int(ende_diff / 0.008)
        
        # Array kürzen
        data = data[anfang_index:ende_index]
    return data

import numpy as np


def describe_array(arr, kanal):
    """
    Berechnet statistische Kennzahlen für ein NumPy-Array oder eine Liste.
    
    :param arr: Liste oder NumPy-Array numerischer Werte
    :return: Dictionary mit statistischen Werten
    """
    arr = np.asarray(arr)  # Sicherstellen, dass es ein NumPy-Array ist

    if arr.size == 0:  # Falls das Array leer ist
        return {"error": "Das Array ist leer."}

    return {
        "count": arr.size,
        "mean": np.mean(arr),
        "std": np.std(arr, ddof=1),  # Standardabweichung (Stichprobe)
        "min": np.min(arr),
        "25%": np.percentile(arr, 25),
        "50% (Median)": np.median(arr),
        "75%": np.percentile(arr, 75),
        "max": np.max(arr),
        "Channel": kanal,
    }

def schreibe_ergebnisse(file, kategorie, ergebnisse):
    with open(file, "a") as f:  # "a" für Anhängen, damit neue Werte hinzugefügt werden
        f.write(f"{kategorie}\n")
        f.write("-" * len(kategorie) + "\n")  # Unterstreichung der Überschrift
        f.write(f"{ergebnisse}\n\n")  # Dictionary direkt schreiben
    return print("Ergebnisse geschrieben!")
