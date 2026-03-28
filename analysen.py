import funktionen
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import filedialog
import json
from scipy.signal import find_peaks
from app import PeakRemoverApp
import BC_Spiro_Functions_flow_npy_general as sm

#plt.style.use(['science', 'notebook', 'grid'])

#Datenliste zum Erstellen der Excel

data_list = None
def run_analysis_1(spiro_resorted, start_time_ruhe, end_time_ruhe, start,  file_name, gewicht, size,on_close = None):
    global data_list  # Zugriff auf die globale Variable
    data_list = []  # Neue leere Liste initialisieren
    start_time_str = start_time_ruhe
    end_time_str = end_time_ruhe
    kanal = 1
    data = np.array(funktionen.get_data_phase(spiro_resorted[kanal], start_time_str, end_time_str,start))
    time = np.array(funktionen.get_data_phase(spiro_resorted[0], start_time_str, end_time_str,start))
    #### 2.
    flow = np.array(data)
    data_corr, insp_int, exp_int, faktor = sm.calc_correction(flow, 0,-1)
    vol_corr = np.cumsum(data_corr[:]) # komplett#
    #Schrittweise Integral
    #vol_intervall = vol_corr.copy()[:] #komplett
    #vol_uncorr = np.cumsum(flow[:])
    
    #Baseline Korrektur nach Halima
    vol_automatic_BC = sm.BC_vol(vol_corr,0, -1)
    #vol_automatic_BC_uncorr = BC_vol(vol_uncorr,0, -1)
    new_vol_automatic_BC = vol_automatic_BC * 8 / 1000
    #new_vol_automatic_BC = vol_corr
    peaks, _ = find_peaks(new_vol_automatic_BC, distance = 200, prominence = 0.3, height = 0.2)  # Maxima
    def main():
        data = np.array(new_vol_automatic_BC)
        time = np.array(funktionen.get_data_phase(spiro_resorted[0], start_time_str, end_time_str,start))
        # Tkinter-Fenster starten
        root = tk.Tk()
        app = PeakRemoverApp(root, time[:-1], data, peaks, kanal = 4, save_fig = True)
        def on_window_close():
            """Wird aufgerufen, wenn das Fenster geschlossen wird."""
            print("Erste Analyse abgeschlossen. Verarbeite Daten...")

            with open("updated_peaks.json", "r") as f:
                peaks = json.load(f)
            # Zeitintervall zwischen zwei Messungen (in Sekunden)
            delta_t = 0.008

            peak_zeiten = np.array(time)[peaks]  # Zeitpunkte der Peaks in Sekunden
            
            minuten_intervall = 60  # 60 Sekunden
            anzahl_minuten = int(np.ceil((time[-1] - time[0]) / minuten_intervall))  # Anzahl der Minuten im Datensatz
            print(anzahl_minuten)
            atemfrequenz = []
            #new Zeit, um Frequenzen für jede Minute auszurechnen
            new_zeit =  np.array([x - time[0] for x in peak_zeiten])
            #Für Jede Minute wird die Atemfrequenz bestimmt
            for minute in range(anzahl_minuten):
                start = minute * minuten_intervall
                ende = (minute + 1) * minuten_intervall
                peaks_in_minute = np.sum((new_zeit >= start) & (new_zeit < ende))
                frequenz = peaks_in_minute  # Peaks pro Minute (Ein- und Ausatmung sind ein Peak)
                atemfrequenz.append(frequenz)
        
            # Umrechnung der Indizes der Peaks in Sekunden
            peaks_times = np.array(peaks) * delta_t
            
            # Gruppierung der Peaks in Minuten
            peaks_by_minute = {}
            for peak_time, peak_index in zip(peaks_times, peaks):
                minute = int(peak_time // 60)  # Bestimme die Minute (0-basiert)
                if minute not in peaks_by_minute:
                    peaks_by_minute[minute] = []
                peaks_by_minute[minute].append(new_vol_automatic_BC[peak_index])  # Speichere Peak-Wert
                
            minuten_liste = [time for time in time[::7500]]
            # Berechnung der Mittelwerte der Peaks für jede Minute
            mean_peaks_by_minute = {minute: np.mean(peaks) for minute, peaks in peaks_by_minute.items()}
            mean_peaks_minute = list(mean_peaks_by_minute.values())
            ziel_laenge = len(minuten_liste)
            mean_peaks_minute.extend([0] * (ziel_laenge - len(mean_peaks_minute)))
            
            # Output
            #kategorie = "Atemfrequenz"
            #atemfrequenz_data = funktionen.describe_array(atemfrequenz, kategorie)
            #data_list.append(atemfrequenz_data)
            #funktionen.schreibe_ergebnisse(file_name, kategorie, atemfrequenz_data)
            
            kategorie = "4.1 mean tidal volume [L]:"
            peak_values = np.array(new_vol_automatic_BC)[peaks]
            volume_data = funktionen.describe_array(peak_values, kategorie)
            data_list.append(volume_data)
            #funktionen.schreibe_ergebnisse(file_name, kategorie, volume_data)

            kategorie = "4.2 Minutevolume per body surface L/qm^2:"
            min_vol1 = [a * b for a, b in zip(mean_peaks_minute, atemfrequenz)]

            print("gewicht = " + str(gewicht))
            weight = gewicht
            height = size
            print(weight)
            bsa = 0.0235 * (height ** 0.42246) * (weight ** 0.51456)
            
            for i in range(len(mean_peaks_minute)):
                min_vol1[i] = min_vol1[i] / bsa

            minutenvolumen_data = funktionen.describe_array(min_vol1, kategorie)
            data_list.append(minutenvolumen_data)
            #funktionen.schreibe_ergebnisse(file_name, kategorie, minutenvolumen_data)            
            print("Results saved.")
            root.destroy()  # Fenster schließen
    
            # Falls eine weitere Analyse folgen soll, diese jetzt starten
            if on_close:
                on_close()  

        root.protocol("WM_DELETE_WINDOW", on_window_close)
        root.mainloop()
    main()


def run_analysis_2(spiro_resorted,  start_time_str, end_time_str, start, file_name, on_close = None):
    # Daten nach Phase: O2
    global data_list  # Zugriff auf die globale Variable
    kanal = 2
    ruhe_spiro = funktionen.get_data_phase(spiro_resorted[kanal],  start_time_str, end_time_str,start)
    ruhe_zeit = funktionen.get_data_phase(spiro_resorted[0], start_time_str, end_time_str,start)
    data = np.array(ruhe_spiro)
    time = np.array(ruhe_zeit)
    # Peaks bestimmen
    
    peaks, _ = find_peaks(ruhe_spiro, **funktionen.find_peaks_params[kanal])  # Maxima
    def main():
        # Tkinter-Fenster starten
        root = tk.Tk()
        app = PeakRemoverApp(root, time, data, peaks, kanal)
        
    
        def on_window_close():
            """Wird aufgerufen, wenn das Fenster geschlossen wird."""
            print("Vierte Analyse abgeschlossen. Verarbeite Daten...")
    
            with open("updated_peaks.json", "r") as f:
                peaks = json.load(f)
            spiro_neg_resorted = [-value for value in ruhe_spiro]
            valleys, _ = find_peaks(spiro_neg_resorted, **funktionen.find_valleys_params[kanal])
            atemfrequenz = len(peaks) / ((len(time) * 0.008)/ 60 )
            print("Durchschnittliche Atemfrequenz: ", np.mean(atemfrequenz))
            peak_values = np.array(ruhe_spiro)[peaks]
            valley_values = np.array(ruhe_spiro)[valleys]
            # Berechnung der üblichen statistischen Werte
            
            # Output
            kategorie = "3.1 O2-Maxima in %:"
            o2_peak_data = funktionen.describe_array(peak_values, kategorie)
            data_list.append(o2_peak_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, o2_peak_data)
            kategorie = "3.2 O2-Minima in %"
            o2_valley_data =  funktionen.describe_array(valley_values, kategorie)
            data_list.append(o2_valley_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, o2_valley_data)
            kategorie = "breathing rate: O2 in 1/min"
            atemfrequenz_o2_data = funktionen.describe_array(atemfrequenz, kategorie)
            #data_list.append(atemfrequenz_o2_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, atemfrequenz_o2_data)
            print("Ergebnisse gespeichert.")
            
            root.destroy()
            if on_close:
                on_close()
        
        root.protocol("WM_DELETE_WINDOW", on_window_close)
        root.mainloop()
    main()

def run_analysis_3(spiro_resorted,  start_time_str, end_time_str, start, file_name, on_close = None):
    # Daten nach Phase: CO2
    global data_list  # Zugriff auf die globale Variable
    kanal = 3
    ruhe_spiro = funktionen.get_data_phase(spiro_resorted[kanal], start_time_str, end_time_str,start)
    ruhe_zeit = funktionen.get_data_phase(spiro_resorted[0], start_time_str, end_time_str,start)
    
    # Peaks bestimmen
    
    peaks, _ = find_peaks(ruhe_spiro, **funktionen.find_peaks_params[kanal])  # Maxima
    def main():
        data = np.array(ruhe_spiro)
        time = np.array(ruhe_zeit)
        # Tkinter-Fenster starten
        root = tk.Tk()
        app = PeakRemoverApp(root, time, data, peaks, kanal)
        def on_window_close():
            """Wird aufgerufen, wenn das Fenster geschlossen wird."""
            print("Vierte Analyse abgeschlossen. Verarbeite Daten...")
    
            with open("updated_peaks.json", "r") as f:
                peaks = json.load(f)
            
            spiro_neg_resorted = [-value for value in ruhe_spiro]
            valleys, _ = find_peaks(spiro_neg_resorted, **funktionen.find_valleys_params[kanal])
            atemfrequenz = len(peaks) / ((len(time) * 0.008)/ 60 )
            peak_values = np.array(ruhe_spiro)[peaks]
            valley_values = np.array(ruhe_spiro)[valleys]
            # Berechnung der üblichen statistischen Werte
            
            
            # Output
            kategorie = "3.3 CO2-Maxima in %:"
            co2_peak_data = funktionen.describe_array(peak_values, kategorie)
            data_list.append(co2_peak_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, co2_peak_data)
            kategorie = "3.4 CO2-Minima in %"
            co2_valley_data =  funktionen.describe_array(valley_values, kategorie)
            data_list.append(co2_valley_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, co2_valley_data)
            kategorie = "Atemfrequenz: CO2 in 1/min"
            atemfrequenz_co2_data = funktionen.describe_array(atemfrequenz, kategorie)
            #data_list.append(atemfrequenz_co2_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, atemfrequenz_co2_data)
            print("Ergebnisse gespeichert.")
            
            root.destroy()
    
            if on_close:
                on_close()
    
        root.protocol("WM_DELETE_WINDOW", on_window_close)
        root.mainloop()
    main()

def run_analysis_4(spiro_resorted,  start_time_str, end_time_str, start, file_name ,on_close = None):
    # Daten nach Phase: Fluss
    global data_list  # Zugriff auf die globale Variable
    kanal = 1
    ruhe_spiro = funktionen.get_data_phase(spiro_resorted[kanal], start_time_str, end_time_str,start)
    ruhe_zeit = funktionen.get_data_phase(spiro_resorted[0], start_time_str, end_time_str,start)
    
    # Peaks bestimmen
    
    peaks, _ = find_peaks(ruhe_spiro, **funktionen.find_peaks_params[kanal])  # Maxima
    def main():
        data = np.array(ruhe_spiro)
        time = np.array(ruhe_zeit)
        # Tkinter-Fenster starten
        root = tk.Tk()
        app = PeakRemoverApp(root, time, data, peaks, kanal)
    
        def on_window_close():
            """Wird aufgerufen, wenn das Fenster geschlossen wird."""
            with open("updated_peaks.json", "r") as f:
                peaks = json.load(f)
            
            spiro_neg_resorted = [-value for value in ruhe_spiro]
            valleys, _ = find_peaks(spiro_neg_resorted, **funktionen.find_valleys_params[kanal])
            atemfrequenz = len(peaks) / ((len(time) * 0.008)/ 60 )
            print("Durchschnittliche Atemfrequenz: ", np.mean(atemfrequenz))
            peak_values = np.array(ruhe_spiro)[peaks]
            valley_values = np.array(ruhe_spiro)[valleys]
            # Berechnung der üblichen statistischen Werte
            
            
            # Output
            
            #kategorie  = "Atemfrequenz in 1/min"
            #atemfrequenz_flow = funktionen.describe_array(atemfrequenz, kategorie)
            #funktionen.schreibe_ergebnisse(file_name, "Atemfrequenz aus Fluss in 1/min", atemfrequenz_flow)
            #data_list.append(atemfrequenz_flow)
            kategorie = "1. Flow in L/s:"
            flow_peak_data = funktionen.describe_array(peak_values, kategorie)
            data_list.append(flow_peak_data)
            funktionen.schreibe_ergebnisse(file_name, kategorie, flow_peak_data)
            print("Vierte Analyse abgeschlossen. Verarbeite Daten...")
            print("Ergebnisse gespeichert.")


            #Ergebnisse in Excel
            # DataFrame erstellen
            df = pd.DataFrame(data_list)
            df = df[['Channel'] + [col for col in df.columns if col != 'Channel']]
            df = df.round(2)
            # Differenz berechnen
            start_time = datetime.strptime(start_time_str, "%H:%M:%S")
            end_time = datetime.strptime(end_time_str, "%H:%M:%S")
            time_diff = end_time - start_time
            total_seconds = time_diff.total_seconds()
            # Minuten und Sekunden extrahieren
            m, s = divmod(int(total_seconds), 60)
            formatted_time = f"{m} min {s} sec"
            time_data = {"Channel":"Length of phase: ", "mean":formatted_time}
            atemfre_data = {"Channel":"Breathing rate 1/min: ", "mean":atemfrequenz}
            df = pd.concat([df, pd.DataFrame([time_data])], ignore_index=True)
            df = pd.concat([df, pd.DataFrame([atemfre_data])], ignore_index=True)
            df = df.rename(columns = {"count": "Number Datapoints"})
            # DataFrame anzeigen
            df.to_excel(file_name + ".xlsx", index=False)
            root.destroy()
    
            if on_close:
                on_close()
    
        root.protocol("WM_DELETE_WINDOW", on_window_close)
        
        root.mainloop()
    main()