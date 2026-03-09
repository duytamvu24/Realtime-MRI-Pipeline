# BC_Spiro_Functions_flow_npy_general liegt hier zwei Ordnerstufen vor dieser Datei. Wenn die Dateien im gleichen Ordner liegen, können die Pünktchen
# und die Backslashes entfernt werden
import funktionen
import numpy as np
import tkinter as tk
import json
from app import PeakRemoverApp 
import analysen
from tkinter import filedialog
from scipy.signal import find_peaks

# Globale Variablen
datei = None
spiro_resorted = None
start_time_str = None
start_time_ruhe = None
end_time_ruhe = None
start_time_belastung = None
end_time_belastung = None
gewicht = None
size = None

def open_channel_selection():

    """Opens window to choose channel."""
    
    kanal_window = tk.Toplevel(root)
    kanal_window.title("Select channel")
analyse
    tk.Label(kanal_window, text="Please select a channel (1-4):").pack()
    legende_text = """
    Channel 1: FLow in L/s
    Channel 2: O2-Volumefraction in %
    Channel 3: CO2-Volumefraction in %
    Channel 4: Volume in L/s
    """

    legende_label = tk.Label(kanal_window, text=legende_text, justify="left", font=("Arial", 10))
    legende_label.pack(padx=10, pady=10)
    kanal_entry = tk.Entry(kanal_window)
    kanal_entry.pack()
    print(gewicht)
    def show_data():
        global spiro_resorted, datei, start_time_str, start_time_ruhe, end_time_ruhe, start_time_belastung, end_time_belastung, size, gewicht

        """Function to look into data and plots!"""
        try:
            kanal = int(kanal_entry.get())
            if kanal not in [1, 2, 3, 4]:
                raise ValueError
            
            # Spiro-Daten auslesen
            data = np.array(spiro_resorted[kanal])
            time = np.array(spiro_resorted[0])

            # Peaks bestimmen
            peaks, _ = find_peaks(data, **funktionen.find_peaks_params[kanal])

            # PeakRemoverApp starten
            peak_window = tk.Toplevel(kanal_window)
            peak_window.title(f"Kanal {kanal} - PeakRemover")
            app = PeakRemoverApp(peak_window, time, data, peaks, kanal, show_atempause=False)
            peak_window.mainloop()

        except ValueError:
            label_status_kanal.config(text="Invalid Channel! Select 1-4:", fg="red")

    
    btn_show_data = tk.Button(kanal_window, text="Look into data", command=show_data)
    btn_show_data.pack()

    label_status_kanal = tk.Label(kanal_window, text="", fg="black")
    label_status_kanal.pack()

    def run_belastung_analyses():
        # Erstellt eine Textdatei mit allen Ergebnissen, um diese darein zu schreiben
        file_name = "Results_Exercisephase"
        with open(file_name, "w") as file:
            file.write("Wrting Results\n")
            file.write("=================\n\n")
        global spiro_resorted, datei, start_time_str, start_time_ruhe, end_time_ruhe, start_time_belastung, end_time_belastung, size, gewicht
        """Starts Analysis going through all channels"""
        def start_second_analysis():
            print("Start 2nd. Analysis...")
            analysen.run_analysis_2(spiro_resorted, 
            start_time_belastung, 
            end_time_belastung, 
            start_time_str, file_name, 
            on_close = start_third_analysis)
    
        def start_third_analysis():
            print("Start 3rd. Analysis...")
            analysen.run_analysis_3(spiro_resorted, 
            start_time_belastung, 
            end_time_belastung, 
            start_time_str, file_name, 
            on_close = start_fourth_analysis)
    
        def start_fourth_analysis():
            print("Start 4th. Analysis...")
            analysen.run_analysis_4(spiro_resorted, 
            start_time_belastung, 
            end_time_belastung, 
            start_time_str, file_name, 
            on_close = None)
            print("All Analysis completed. Results of last analysis:")
    
        # Starte erste Analyse
        analysen.run_analysis_1(spiro_resorted, 
            start_time_belastung, 
            end_time_belastung, 
            start_time_str, file_name, 
            gewicht, size, 
            on_close = start_second_analysis)

    

    def run_ruhe_analyses():
        # Erstellt eine Textdatei mit allen Ergebnissen, um diese darein zu schreiben
        file_name = "Results_Restingphase"
        with open(file_name, "w") as file:
            file.write("Wrting Results\n")
            file.write("=================\n\n")
        global spiro_resorted, datei, start_time_str, start_time_ruhe, end_time_ruhe, start_time_belastung, end_time_belastung, gewicht, size
        """Starts Analysis going through all channels"""
        def start_second_analysis():
            print("Start 2nd Analysis...")
            analysen.run_analysis_2(spiro_resorted, 
            start_time_ruhe, 
            end_time_ruhe, 
            start_time_str, file_name, 
            on_close = start_third_analysis)
    
        def start_third_analysis():
            print("Start §TH Analysis...")
            analysen.run_analysis_3(spiro_resorted, 
            start_time_ruhe, 
            end_time_ruhe, 
            start_time_str, file_name, 
            on_close = start_fourth_analysis)
    
        def start_fourth_analysis():
            print("Start 4th. analysis...")
            analysen.run_analysis_4(spiro_resorted, 
            start_time_ruhe, 
            end_time_ruhe, 
            start_time_str, file_name, 
            on_close = None)
            print("All Analysis completed. Results of last analysis:")
    
        # Starte erste Analyse
        print(gewicht, size)
        analysen.run_analysis_1(spiro_resorted, 
            start_time_ruhe, 
            end_time_ruhe, 
            start_time_str, file_name, 
            gewicht, size, 
            on_close = start_second_analysis)

    btn_run_analysis = tk.Button(kanal_window, text="Start Rest-Analysis", command=run_ruhe_analyses)
    btn_run_analysis.pack()
    btn_run_analysis = tk.Button(kanal_window, text="Start Exercise-Analysis", command=run_belastung_analyses)
    btn_run_analysis.pack()
    

def start_app():
    global spiro_resorted, datei, start_time_str, start_time_ruhe, end_time_ruhe, start_time_belastung, end_time_belastung, gewicht, size

    datei = filedialog.askopenfilename(title="Select File")
    print(datei)
    spiro_resorted = funktionen.read_spiro_data(datei)
    start_time_str = start_time.get()
    start_time_ruhe = start_time_ruhe.get()
    end_time_ruhe = end_time_ruhe.get()
    start_time_belastung = start_time_belastung.get()
    end_time_belastung = end_time_belastung.get()
    gewicht = int(gewicht_entry.get())
    size = int(size_entry.get())

    #Volumen berechnen und in Kanal 4 speichern
    volume = funktionen.get_volume_from_flow(spiro_resorted[1])
    volume = np.append(volume, 4)
    spiro_resorted.append(volume)

    if not datei:
        label_status.config(text="Error: No file selected!", fg="red")
        return
    # Öffne Fenster zur Kanalauswahl
    open_channel_selection()


# GUI zur Parameter-Eingabe
root = tk.Tk()
root.title("Spiro-Module")

#Design des ersten Fensters:
tk.Label(root, text = "Please enter the Measurement-Times! Exp: 18:10:40").pack()


tk.Label(root, text="Start Measurement:").pack()
start_time = tk.Entry(root)
start_time.pack()
start_time.insert(0, "00:00:00")

tk.Label(root, text="Start Restingphase:").pack()
start_time_ruhe = tk.Entry(root)
start_time_ruhe.pack()

tk.Label(root, text="End Restingphase:").pack()
end_time_ruhe = tk.Entry(root)
end_time_ruhe.pack()

tk.Label(root, text="Start Exercisephase:").pack()
start_time_belastung = tk.Entry(root)
start_time_belastung.pack()

tk.Label(root, text="End Exercisephase:").pack()
end_time_belastung = tk.Entry(root)
end_time_belastung.pack()

tk.Label(root, text="Weight (kg):").pack()
gewicht_entry = tk.Entry(root)
gewicht_entry.pack()
print(gewicht_entry)
tk.Label(root, text="Height (cm):").pack()
size_entry = tk.Entry(root)
size_entry.pack()

btn_datei = tk.Button(root, text="Choose spirofile and start!", command=start_app)
btn_datei.pack()

label_status = tk.Label(root, text="", fg="black")
label_status.pack()

root.mainloop()
