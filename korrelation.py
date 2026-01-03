import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from korrelation_functions import *

# Deine Datenladefunktionen (angenommen bereits importiert)
# from korrelation_functions import read_ecg, read_spiro_data

path = "D:/Uni/Master-Thesis/daten/28022025/daten/"
filename_ecg = f'{path}echtzeit-spiroergo-04-20250228-1721.ecg'
filename_spiro = f"{path}Echtzeit-spiroergo04_28_2_2025_18_19raw.log"

# Variablen für die Daten
data = None
y_data = None
x_data = None
log_start_time = None
spiro_resorted = None
# Liste der Datensätze (Pfad und Typ)
datasets = [
    {"type": "ecg", "file": filename_ecg},
    {"type": "spiro", "file": filename_spiro},
]

points_per_dataset = [[], []]  # Punkte für ECG und Spiro
current_dataset_idx = 0
clicks_enabled = False

# Parameter zur Signalbeschneidung (wie von dir)
perc = 0
perc_end = 0.15

def load_ecg_dataset(filepath, kanal = 2):
    global data, y_data, x_data, log_start_time
    data, log_start_time = read_ecg(filepath)
    kanal_dat = kanal + 4
    y = data[kanal_dat::4]
    start_idx = int(perc * len(y))
    end_idx = int(perc_end * len(y))
    y_data = y[start_idx:end_idx]
    x_data = np.arange(len(y_data)) * 0.0025

def load_spiro_dataset(filepath):
    global data, y_data, x_data, spiro_resorted
    spiro_resorted = read_spiro_data(filepath)
    # Beispiel: Hier musst du die x_data und y_data passend zu deinem Spiro-Datensatz setzen.
    y_data = spiro_resorted[1]
    x_data = np.arange(len(y_data)) * 0.008

def load_dataset(index):
    if datasets[index]["type"] == "ecg":
        load_ecg_dataset(datasets[index]["file"])
    elif datasets[index]["type"] == "spiro":
        load_spiro_dataset(datasets[index]["file"])
    else:
        raise ValueError("Unbekannter Datensatztyp")

def find_nearest_point(x_click):
    idx = (np.abs(x_data - x_click)).argmin()
    return idx, y_data[idx]

def on_plot_click(event):
    global clicks_enabled
    if not clicks_enabled:
        return

    if event.xdata is not None:
        idx, val = find_nearest_point(event.xdata)
        points_per_dataset[current_dataset_idx].append((idx, val))
        print(f"Punkt gesetzt bei ({x_data[idx]:.2f}, {val:.2f}) im Datensatz {current_dataset_idx+1}")
        update_plot()
        clicks_enabled = False
        print("Klick deaktiviert nach Punktsetzung")

def update_plot():
    ax.clear()
    ax.plot(x_data, y_data, label=f"Datensatz {current_dataset_idx+1}: {datasets[current_dataset_idx]['type'].upper()}")

    if points_per_dataset[current_dataset_idx]:
        x_pts, y_pts = zip(*points_per_dataset[current_dataset_idx])
        ax.scatter([x_data[i] for i in x_pts], y_pts, color="red", s=100, label="Manuelle Punkte")

    ax.set_xlabel("Zeit")
    ax.set_ylabel("Amplitude")
    ax.set_title(f"Manuelle Punktwahl - Datensatz {current_dataset_idx+1}")
    ax.legend()
    ax.grid(True)
    canvas.draw_idle()

def enable_clicks():
    global clicks_enabled
    clicks_enabled = True
    print("Klick aktiviert")

def remove_last_point():
    if points_per_dataset[current_dataset_idx]:
        removed = points_per_dataset[current_dataset_idx].pop()
        print(f"Letzten Punkt entfernt: ({x_data[removed[0]]:.2f}, {removed[1]:.2f}) im Datensatz {current_dataset_idx+1}")
        update_plot()

def switch_to_next_dataset():
    global current_dataset_idx, clicks_enabled
    if current_dataset_idx < len(datasets) - 1:
        current_dataset_idx += 1
        load_dataset(current_dataset_idx)
        update_plot()
        clicks_enabled = False
        print(f"Zu Datensatz {current_dataset_idx+1} gewechselt: {datasets[current_dataset_idx]['type'].upper()}")
        btn_save_points.pack_forget()
    else:
        print("Kein weiterer Datensatz.")

    check_save_button()

def check_save_button():
    if all(len(points) > 0 for points in points_per_dataset):
        btn_save_points.pack(side=tk.LEFT, padx=5, pady=5)
    else:
        btn_save_points.pack_forget()

def save_all_points():
    print("Alle Punkte speichern:")
    for i, pts in enumerate(points_per_dataset):
        print(f"Datensatz {i+1} ({datasets[i]['type']}):")
        for idx, val in pts:
            print(f"  ({x_data[idx]:.2f}, {val:.2f})")
    print("Speicherung abgeschlossen (hier kannst du Datei schreiben etc.)")

def korreliere_und_plotten():
    print(points_per_dataset)
    time_stamp_spiro = points_per_dataset[0][0][0]  # Index aus Spiro
    timestamp_ecg = points_per_dataset[1][0][0]     # Index aus ECG

    print(f"Spiro-Timestamp: {time_stamp_spiro}, EKG-Timestamp: {timestamp_ecg}")
    print(f"EKG Startzeit: {log_start_time}")
    
    real_time_ecg = timeConverter(log_start_time)
    print(f"Uhrzeit EKG Beginn: {real_time_ecg}")

    # ECG-Zeitstempel berechnen
    ecg_real_time_time_stamp = (2.5 / 1000) * timestamp_ecg
    print(f"Zeit nach Beginn bis EKG-Timestamp: {ecg_real_time_time_stamp:.3f} s")

    # Berechne tatsächliche Uhrzeit des EKG-Timestamps
    time_object = datetime.strptime(real_time_ecg, "%H:%M:%S")
    new_time = time_object + timedelta(seconds=ecg_real_time_time_stamp)
    print("Zeit EKG-Timestamp:", new_time.strftime("%H:%M:%S"))

    # Ausschnitte aus den Rohdaten
    sr_ecg = 0.0025
    sr_spiro = 0.008
    ecg_from_timestamp = data[6::4][timestamp_ecg:]
    spiro_from_timestamp = spiro_resorted[1][time_stamp_spiro:]

    # Zeitachsen
    time_ecg = np.arange(0, len(ecg_from_timestamp) * sr_ecg, sr_ecg)
    time_spiro = np.arange(0, len(spiro_from_timestamp) * sr_spiro, sr_spiro)

    # Längen angleichen
    min_len = min(len(time_ecg), len(time_spiro))
    time_ecg = time_ecg[:min_len]
    time_spiro = time_spiro[:min_len]
    ecg_from_timestamp = ecg_from_timestamp[:min_len]
    spiro_from_timestamp = spiro_from_timestamp[:min_len]

    # EKG skalieren
    scaling_factor = 1500
    mean_ecg = np.mean(ecg_from_timestamp)
    mod_ecg_data = [((y - mean_ecg)/scaling_factor) + 1 for y in ecg_from_timestamp]

    # Plot
    fig_corr, ax_corr = plt.subplots(figsize=(12, 4))
    ax_corr.plot(time_spiro, spiro_from_timestamp, label="Spiro", alpha=0.8)
    ax_corr.plot(time_ecg, mod_ecg_data, label="EKG (skaliert)", alpha=0.8)
    ax_corr.set_xlabel("Zeit in Sekunden")
    ax_corr.set_ylabel("Amplitude / Fluss")
    ax_corr.set_title("Korrelation von Spiro- und EKG-Signal")
    ax_corr.legend()
    ax_corr.grid(True)
    plt.tight_layout()
    plt.show()

# --- Tkinter GUI Setup ---
root = tk.Tk()
root.title("Datensatz Punkteditor")

fig, ax = plt.subplots(figsize=(10, 4))
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
toolbar.pack(side=tk.TOP, fill=tk.X)

canvas.mpl_connect("button_press_event", on_plot_click)

button_frame = tk.Frame(root)
button_frame.pack(fill=tk.X)

btn_activate = tk.Button(button_frame, text="Punkte setzen aktivieren", command=enable_clicks)
btn_activate.pack(side=tk.LEFT, padx=5, pady=5)

btn_remove = tk.Button(button_frame, text="Letzten Punkt entfernen", command=remove_last_point)
btn_remove.pack(side=tk.LEFT, padx=5, pady=5)

btn_next_dataset = tk.Button(button_frame, text="Zum nächsten Datensatz wechseln", command=switch_to_next_dataset)
btn_next_dataset.pack(side=tk.LEFT, padx=5, pady=5)

btn_save_points = tk.Button(button_frame, text="Alle Punkte speichern", command=save_all_points)
# Anfangs versteckt
btn_save_points.pack_forget()

btn_exit = tk.Button(button_frame, text="Beenden", command=root.destroy)
btn_exit.pack(side=tk.RIGHT, padx=5, pady=5)

btn_corr = tk.Button(button_frame, text="Signale korrelieren & plotten", command=korreliere_und_plotten)
btn_corr.pack(side=tk.LEFT, padx=5, pady=5)

# Erster Datensatz laden und plotten
load_dataset(current_dataset_idx)
update_plot()

from datetime import datetime, timedelta



root.mainloop()
