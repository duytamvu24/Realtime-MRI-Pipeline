import numpy as np
import matplotlib.pyplot as plt
import funktionen
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import json

class PeakRemoverApp:
    def __init__(self, root, time, data, peaks, kanal, show_atempause=False, save_fig=False):
        self.root = root
        self.time = time / 60
        self.data = data
        self.peaks = peaks
        self.removed_peaks = []
        self.kanal = kanal
        self.show_atempause = show_atempause
        self.save_fig = save_fig
        self.start_index = 0
        self.window_size = 10000
        self.allow_removal = False  # Steuerung für Peak-Entfernung
        self.allow_manual_addition = False # Steuerung für mannuelles Hinzufügen
        
        self.root.title("Interactive Peak-Removal")
        
        self.fig, self.ax = plt.subplots(figsize=(10, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.slider = tk.Scale(self.root, from_=0, to=len(self.time) - self.window_size, orient=tk.HORIZONTAL, 
                               label="Timescale", command=self.update_window)
        self.slider.pack(fill=tk.X)
        
        self.plot_data()
        
        self.canvas.mpl_connect("button_press_event", self.on_click)
        
        self.toggle_button = tk.Button(self.root, text="Remove Peaks", command=self.toggle_removal)
        self.toggle_button.pack()

        self.manual_add_button = tk.Button(self.root, text="Add peaks manually", command=self.toggle_manual_addition)
        self.manual_add_button.pack()

        self.restore_button = tk.Button(self.root, text="Restore last removed peak", command=self.restore_last_peak)
        self.restore_button.pack()
        
        self.save_button = tk.Button(self.root, text="Save changes", command=self.save_peaks)
        self.save_button.pack()
        
        self.show_button = tk.Button(self.root, text="Show saved peaks", command=self.show_saved_peaks)
        self.show_button.pack()
        
        if self.save_fig:
            self.save_full_plot()

        self.root.protocol("WM_DELETE_WINDOW", self.on_close_window)
    
    def plot_data(self):
        self.ax.clear()
        end_index = self.start_index + self.window_size
        time_window = self.time[self.start_index:end_index]
        data_window = self.data[self.start_index:end_index]

        peak_mask = (self.peaks >= self.start_index) & (self.peaks < end_index)
        peaks_window = self.peaks[peak_mask] - self.start_index

        if self.show_atempause:
            atempause = funktionen.bestimme_atempausen(self.data, 200)
            self.ax.plot(time_window, atempause[self.start_index:end_index], label="Atempausenflag", color="green")
        
        self.ax.plot(time_window, data_window, label="Signal", color="blue")
        self.ax.scatter(time_window[peaks_window], data_window[peaks_window], color="red", label="Peaks", picker=True)
        
        y_label = funktionen.kanal_infos.get(self.kanal, {}).get("y_label", "Messwert")
        self.ax.set_xlabel("Time in minutes")
        self.ax.set_ylabel(y_label)
        self.ax.legend()
        self.canvas.draw()
    
    def save_full_plot(self):
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(self.time, self.data, label="Signal [a.u.]", color="blue")
        ax.scatter(self.time[self.peaks], self.data[self.peaks], color="red", label="Peaks")
        if self.show_atempause:
            atempause = funktionen.bestimme_atempausen(self.data, 200)
            ax.plot(self.time, atempause, label="Atempausenflag", color="green")
        ax.set_xlabel("Time in minutes")
        ax.set_ylabel(funktionen.kanal_infos.get(self.kanal, {}).get("y_label", "Measurement value [a.u.]"))
        ax.legend()
        fig.savefig("full_plot.png")
        plt.close(fig)
        print("Vollständiger Plot als 'full_plot.png' gespeichert.")

    
    def update_window(self, value):
        self.start_index = int(value)
        self.plot_data()
    

    def on_click(self, event):
        if event.xdata is None or event.ydata is None:
            return  
        
        if self.allow_removal:
            distances = np.sqrt((self.time[self.peaks] - event.xdata) ** 2 + (self.data[self.peaks] - event.ydata) ** 2)
            closest_index = np.argmin(distances)
        
            self.removed_peaks.append(self.peaks[closest_index])
        
            self.peaks = np.delete(self.peaks, closest_index)
            self.plot_data()
        
        elif self.allow_manual_addition:
            distances = np.sqrt((self.time - event.xdata) ** 2 + (self.data - event.ydata) ** 2)
            new_peak_index = np.argmin(distances)
            if new_peak_index not in self.peaks:
                self.peaks = np.append(self.peaks, new_peak_index)
                self.peaks.sort()
                self.plot_data()


    def restore_last_peak(self):
        if self.removed_peaks:
            last_peak = self.removed_peaks.pop()
            self.peaks = np.append(self.peaks, last_peak)
            self.peaks.sort()  # Peaks sortieren, damit sie an der richtigen Stelle bleiben
            self.plot_data()
        else:
            print("Keine entfernten Peaks zum Wiederherstellen.")

    
    def toggle_removal(self):
        self.allow_removal = not self.allow_removal
        self.allow_manual_addition = False  # Deaktiviert das manuelle Hinzufügen
        self.toggle_button.config(text="Remove peaks (activated)" if self.allow_removal else "Remove peaks")
        self.manual_add_button.config(text="Add peaks manually")
    
    def toggle_manual_addition(self):
        self.allow_manual_addition = not self.allow_manual_addition
        self.allow_removal = False  # Deaktiviert die Peak-Entfernung
        self.manual_add_button.config(text="Add peaks manually (activated)" if self.allow_manual_addition else "Add peaks manually")
        self.toggle_button.config(text="Remove Peaks")
    
    def save_peaks(self):
        with open("updated_peaks.json", "w") as f:
            json.dump(self.peaks.tolist(), f)
        print("Geänderte Peaks gespeichert!")
    
    def show_saved_peaks(self):
        try:
            with open("updated_peaks.json", "r") as f:
                self.peaks = np.array(json.load(f))
            print("Gespeicherte Peaks geladen!")
            self.plot_data()
        except FileNotFoundError:
            print("Keine gespeicherten Peaks gefunden.")

    def on_close_window(self):
        """Wird aufgerufen, wenn das Fenster geschlossen wird."""
        print("Fenster wird geschlossen. Speichern oder Cleanup hier möglich.")
        self.root.destroy()  # Schließt das Tkinter-Fenster