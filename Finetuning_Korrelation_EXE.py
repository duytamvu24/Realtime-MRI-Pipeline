# Import Helperfunctions
from finetuning_functions import *
from korrelation_functions_exe import *
#Import libraries
import tkinter as tk
from tkinter import filedialog, messagebox
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Slider
import matplotlib.patches as patches
from matplotlib.path import Path

# --- GUI-Logic ---
params = {}
# Global variables
# paraters for saving data of roi
roi_points = []
roi_polygon = None
roi_plot = None
point_plots = []
# parameters to save dicom data
ArrayDicomsort = None
AT = None
vol = None

plot_initialized = False
shift_entry = None
shift_button = None
hover_label = None

def select_dicom_folder():
    # Select Folder where Sclies with Dicom Images are stored
    folder = filedialog.askdirectory(title="Wähle DICOM-Ordner")
    if folder:
        params["dicom_folder"] = folder
        dicom_label.config(text=f"Ausgewählt: {folder}")

def select_spiro():
    # Select Spiro-log file
    file = filedialog.askopenfilename(
        title="Wähle Spiro-Datei mit raw log am Ende",
        filetypes=[("Log-Dateien", "*.log")]
    )
    if file:
        params["spiro_file"] = file
        spiro_label.config(text=f"Ausgewählt: {file}")

def select_indizes():
    # indices from the first korrelation as npy file
    file = filedialog.askopenfilename(
        title="Wähle Indize-Datei",
        filetypes=[("NumPy-Dateien", "*.npy")]
    )
    if file:
        params["indizes_file"] = file
        indizes_label.config(text=f"Ausgewählt: {file}")

def select_zeitstempel():
    # select timestamps 
    file = filedialog.askopenfilename(
        title="Wähle Zeitstempel-Datei",
        filetypes=[("NumPy-Dateien", "*.npy")]
    )
    if file:
        params["zeitstempel_file"] = file
        zeitstempel_label.config(text=f"Ausgewählt: {file}")

def select_bauchgurt():
    # select respiratory bellow file
    file = filedialog.askopenfilename(
        title="Wähle Bauchgurt-Datei",
        filetypes=[("Resp-Dateien", "*.resp")]
    )
    if file:
        params["bauchgurt_file"] = file
        bauchgurt_label.config(text=f"Ausgewählt: {file}")


def start_analysis():
    # function to start the analysis after every file has been selected
    global ArrayDicomsort, img, canvas, fig, ax, AT
    if "dicom_folder" not in params or "bauchgurt_file" not in params:
        messagebox.showerror("Fehler", "Bitte DICOM-Ordner und Bauchgurtdatei auswählen!")
        return
    
    try:
        selected_slice = int(slice_entry.get())
    except ValueError:
        messagebox.showerror("Fehler", "Bitte eine gültige Slice-Zahl eingeben!")
        return

    # Read out the dicom data from the dicom folder
    parent_dir = params["dicom_folder"]
    slice_liste = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    path = os.path.join(parent_dir, slice_liste[selected_slice]) + "/"
    sortedDicom,sortedArr,AT,amountFrames = read_dicomDir(path)
    # Store dicom data
    ArrayDicomsort = sortedDicom[:, :, AT[:,1].argsort()]

    slice_init = 100
    max_slice = ArrayDicomsort.shape[2] - 1

    # --- Figure für Tkinter erstellen ---
    fig, ax = plt.subplots()
    img = ax.imshow(ArrayDicomsort[:, :, slice_init], cmap='grey')

    # Klick-Handler für ROI
    def onclick(event):
        # function to enable clicking in the images to create a roi
        global roi_points, roi_polygon, roi_plot, point_plots
        if event.inaxes == ax and collecting_roi:
            roi_points.append((event.xdata, event.ydata))
            p, = ax.plot(event.xdata, event.ydata, "ro")  # Punkte speichern
            point_plots.append(p)
            fig.canvas.draw()

            if len(roi_points) == 4:
                # Altes ROI-Polygon löschen
                print(roi_points, roi_plot)

                # Polygon zeichnen
                xs, ys = zip(*roi_points)
                roi_plot, = ax.plot(list(xs) + [xs[0]], list(ys) + [ys[0]], "r-")
                fig.canvas.draw()

                # ROI als Path speichern
                roi_polygon = Path(roi_points)

                # Automatisch Kurve berechnen
                plot_signal_curve()

    fig.canvas.mpl_connect("button_press_event", onclick)

    # show the plots in the tkinter window
    for widget in plot_frame.winfo_children():
        widget.destroy()
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Slider in Tkinter
    def change_slice(val):
        # slider to watch every image of one slice
        idx = int(slice_slider.get())
        img.set_data(ArrayDicomsort[:, :, idx])
        canvas.draw_idle()
        
    slice_slider = tk.Scale(plot_frame, from_=0, to=max_slice, orient=tk.HORIZONTAL,
                            label="Slice auswählen", command=change_slice, length=400)
    slice_slider.set(slice_init)
    slice_slider.pack(pady=5)

    # Click to create the ROI
    btn_roi = tk.Button(plot_frame, text="Setze ROI", command=reset_roi)
    btn_roi.pack(pady=5)

    btn_finetuning = tk.Button(plot_frame, text = "Speichere Signalkurve!", command = finetuning)
    btn_finetuning.pack(pady=5)

def finetuning():
    # Finetuning: Correlation of Spirodata with Signalintensivity curve of one slice
    import numpy as np
    global AT, vol, curve_frame, hover_label
    # Readout data
    filename = params["bauchgurt_file"]
    sig_rb, real_time_bellow = read_bellow(filename)
    indizes_file = params["indizes_file"]
    indizes = np.load(indizes_file)
    timestamps_file = np.load(params["zeitstempel_file"], allow_pickle=True)
    timestamps = timestamps_file.item() 
    print(timestamps)
    spiro_file = params["spiro_file"]
    spiro_resorted = read_spiro_data(spiro_file)

    # time from timestamp to first image of slice
    from datetime import datetime, timedelta
    time_first_slice = str(AT[0][1])
    time_first_slice = float(time_first_slice)
    print(time_first_slice)
    
    hms = int(time_first_slice)  
    ms = round((time_first_slice - hms) * 1000)  # weil 4 Nachkommastellen
    hours = hms // 10000
    minutes = (hms % 10000) // 100
    seconds = hms % 100
    
    
    time_first_slice = datetime(1900, 1, 1, hours, minutes, seconds) + timedelta(milliseconds=ms)
    diff_ms_timestamp_first_slice = (time_first_slice - timestamps).total_seconds() * 1000
    diff_ms_timestamps = int(diff_ms_timestamp_first_slice/8)
    print(
        f"Differenz vom Zeitstempel bis zum ersten Bild: {diff_ms_timestamp_first_slice:.2f} ms und {(diff_ms_timestamp_first_slice/60000):.2f} min \n"
        f"Startzeit Timestamp: {timestamps}\n"
        f"Uhrzeit vom ersten Bild: {time_first_slice}"
    )
    
    # determine duration of slice in time and timesteps
    time_last_image = str(AT[-1][1])
    time_last_image = float(time_last_image)
    hms = int(time_last_image)  
    ms = round((time_last_image - hms) * 1000)  # weil 4 Nachkommastellen
    # Zerlegen in h, m, s
    hours = hms // 10000
    minutes = (hms % 10000) // 100
    seconds = hms % 100
    # new timeline to adjust the different acquisition times
    time_last_image = datetime(1900, 1, 1, hours, minutes, seconds) + timedelta(milliseconds=ms)
    duration_slice = int((time_last_image - time_first_slice).total_seconds() * 1000)
    
    print(f"Dauer des Slices: {duration_slice}")

    # Cut spiro data
    shift = 0
    if shift_entry is not None:
        val = shift_entry.get()
        print("val in finetuning: " + str(val))
        if val.strip() != "":
            try:
                shift = float(val)
            except ValueError:
                print("Ungültiger Wert, bitte Zahl eingeben.")
                shift = 0  
    

    # cut bellow data
    time_bellow_start = datetime.strptime(real_time_bellow, "%H:%M:%S.%f").replace(year=1900, month=1, day=1)
    # offset, manual offset to add onto the timestamp shift if a shift is needed for the finetuning
    offset_finetuning = int(shift) # In timestamps
    duration_slice_steps = int(duration_slice/8)
    spiro_flow = spiro_resorted[1][(indizes[1] + diff_ms_timestamps + offset_finetuning) : (indizes[1] + diff_ms_timestamps + offset_finetuning) + duration_slice_steps]
    spiro_vol = get_volume_from_flow(spiro_flow)
    print(
        f"Schritten vom Timestamp bis slice: {diff_ms_timestamps} und in min {(diff_ms_timestamps*8/60000):.2f}\n"
        f"  Schritte durch Timestamp: {indizes[1]}\n"
        f"  Insgesamtes Schneiden:  {(indizes[1] + diff_ms_timestamps) * 8 / 60000:.2f}"
    )
    diff_bellow_start_to_first_image = (time_first_slice - time_bellow_start).total_seconds() * 1000
    diff_bellow_start_to_first_image_steps = int(diff_bellow_start_to_first_image/2.5)
    duration_slice_bellow_steps = int(duration_slice/2.5)
    bellow_timestamp_and_to_first_slice = diff_ms_timestamp_first_slice + indizes[0] * 2.5
    bellow_timestamp_and_to_first_slice_steps = int(bellow_timestamp_and_to_first_slice / 2.5)
    # cut all the data
    sig_rb_cut = sig_rb[diff_bellow_start_to_first_image_steps : diff_bellow_start_to_first_image_steps + duration_slice_bellow_steps ]
    sig_rb_cut = sig_rb_cut / np.mean(sig_rb_cut) 
    print("sig_rb_cut:" + str(len(sig_rb_cut)))
    # Zeitskalen mit den verschiedenen Kanälen erstellen
    import numpy as np
    # Zeitskala der Bilder
    # zweite Spalte extrahieren
    times = AT[:, 1]
    
    
    from datetime import time
    
    
    def to_seconds(val):
        hms = int(val)
        frac = val - hms
        micro = int(round(frac * 1_000_000))
        
        hours = hms // 10000
        minutes = (hms % 10000) // 100
        seconds = hms % 100
        
        return hours*3600 + minutes*60 + seconds + micro/1_000_000
    
    seconds_list = [to_seconds(t) for t in times]
    
    # new timeline to show cut data of spiro, bellow and signalintensivity curve
    diffs = [s - seconds_list[0] for s in seconds_list]
    timescale_images = [diff*1000 for diff in diffs]
    # Zeitskala bellow
    import numpy as np
    scale_bellow = np.arange(0, int(len(sig_rb_cut) * 2.5), 2.5)
    scale_flow = np.arange(0, int(len(spiro_flow) * 8), 8)
    scale_spiro = np.arange(0, len(spiro_vol) * 8, 8) 
    vol = (vol / np.mean(vol))

    # Plot curve
    for widget in curve_frame.winfo_children():
            widget.destroy()
    # neue Kurve erstellen
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(scale_spiro,spiro_vol, label='Volumen aus der Volumetrie', linewidth=2)
    ax.plot(scale_bellow,sig_rb_cut, label='Bauchkurve', alpha=0.7)
    ax.plot(timescale_images,vol, label='Volumetrie aus der Signalkurve',)

    # Layout
    ax.set_title("Finetuning")
    ax.set_xlabel("Zeit in ms")
    ax.set_ylabel("Signalwert")
    ax.legend()
    ax.grid(True)
    if hover_label is None:
        hover_label = tk.Label(root, text="x: -", font=("Arial", 10))
        hover_label.pack(pady=8)
    
    hover_line = ax.axvline(x=0, color="gray", linestyle="--", linewidth=1, alpha=0.6)

    def on_hover(event):
        if event.inaxes == ax and event.xdata is not None:
            x_val = event.xdata
            hover_label.config(text=f"x: {x_val:.1f}")
            
            # x-Wert muss als Sequenz übergeben werden:
            hover_line.set_xdata([x_val])
            fig.canvas.draw_idle()
        else:
            hover_label.config(text="x: -")

    fig.canvas.mpl_connect("motion_notify_event", on_hover)
    # Canvas in Tkinter einfügen
    canvas = FigureCanvasTkAgg(fig, master=curve_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

def plot_everything():
    # function to correlate the whole data after reading out all the data from files
    import numpy as np
    global AT, vol, curve_frame, shift_entry, plot_initialized, hover_label
    # read out data

    parent_dir = params["dicom_folder"]
    filename = params["bauchgurt_file"]
    sig_rb, real_time_bellow = read_bellow(filename)
    print("realtimebellow:" +str(real_time_bellow))
    indizes_file = params["indizes_file"]
    indizes = np.load(indizes_file)
    timestamps_file = np.load(params["zeitstempel_file"], allow_pickle=True)
    timestamps = timestamps_file.item() 
    spiro_file = params["spiro_file"]
    spiro_resorted = read_spiro_data(spiro_file)
    slice_liste = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    slice_liste = sorted(slice_liste)
    selected_slice = 0
    path = os.path.join(parent_dir, slice_liste[selected_slice]) + "/"
    
    sortedDicom,sortedArr,AT,amountFrames = read_dicomDir(path)
    # time from timestamp to first image of slice
    from datetime import datetime, timedelta
    time_first_slice = str(AT[0][1])
    time_first_slice = float(time_first_slice)
    print(time_first_slice)
    
    hms = int(time_first_slice)  
    ms = round((time_first_slice - hms) * 1000)  # weil 4 Nachkommastellen
    hours = hms // 10000
    minutes = (hms % 10000) // 100
    seconds = hms % 100

    # manual shift
    shift = 0
    if shift_entry is not None:
        val = shift_entry.get()
        print("val: " + str(val)) 
        if val.strip() != "":
            try:
                shift = float(val)
            except ValueError:
                print("Ungültiger Wert, bitte Zahl eingeben.")
                shift = 0
    
    time_first_slice = datetime(1900, 1, 1, hours, minutes, seconds) + timedelta(milliseconds=ms)
    diff_ms_timestamp_first_slice = (time_first_slice - timestamps).total_seconds() * 1000
    print(
        f"Differenz vom Zeitstempel bis zum ersten Bild: {diff_ms_timestamp_first_slice:.2f} ms und {(diff_ms_timestamp_first_slice/60000):.2f} min \n"
        f"Startzeit Timestamp: {timestamps}\n"
        f"Uhrzeit vom ersten Bild: {time_first_slice}"
    )
    

    time_last_image = str(AT[-1][1])
    time_last_image = float(time_last_image)
    hms = int(time_last_image)  
    ms = round((time_last_image - hms) * 1000)  # weil 4 Nachkommastellen
    # Zerlegen in h, m, s
    hours = hms // 10000
    minutes = (hms % 10000) // 100
    seconds = hms % 100
    time_last_image = datetime(1900, 1, 1, hours, minutes, seconds) + timedelta(milliseconds=ms)
    duration_slice = int((time_last_image - time_first_slice).total_seconds() * 1000)
    
    print(f"Dauer des Slices: {duration_slice}")

    offset_finetuning = int(shift)
    print("offset_ ist:" + str(offset_finetuning))
    diff_ms_timestamps = int(diff_ms_timestamp_first_slice/8)
    
    spiro_flow = spiro_resorted[1][(indizes[1] + diff_ms_timestamps + offset_finetuning) :]
    time_bellow_start = datetime.strptime(real_time_bellow, "%H:%M:%S.%f").replace(year=1900, month=1, day=1)
    diff_bellow_start_to_first_image = (time_first_slice - time_bellow_start).total_seconds() * 1000
    diff_bellow_start_to_first_image_steps = int(diff_bellow_start_to_first_image/2.5)
    
    bellow_timestamp_and_to_first_slice = diff_ms_timestamp_first_slice + indizes[0] * 2.5
    bellow_timestamp_and_to_first_slice_steps = int(bellow_timestamp_and_to_first_slice / 2.5)
    sig_rb_cut = sig_rb[diff_bellow_start_to_first_image_steps : ]
    sig_rb_cut = sig_rb_cut / np.mean(sig_rb_cut) 
    
    spiro_vol = get_volume_from_flow(spiro_flow)
    import numpy as np
    # Zeitskala der Bilder
    # zweite Spalte extrahieren
    times = AT[:, 1]
    
    
    from datetime import time
    
    
    def to_seconds(val):
        hms = int(val)
        frac = val - hms
        micro = int(round(frac * 1_000_000))
        
        hours = hms // 10000
        minutes = (hms % 10000) // 100
        seconds = hms % 100
        
        return hours*3600 + minutes*60 + seconds + micro/1_000_000
    # recalculte from seconds to time
    seconds_list = [to_seconds(t) for t in times]
    
    # new timeline for plotting
    diffs = [s - seconds_list[0] for s in seconds_list]
    timescale_images = [diff*1000 for diff in diffs]
    # Zeitskala bellow
    import numpy as np
    scale_bellow = np.arange(0, int(len(sig_rb_cut) * 2.5), 2.5)
    scale_flow = np.arange(0, int(len(spiro_flow) * 8), 8)
    scale_spiro = np.arange(0, len(spiro_vol) * 8, 8) 
    #vol = (vol / np.mean(vol) +1)
    print(scale_spiro[-1], scale_bellow[-1], timescale_images[-1])
    
    
    # Plot curve
    from matplotlib.figure import Figure
    for widget in curve_frame.winfo_children():
            widget.destroy()
    # --- Neue Figure erstellen ---
    fig = Figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    line1, = ax.plot(scale_spiro, spiro_vol, label='Volumen aus Fluss', linewidth=2)
    line2, = ax.plot(scale_bellow, sig_rb_cut, label='Bauchkurve1')

    ax.set_title("Finetuning")
    ax.set_xlabel("Zeit in ms")
    ax.set_ylabel("Signalwert")
    ax.legend()
    ax.grid(True)


    # --- Hover-function ---
    if hover_label is None:
        hover_label = tk.Label(root, text="x: -", font=("Arial", 10))
        hover_label.pack(pady=8)

    def on_hover(event):
        if event.inaxes == ax and event.xdata is not None:
            x_val = event.xdata
            hover_label.config(text=f"x: {x_val:.1f}")
            
            # x-Wert muss als Sequenz übergeben werden:
            hover_line.set_xdata([x_val])
            fig.canvas.draw_idle()
        else:
            hover_label.config(text="x: -")

    fig.canvas.mpl_connect("motion_notify_event", on_hover)
    # --- Startbereich X-Achse ---
    window = 40000
    ax.set_xlim(0, window)

    # --- Canvas in Tkinter einbetten ---
    canvas = FigureCanvasTkAgg(fig, master=curve_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

    # --- Scroll-Funktion ---
    def update(val):
        pos = slider.get()
        ax.set_xlim(pos, pos + window)
        canvas.draw_idle()

    # --- Slider hinzufügen ---
    slider = tk.Scale(
        slider_frame,
        from_=0,
        to=max(scale_spiro) - window,
        orient="horizontal",
        label="Scroll durch Zeitachse",
        length=800,
        command=lambda v: update(v)
    )
    slider.pack()

    if not plot_initialized:
            shift_entry = tk.Entry(root)
            shift_entry.pack(pady=6)
            shift_entry.insert(0, "0")  # optionaler Startwert
            plot_initialized = True
    
collecting_roi = False
def reset_roi():
    """resets ROI"""
    global roi_points, roi_polygon, roi_plot, point_plots, collecting_roi

    # Reset
    roi_points = []
    roi_polygon = None
    collecting_roi = True

    # Alte Punkte löschen
    if roi_plot is not None:
        roi_plot.remove()
    for p in point_plots:
        p.remove()
    point_plots = []

    # Alte Kurve löschen
    for widget in curve_frame.winfo_children():
        widget.destroy()

    if canvas:
        canvas.draw_idle()

    print("ROI zurückgesetzt – bitte 4 Punkte setzen!")


def plot_signal_curve():
    global roi_polygon, ArrayDicomsort, collecting_roi, vol, hover_label
    if ArrayDicomsort is None:
        messagebox.showerror("Fehler", "Bitte zuerst Analyse starten!")
        return
    if roi_polygon is None:
        return  # Noch keine ROI gesetzt

    # Maske für ROI erstellen
    ny, nx = ArrayDicomsort.shape[0], ArrayDicomsort.shape[1]
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    coords = np.vstack((X.flatten(), Y.flatten())).T
    mask = roi_polygon.contains_points(coords).reshape((ny, nx))

    # Mittelwert pro Frame berechnen
    vol = [np.mean(ArrayDicomsort[:, :, i][mask]) for i in range(ArrayDicomsort.shape[2])]

    # --- alte Kurve löschen ---
    for widget in curve_frame.winfo_children():
        widget.destroy()

    # neue Kurve plotten
    fig2, ax2 = plt.subplots()
    ax2.plot(vol, color="blue")
    ax2.set_title("Signalintensitätskurve (ROI)")
    ax2.set_xlabel("Bildindex")
    ax2.set_ylabel("Durchschnittsintensität")

    canvas2 = FigureCanvasTkAgg(fig2, master=curve_frame)
    canvas2.draw()
    canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    collecting_roi = False  # Nach 4 Punkten fertig

# --- Tkinter Fenster ---
root = tk.Tk()
root.title("Daten-Loader mit ROI")
root.geometry("1200x700")

control_frame = tk.Frame(root)
control_frame.pack(side=tk.TOP, fill=tk.X, pady=10)

plot_frame = tk.Frame(root)
plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

curve_frame = tk.Frame(root)
curve_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

# Buttons und Eingaben
tk.Button(control_frame, text="Wähle DICOM-Ordner", command=select_dicom_folder).pack(pady=2)
dicom_label = tk.Label(control_frame, text="Noch kein Ordner gewählt")
dicom_label.pack()


tk.Label(control_frame, text="Slice Nummer:").pack()
slice_entry = tk.Entry(control_frame)
slice_entry.insert(0, "0")
slice_entry.pack()

tk.Button(control_frame, text="Wähle Bauchgurt-Datei", command=select_bauchgurt).pack(pady=2)
bauchgurt_label = tk.Label(control_frame, text="Noch keine Datei gewählt")
bauchgurt_label.pack()

tk.Button(control_frame, text="Wähle zeitstempel.npy", command=select_zeitstempel).pack(pady=2)
zeitstempel_label = tk.Label(control_frame, text="Noch keine Datei gewählt")
zeitstempel_label.pack()

tk.Button(control_frame, text="Wähle spiro.npy", command=select_spiro).pack(pady=2)
spiro_label = tk.Label(control_frame, text="Noch keine Datei gewählt")
spiro_label.pack()

tk.Button(control_frame, text="Wähle indizes.npy", command=select_indizes).pack(pady=2)
indizes_label = tk.Label(control_frame, text="Noch keine Datei gewählt")
indizes_label.pack()

tk.Button(control_frame, text="2. Betrachte einzelne Schichten!", command=start_analysis).pack(pady=5)

tk.Button(control_frame, text="1. Betrachte ganze Kurve!", command=plot_everything).pack(pady=5)

hover_label = tk.Label(root, text="x: -", font=("Arial", 10))
hover_label.pack(pady=9)

slider_frame = tk.Frame(root)
slider_frame.pack(pady=10)

root.mainloop()