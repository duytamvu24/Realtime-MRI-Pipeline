# import libraries
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from datetime import datetime, timedelta
from korrelation_functions_exe import *


def read_ecg_file(filename):
    data, log_start_time = read_ecg(filename)
    return data, log_start_time

def read_spiro_file(filename_spiro):
    spiro_resorted = read_spiro_data(filename_spiro)
    return spiro_resorted


# ---------------- File selection window ----------------
def choose_files():
    def select_ecg():
        filename = filedialog.askopenfilename(title="Select ECG-File", filetypes=[("All Files", "*.*")])
        if filename:
            entry_ecg.delete(0, tk.END)
            entry_ecg.insert(0, filename)

    def select_spiro():
        filename = filedialog.askopenfilename(title="Select Spirofile", filetypes=[("All Files", "*.*")])
        if filename:
            entry_spiro.delete(0, tk.END)
            entry_spiro.insert(0, filename)

    def start_program():
        nonlocal root
        ecg_file = entry_ecg.get()
        spiro_file = entry_spiro.get()
        ecg_kanal = spin_ecg_channel.get()
        if not ecg_file or not spiro_file:
            tk.messagebox.showerror("Error", "Please select both files!")
            return
        root.destroy()
        launch_main_window(ecg_file, spiro_file, ecg_kanal)

    root = tk.Tk()
    root.title("Select Files")

    frm = ttk.Frame(root, padding=20)
    frm.pack(fill=tk.BOTH, expand=True)

    # ECG file
    ttk.Label(frm, text="ECG-File:").grid(row=0, column=0, sticky="w")
    entry_ecg = ttk.Entry(frm, width=50)
    entry_ecg.grid(row=0, column=1, padx=5)
    ttk.Button(frm, text="Select", command=select_ecg).grid(row=0, column=2)

    # ECG channel selection
    ttk.Label(frm, text="ECG-File:").grid(row=1, column=0, sticky="w")
    spin_ecg_channel = ttk.Spinbox(frm, from_=1, to=4, width=5)
    spin_ecg_channel.set(2)  # default value
    spin_ecg_channel.grid(row=1, column=1, sticky="w", padx=5)
    
    # Spiro file
    ttk.Label(frm, text="Spiro-File:").grid(row=2, column=0, sticky="w")
    entry_spiro = ttk.Entry(frm, width=50)
    entry_spiro.grid(row=2, column=1, padx=5)
    ttk.Button(frm, text="Select Spiro-File", command=select_spiro).grid(row=2, column=2)

    ttk.Button(frm, text="Start", command=start_program).grid(row=3, column=1, pady=10)

    root.mainloop()


# Main window
def launch_main_window(filename, filename_spiro, ecg_kanal):
    global root, signals, x_datas, clicked_points, time_stamp, current_signal_index, data, spiro_resorted, log_start_time

    # Load selected files
    data, log_start_time = read_ecg_file(filename)
    spiro_resorted = read_spiro_file(filename_spiro)
    
    # Different ECG channels can be selected
    signals = [
        data[4+int(ecg_kanal)::4],              # ECG data
        spiro_resorted[1]                       # Spiro data
    ]

    x_datas = [
        np.arange(len(signals[0])) * 0.0025,   # ECG timeline
        np.arange(len(signals[1])) * 0.008     # Spiro timeline
    ]

    # Clicked peak storage
    clicked_points = []
    time_stamp = []
    current_signal_index = 0

    # Tkinter plot
    def find_nearest_point(x_click, x_data, y_data):
        idx = (np.abs(x_data - x_click)).argmin()
        return idx, y_data[idx]

    def onclick(event):
        if point_mode.get():
            if event.xdata is not None and event.ydata is not None:
                x_data = x_datas[current_signal_index]
                y_data = signals[current_signal_index]
                nearest_x, nearest_y = find_nearest_point(event.xdata, x_data, y_data)
                time_stamp.clear()
                time_stamp.append((nearest_x, nearest_y))
                print(f"Selected point at signal {current_signal_index+1}: ({nearest_x}, {nearest_y:.2f})")
                update_plot()

    def update_plot(*args):
        ax.clear()
        x_data = x_datas[current_signal_index]
        y_data = signals[current_signal_index]
        ax.plot(x_data, y_data, label=f"Signal {current_signal_index+1}")

        if time_stamp:
            x_pts, y_pts = zip(*time_stamp)
            ax.scatter(x_data[list(x_pts)], y_pts, color="red", s=100, label="Point")

        xmin = slider_min.get()
        xmax = slider_max.get()
        ax.set_xlim([xmin, xmax])

        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [a.u.]")
        ax.legend()
        ax.grid(True)
        canvas.draw()

    def save_point_and_next():
        global current_signal_index
        if time_stamp:
            clicked_points.append(time_stamp.copy())
            print(f"Signal {current_signal_index+1}: Point saved {time_stamp[0]}")

            current_signal_index += 1
            if current_signal_index < len(signals):
                load_signal(current_signal_index)
            else:
                print("All signals processed!")
                print("Saved Points:", clicked_points)
                final_correlation()

    def load_signal(index):
        time_stamp.clear()
        x_data = x_datas[index]
        slider_min.config(to=int(x_data[-1]))
        slider_max.config(to=int(x_data[-1]))
        slider_min.set(0)
        slider_max.set(min(3000, int(x_data[-1])))
        update_plot()

    def final_correlation():
        print("\n--- Start temporal correlation ---")
    
        # Remove previous plot
        canvas.get_tk_widget().destroy()
    
        time_stamp_ecg = clicked_points[0][0]  # (index, value)
        time_stamp_spiro = clicked_points[1][0]
        
        
        sr_ecg = 0.0025
        sr_spiro = 0.008
    
        ecg_from_timestamp = data[6::4][time_stamp_ecg[0]:]
        spiro_from_timestamp = spiro_resorted[1][time_stamp_spiro[0]:]
    
        time_ecg = np.arange(0, len(ecg_from_timestamp) * sr_ecg, sr_ecg)[:len(ecg_from_timestamp)]
        time_spiro = np.arange(0, len(spiro_from_timestamp) * sr_spiro, sr_spiro)[:len(spiro_from_timestamp)]
    
        scaling_factor = 1500
        mean_ecg = np.mean(ecg_from_timestamp)
        mod_ecg_data = [((y - mean_ecg)/scaling_factor) + 1 for y in ecg_from_timestamp]
    
        # Frame for final plot
        result_frame = ttk.LabelFrame(root, text="Overlapping signals")
        result_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
    
        fig2, ax2 = plt.subplots(figsize=(10, 4))
        line_ecg, = ax2.plot(time_ecg, mod_ecg_data, label="ECG scaled")
        line_spiro, = ax2.plot(time_spiro, spiro_from_timestamp, label="Spiro")
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel("Amplitude [a.u.]")
        ax2.set_title("Overlapping Signals")
        ax2.legend()
        ax2.grid(True)
    
        canvas2 = FigureCanvasTkAgg(fig2, master=result_frame)
        canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas2.draw()
    
        # slider
        def update_final_xlim(event=None):
            xmin = slider_min.get()
            xmax = slider_max.get()
            ax2.set_xlim([xmin, xmax])
            canvas2.draw()
    
        # function to save correlated data
        def save_correlated_data():

            # calculate volume from flow
            spiro_vol = get_volume_from_flow(spiro_resorted[1])
            
            # indices of selected points
            time_stamp_ecg = clicked_points[0][0][0]
            time_stamp_spiro = clicked_points[1][0][0]
            indizes = [time_stamp_ecg, time_stamp_spiro]

            ecg_real_time_time_stamp = (2.5 / 1000) * time_stamp_ecg

            real_time_ecg1 = timeConverter(log_start_time)

            time_format = "%H:%M:%S"
            print(real_time_ecg1)

            real_time_ecg1 = real_time_ecg1[:8]
            time_object = datetime.strptime(real_time_ecg1, time_format)

            new_time = time_object + timedelta(seconds=ecg_real_time_time_stamp)
            
            np.save("timestamp.npy", new_time)
            np.save('indizes.npy', indizes)
            
            print("Data saved: ecg.npy, spiro.npy, indizes.npy")
    
        control_frame = ttk.Frame(result_frame)
        control_frame.pack(fill=tk.X, pady=5)
    
        slider_min = tk.Scale(control_frame, from_=0, to=time_ecg[-1], resolution=0.01,
                              orient="horizontal", label="X-Min", command=update_final_xlim, length=400)
        slider_min.set(0)
        slider_min.pack(side=tk.LEFT, padx=5)
    
        slider_max = tk.Scale(control_frame, from_=0, to=time_ecg[-1], resolution=0.01,
                              orient="horizontal", label="X-Max", command=update_final_xlim, length=400)
        slider_max.set(min(10, time_ecg[-1]))
        slider_max.pack(side=tk.LEFT, padx=5)
    
        save_data_btn = ttk.Button(result_frame, text="Save data", command=save_correlated_data)
        save_data_btn.pack(pady=5)

    # --- GUI ---
    root = tk.Tk()
    root.title("Temporal-Correlation")

    point_mode = tk.BooleanVar(value=False)

    fig, ax = plt.subplots(figsize=(10, 4))
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    fig.canvas.mpl_connect("button_press_event", onclick)

    frame = ttk.Frame(root)
    frame.pack(fill=tk.X, padx=10, pady=5)

    slider_min = tk.Scale(frame, from_=0, to=1000, orient=tk.HORIZONTAL, label="X-Min", command=update_plot)
    slider_min.pack(side=tk.LEFT, fill=tk.X, expand=True)

    slider_max = tk.Scale(frame, from_=0, to=1000, orient=tk.HORIZONTAL, label="X-Max", command=update_plot)
    slider_max.pack(side=tk.LEFT, fill=tk.X, expand=True)

    button_frame = ttk.Frame(root)
    button_frame.pack(pady=5)

    toggle_btn = ttk.Checkbutton(button_frame, text="Activate selecting points", variable=point_mode)
    toggle_btn.pack(side=tk.LEFT, padx=5)

    save_btn = ttk.Button(button_frame, text="Save points and next plot", command=save_point_and_next)
    save_btn.pack(side=tk.LEFT, padx=5)

    # Start with the first signal
    load_signal(0)

    root.mainloop()


# ---------------- Start ----------------
if __name__ == "__main__":
    choose_files()
