# import helper functions
from DataIntoDicom import *
from finetuning_functions import *
from korrelation_functions_exe import *
# import libraries
import tkinter as tk
from tkinter import filedialog, ttk

def start_analysis(params, log_widget):
    # main function of tkinter window
    def log(msg):
        log_widget.insert(tk.END, msg + "\n")
        log_widget.see(tk.END)  # immer nach unten scrollen
        log_widget.update_idletasks()

    # readout all files
    log("Start analysis...")
    print(params)
    main_path = params["main_path"]
    print(main_path)
    if not main_path.endswith("/"):
        main_path = main_path + "/"
    dicom_folder = params["dicom_path"]
    if not dicom_folder.endswith("/"):
        dicom_folder = dicom_folder + "/"
    filename_spiro = params["spiro_file"]
    rr_intervall = params["rr_intervall"]
    mean_volume = params["mean_volume"]
    threshold_volume = params["threshold_volume"]
    timestamp_numpy = params["timestamp_npy"]
    indizes_numpy = params["indizes_npy"]
    manuelle_verschiebung = params["manuelle_verschiebung"]
    n_ecg_phasen = params["n_ecg_phasen"]
    divider_ruhe = int(rr_intervall / n_ecg_phasen)

    # create folder and read out all slices in dicom folder
    # as long as there is only one folder add this folder to the path
    path = dicom_folder
    while len(os.listdir(path))==1: 
        path += os.listdir(path)[0]+'/'
    log(f"Dicom-Pfad: {path}")
    log("Following slices were found in DICOM-folder:")
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
        
        log(folder)
        for file in os.listdir(current):
            if file=='.DS_Store':
                continue
            i += 1
            current_file = current+'/'+file
            scan[file] = pydicom.dcmread(current_file, force=True)
        log(f"Read files from {current}")
        scans[folder] = scan
    print("Read", i, "files")

    # count imported dicom files
    scans = {k: scans[k] for k in sorted(scans.keys())}

    # sort folders to not mix up the slice order
    total = 0
    for k in scans.keys():  
        print(k)
        total += len(scans[k])
        print(len(scans[k]))
    print('total:', total)
    
    from natsort import natsorted
    # sort dicom images within each slice
    for k in scans.keys():
        scans[k] = {kk: scans[k][kk] for kk in natsorted(scans[k].keys())}
    print("Slices and images are sorted!")

    from math import floor
    import datetime 
    def convert_dicom_time(dicomtime):
        """
        Convert timestamp from dicom
        """
        dicomtime = dicomtime.strip('b\'')
        dicomhours = int(str(dicomtime)[0:2])   
        dicomminutes = int(str(dicomtime)[2:4])
        dicomseconds = int(str(dicomtime)[4:6])
        dicommilliseconds = int(float(str(dicomtime)[6:])*1000)
        return datetime.datetime(2000, 1, 1, dicomhours, dicomminutes, dicomseconds, dicommilliseconds*1000)

    # read ecg from dicom and dicom tags
    log("Read out log-time:")
    keys_slice = list(scans.keys())
    keys_images_first_slice = list(scans[keys_slice[0]])
    current_file = scans[keys_slice[0]][keys_images_first_slice[0]]
    #ds = pydicom.dcmread(current_file)
    # inspect meta information of a random scan
    current_file.add_new((0x08, 0x04),'LO', 1234)
    current_file.get_item((0x08, 0x04)).value
    
    
    dicomtime = str(current_file.get_item((0x08, 0x32)).value)
    dicomtime = convert_dicom_time(dicomtime)
    # Millisekunden seit Mitternacht berechnen
    msnmn = (dicomtime.hour * 3600000 + dicomtime.minute * 60000 + dicomtime.second * 1000 + dicomtime.microsecond / 1000)
    
    log(f"Das entspricht {msnmn} ms nach Mitternacht!")
    # MRI STart in ms nach Mitternacht für Start MRI

    # calculate timestamps for the respiratory flow measurements 
    # respiratory flow was measured every 10ms 
    # time from timestamp to first image of all selected images
    time_stamp_ecg_file = np.load(timestamp_numpy, allow_pickle=True)
    indizes = np.load(indizes_numpy, allow_pickle = True)
    spiro_resorted = read_spiro_data(filename_spiro)
    time_stamp_ecg = time_stamp_ecg_file.item()

    time_stamp_ecg_time = datetime.timedelta(
        hours=time_stamp_ecg.hour, minutes=time_stamp_ecg.minute,
        seconds=time_stamp_ecg.second, microseconds=time_stamp_ecg.microsecond
    )
    first_image_time = datetime.timedelta(
        hours=dicomtime.hour, minutes=dicomtime.minute,
        seconds=dicomtime.second, microseconds=dicomtime.microsecond
    )
    # calculate time to cut of of spiro to correlate with images
    diff = first_image_time - time_stamp_ecg_time
    print(diff)              # timedelta-Objekt
    print(diff.total_seconds())  # Unterschied in Sekunden (inkl. ms)
    
    steps_timestamp_to_first_image = int(diff.total_seconds() * 1000) / 8
    
    timestamp_and_first_slice_steps = int(steps_timestamp_to_first_image + indizes[1])
    print(f"Steps mit Betracht auf Zeitstempel bis zum ersten Bild: {timestamp_and_first_slice_steps }")
    
    # cutoff of timestamp and manual shift from finetuning
    total_cutoff = timestamp_and_first_slice_steps + manuelle_verschiebung
    # Schneiden der Spirodatei:
    data_flow_cutoff =  np.array(spiro_resorted[1][total_cutoff:])
    data_vol_cutoff = get_volume_from_flow(data_flow_cutoff)


    # new timestamps to get the same timestaps for spiro and images
    SR_spiro = 8 # in ms
    timestamps = [dicomtime]
    while len(timestamps) < len(data_vol_cutoff):
        timestamps.append(timestamps[-1] + datetime.timedelta(microseconds=(SR_spiro*1000)))
        # timestamps der Bilder
    timestamps = np.array(timestamps)

    log("Use timestamp to temporally correlate Spirodata with image data:")
    # iterate through every image and add 0x19, 0x01 dicom tag
    k = []
    for i in scans:
        l = []
        for j in scans[i]:
            scans[i][j].add_new((0x19, 0x01),'LO', str())
            # volume
            # infos for plotting
        k += l

    # copy ecg information (0018-1060) and overwrite instance number (0020-0013)
    # this is required for further post-processing of the files using the Software Circle Cvi 42
    log(f"Overwrite instance number for sorting files!")
    count = 0
    for slice_name in scans:
        for image_name in scans[slice_name]:
            current_scan = scans[slice_name][image_name]
            ecg = str(current_scan.get_item((0x18, 0x1060)).value).strip('b\'')
            current_scan[0x20, 0x0013].value = floor(float(ecg)) # truncate decimal places
            count += 1
    log(f'{count} images processed')

    # sort all imported scans by their filename
    scans_sort = OrderedDict(sorted(scans.items()))
    list(scans_sort.keys())
    
    for scan in scans_sort:
        d = scans_sort[scan]
        scans_sort[scan] = OrderedDict(sorted(d.items()))
        
    # check if sorted
    for scan in scans_sort:
        l = list(scans_sort[scan].keys())
        for i in range(0, len(l)-1):
            assert l[i]<l[i+1]

    # start adding the spiro data into dicom tags
    log("Fill in spirometry values into Dicom-Tags:")
    k = []
    for i in scans_sort:
        l = []
        log(f"Reinschreiben in {i}")
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

    # binning by inspiration/expiration and ecg groups 
    # the group id is calculated as (time after r-peak + 16.5)//37
    log("Preparation for binning: Every image will be labeled into an ECG class and into inspiration or expiration.")
    count = 0
    scans_sort_by_ecg = {}
    for i in scans_sort:
        # iteration through slices
        scans_sort_by_ecg[i] = {}
        for j in scans_sort[i]:
            # iterate through images
            count += 1
            current_scan = scans_sort[i][j]
            ecg = float(str(current_scan.get_item((0x18, 0x1060)).value).strip('b\''))
            ecg_group = int((ecg+16.5)//divider_ruhe) # calculate group id
            if ecg_group == n_ecg_phasen:
                ecg_group = 0
            if ecg_group not in scans_sort_by_ecg[i].keys():
                scans_sort_by_ecg[i][ecg_group] = [[], []]
            if float(current_scan.get_item((0x08, 0x03)).value) >= 0:
                # inspiration, fügt in 
                scans_sort_by_ecg[i][ecg_group][0].append(current_scan)
            else:
                # exspiration
                scans_sort_by_ecg[i][ecg_group][1].append(current_scan)
    log(f"{count }images processed")
    resp = 1 # 0 for inspiration, 1 for expiration

    # count all scans in scans_sort
    count = 0
    for i in scans_sort:
        # slice i
        for j in scans_sort[i]:
            # scan j
            count += 1
    log(f"TEST: Anzahl Bilder: {count}")
    log("Aufteilen der Bilder in Endexspiration und jeweils ein Bild aus jeder EKG-Klasse, wo das Volumen unter dem angegebenen threshold ist!")
    log("Selecting image in desired respiratory phase. For each ECG-class one image will be selected with a volume closest to the threshold!")
    min_vol, max_vol = np.min(ecg_v_f.reshape(-1,3)[:,1]), np.max(ecg_v_f.reshape(-1,3)[:,1])
    
    scans_sort_by_ecg_and_vol = {}
    
    # initialize empty dictionary
    for i in scans_sort:
        for ecg_group in range(n_ecg_phasen):
            scans_sort_by_ecg[i] = {}
            scans_sort_by_ecg[i][ecg_group] = [{'min_th':[], 'over_th':[],},
                                              {'min_th':[], 'over_th':[]}]
    
    count = 0
    count_less_24 = 0
    # binning ecg group in exspiration
    for i in scans_sort:
        for j in scans_sort[i]:
            count += 1
            current_scan = scans_sort[i][j]
            ecg = float(str(current_scan.get_item((0x18, 0x1060)).value).strip('b\''))
            # +16.5ms, da trigger time mittig vom scan ist 
            ecg_group = int((ecg+16.5)//divider_ruhe)
            if ecg_group == n_ecg_phasen:
                ecg_group = 0
            if ecg_group <= n_ecg_phasen - 1:
                # Zählt wie viele in 24 oder darunter sind
                count_less_24 += 1
                # aber es gibt doch bereits die ecg gruppe
                if ecg_group not in scans_sort_by_ecg[i].keys():
                    scans_sort_by_ecg[i][ecg_group] = [{'min_th':[], 'over_th':[],},
                                                      {'min_th':[], 'over_th':[]}]
                    #print(f"extra_bin hinzugefügt {i}, {ecg_group}")
                if float(current_scan.get_item((0x08, 0x03)).value) <= 0:
                    # exspiration
                    vol = float(current_scan.get_item((0x08, 0x04)).value)
                    if vol<=(mean_volume * threshold_volume):
                        
                        scans_sort_by_ecg[i][ecg_group][1]['min_th'].append(current_scan)
                    elif vol < mean_volume and vol > (mean_volume * threshold_volume):
                        scans_sort_by_ecg[i][ecg_group][1]['over_th'].append(current_scan)
                        
    log(f'{count_less_24}images kept')
    log(f'{count} images processed')
    log("Exporting data: For each ECG class, images at end-expiration with the smallest volume are selected. Missing entries will be filled.")
    # only export scans that have a respiratory volume clostest to the center of the bin
    # binning of endexspiration depending on mean tidal volume 
    for slc in scans_sort_by_ecg:
        # iterates through slice
        for ecg_group in scans_sort_by_ecg[slc]:
            #iterates through ecg_grouü
            for vol_group in scans_sort_by_ecg[slc][ecg_group][1]:
                # iterates through volume group
                group = scans_sort_by_ecg[slc][ecg_group][1][vol_group].copy()
                if vol_group=='min_th':
                    mid = (min_vol + (mean_volume * threshold_volume))/2
                elif vol_group=='over_th':
                    mid = mean_volume/2
                else:
                    log('invalid group')
                
                # pic scan with volume closest to center
                volumes = [float(scan.get_item((0x08, 0x04)).value) for scan in group]
                distances = np.array(volumes)-mid
                
                if len(distances)>0:
                    idx = np.argmin(distances)
                    to_export = group[idx]
                else:
                    log(f'empty bin: {vol_group}, in Schicht {slc},für die EKG-Gruppe: {ecg_group}')
                    to_export = find_closest_from_other_bin(mid, scans_sort_by_ecg, 
                                                              vol_group, slc, ecg_group)
                    
                    if to_export is None:
                        log('cannot fill bin')
                        continue
            
                        
                path_to_save = f'{main_path}test/output_ruhe/'+vol_group+'/'+str(1)+'/'+str(slc)+'/'+str(ecg_group)+'/test/'
                ensure_dir(path_to_save)
                pydicom.filewriter.dcmwrite(path_to_save+'scan'+'.dcm', to_export)
    log("Analysis completed ✅")


def ask_user_inputs():
    root = tk.Tk()
    root.title("Spiro-indo-DICOM-Module")

    results = {}

    def choose_file(var, filetypes=(("All Files", "*.*"),)):
        path = filedialog.askopenfilename(filetypes=filetypes)
        if path:
            var.set(path)

    def choose_dir(var):
        path = filedialog.askdirectory()
        if path:
            var.set(path)

    # --- Variables ---
    main_path = tk.StringVar()
    dicom_path = tk.StringVar()
    spiro_file = tk.StringVar()
    spiro_npy = tk.StringVar()
    timestamp_npy = tk.StringVar()
    indizes_npy = tk.StringVar()

    rr_intervall = tk.StringVar(value="1000")
    mean_volume = tk.StringVar(value="0")
    threshold_volume = tk.StringVar(value="30")  # Prozent
    manuelle_verschiebung = tk.StringVar(value = "0") # in zeitschritten
    n_ecg_phasen = tk.StringVar(value = "25")
    # --- GUI Layout ---
    frame = ttk.Frame(root, padding=10)
    frame.pack(fill="both", expand=True)

    row = 0
    ttk.Label(frame, text="Main path:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=main_path, width=50).grid(row=row, column=1)
    ttk.Button(frame, text="Select", command=lambda: choose_dir(main_path)).grid(row=row, column=2)
    row += 1

    ttk.Label(frame, text="DICOM path:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=dicom_path, width=50).grid(row=row, column=1)
    ttk.Button(frame, text="Select", command=lambda: choose_dir(dicom_path)).grid(row=row, column=2)
    row += 1

    ttk.Label(frame, text="Spirofile:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=spiro_file, width=50).grid(row=row, column=1)
    ttk.Button(frame, text="Select", command=lambda: choose_file(spiro_file)).grid(row=row, column=2)
    row += 1

    ttk.Label(frame, text="timestamp.npy:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=timestamp_npy, width=50).grid(row=row, column=1)
    ttk.Button(frame, text="Select", command=lambda: choose_file(timestamp_npy, (("NumPy Datei", "*.npy"),))).grid(row=row, column=2)
    row += 1

    ttk.Label(frame, text="indizes.npy:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=indizes_npy, width=50).grid(row=row, column=1)
    ttk.Button(frame, text="Select", command=lambda: choose_file(indizes_npy, (("NumPy Datei", "*.npy"),))).grid(row=row, column=2)
    row += 1

    ttk.Label(frame, text="rr_intervall in ms:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=rr_intervall, width=20).grid(row=row, column=1, sticky="w")
    row += 1

    ttk.Label(frame, text="mean_volume in L:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=mean_volume, width=20).grid(row=row, column=1, sticky="w")
    row += 1

    ttk.Label(frame, text="threshold_volume (%):").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=threshold_volume, width=20).grid(row=row, column=1, sticky="w")
    row += 1

    ttk.Label(frame, text="manual time shift in timesteps:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=manuelle_verschiebung, width=20).grid(row=row, column=1, sticky="w")
    row += 1

    ttk.Label(frame, text="Number ECG-Phases:").grid(row=row, column=0, sticky="w")
    ttk.Entry(frame, textvariable=n_ecg_phasen, width=20).grid(row=row, column=1, sticky="w")
    row += 1

    # --- Textfeld für Logs ---
    log_box = tk.Text(frame, height=12, width=70, wrap="word")
    log_box.grid(row=row, column=0, columnspan=3, pady=10)
    row += 1

    def start():
        results.update({
            "main_path": main_path.get(),
            "dicom_path": dicom_path.get(),
            "spiro_file": spiro_file.get(),
            "timestamp_npy": timestamp_npy.get(),
            "indizes_npy": indizes_npy.get(),
            "rr_intervall": float(rr_intervall.get()),
            "mean_volume": float(mean_volume.get()),
            "threshold_volume": float(threshold_volume.get()) / 100.0,  # Prozent in Dezimal
            "manuelle_verschiebung": int(manuelle_verschiebung.get()),
            "n_ecg_phasen": int(n_ecg_phasen.get())
        })
        start_analysis(results, log_box)

    ttk.Button(frame, text="START Analysis", command=start).grid(row=row, column=0, columnspan=3, pady=10)

    root.mainloop()
    return results

if __name__ == "__main__":
    params = ask_user_inputs()
    print(params)