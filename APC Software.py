#####################################################################
#                                                                   #
#  Developed by Francisco Rogé Vallone and Ivan Kaspierowicz        #
#  Contact: rogevallone@gmail.com                                    #
#  Date: 4/21                                                       #
#  Libraries needed: Matplotlib, Tkinter, Scipy, Numpy, Soundfile   #
#                                                                   #
#####################################################################



from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.ttk import Progressbar, Combobox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import soundfile as sf
from scipy import ndimage
from scipy import stats
import scipy.signal as sc
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def func_calculate(): # Main function used for parameters calculation. Callback of "Calculate" button 
    if len(RIRs) == 0:
        messagebox.showerror("Error", message="No data loaded")
        return
    if cmbx_bandfilter.current() == -1:
        messagebox.showerror("Error", message="No filter selected")
        return

    try:
        global results
        global signals
        global signalsdb
        global envelopes
        global averages
        global channels
        global fs
        
        changestate(0) # Disabled settings while calculating
        # Creation of vectors and band frequencies.
        channels = RIRchannels[:]
        strbands = bands_labels[cmbx_bandfilter.current()][
                   cmbx_minband.current():cmbx_minband.current() + cmbx_maxband.current() + 1]
        bands = bands_centers[cmbx_bandfilter.current()][
                cmbx_minband.current():cmbx_minband.current() + cmbx_maxband.current() + 1]
        # Creation of data allocation.
        results = np.zeros([len(parameters), len(bands) + 1, len(RIRs), 2])
        signals = []
        signalsdb = []
        envelopes = []
        fs = samplerate[:]
        # One-Third Octave Frequencies limits.
        if cmbx_bandfilter.current() == 0:
            freclim = (2 ** (1 / 2))
        elif cmbx_bandfilter.current() == 1:
            freclim = (2 ** (1 / 6))
        # Set progression bar
        total_steps = len(bands) * sum(channels)
        bar_step = 0
        pbar.grid()
        btn_calculate.grid_remove()
        
        # Signals 
        for s in range(len(RIRs)):
    
            list.append(signals, [])
            list.append(signalsdb, [])
            list.append(envelopes, [])
            # Channels
            for c in range(channels[s]):
    
                list.append(signals[s], [])
                list.append(signalsdb[s], [])
                list.append(envelopes[s], [])
                # Bands
                for b in range(len(bands) + 1):
                    # Inferior and superior frequencies and Window Lenght.
                    if b == 0: #global 
                        finf = (2 / fs[s]) * bands[0] / freclim
                        fsup = (2 / fs[s]) * bands[-1] * freclim
                        winlen = int(2000 * np.log10(fs[s] / (bands[0] / freclim)))
                    else: # band filter
                        finf = (2 / fs[s]) * bands[b - 1] / freclim
                        fsup = (2 / fs[s]) * bands[b - 1] * freclim
                        winlen = int(2000 * np.log10(fs[s] / (bands[b - 1] / freclim)))
                    if fsup > 1:
                        fsup = (fs[s] - 1) / fs[s]
                    # Filter bank order 8.
                    sos = sc.butter(8, [finf, fsup], btype="band", output="sos", analog=False)
                    # Filtering.
                    sgn = sc.sosfiltfilt(sos, RIRs[s][:, c])
                    list.append(signals[s][c], sgn)
                    # Filtered IR to dB.
                    list.append(signalsdb[s][c], np.abs(signals[s][c][b]))
                    # Determination of start sample as -20dB before the maximum  
                    IRstart = \
                    np.where(signalsdb[s][c][b][0:np.argmax(signalsdb[s][c][b])] <= np.max(signalsdb[s][c][b]) / 10)[0]
                    if len(IRstart) == 0:
                        IRstart = 0
                    else:
                        IRstart = IRstart[-1]
                    # Truncation from IR start
                    signals[s][c][b] = signals[s][c][b][IRstart:]
                    signalsdb[s][c][b] = signalsdb[s][c][b][IRstart:]
                    signalsdb[s][c][b][signalsdb[s][c][b] == 0] = 0.00000001
                    signalsdb[s][c][b] = 20 * np.log10(signalsdb[s][c][b])
                    # Truncation method
                    if cmbx_chunk.current() == 0:
                        [IRend, env] = lundeby(signalsdb[s][c][b], fs[s])
                    elif cmbx_chunk.current() == 1:
                        [IRend, env] = pepino(signalsdb[s][c][b], fs[s], winlen)
                    elif cmbx_chunk.current() == 2:
                        [IRend, env] = metodo_propio(signalsdb[s][c][b], fs[s], winlen)
                    # Smooth method
                    if cmbx_envelope.current() == 0:
                        h2 = 10 ** (signalsdb[s][c][b][0:IRend] / 10)
                        env = 10 * np.log10(np.flip(np.cumsum(np.flip(h2))) / np.sum(h2))
                    elif cmbx_envelope.current() == 1:
                        if IRend + 1 + winlen < len(signalsdb[s][c][b]):
                            env = ndimage.median_filter(signalsdb[s][c][b][0:IRend + 1 + winlen], winlen)
                        else:
                            env = ndimage.median_filter(signalsdb[s][c][b], winlen)
                        env = env[int(winlen / 2):int(len(env) - winlen / 2)]
                    elif cmbx_envelope.current() == 2:
                        if cmbx_chunk.current() == 2:
                            env = env[0:IRend + 1]
                        else:
                            if IRend + 1 + winlen < len(signalsdb[s][c][b]):
                                env = mediamovil(signalsdb[s][c][b][0:IRend + winlen], winlen)
                            else:
                                env = mediamovil(signalsdb[s][c][b], winlen)
                    list.append(envelopes[s][c], env)
                    # Acoustic parameters
                    
                    #EDT
                    edtsamples = getsamplesbetween(envelopes[s][c][b], 0, -10)
                    results[0, b, s, c] = -60 / (slope(edtsamples) * fs[s])
                    #T20
                    t20samples = getsamplesbetween(envelopes[s][c][b], -5, -25)
                    results[1, b, s, c] = -60 / (slope(t20samples) * fs[s])
                    #T30
                    t30samples = getsamplesbetween(envelopes[s][c][b], -5, -35)
                    results[2, b, s, c] = -60 / (slope(t30samples) * fs[s])
                    #C50
                    c50i = np.int(0.05 * fs[s] + 1)
                    num = np.sum(signals[s][c][b][0:c50i] ** 2)
                    den = np.sum(signals[s][c][b][c50i:] ** 2)
                    results[3, b, s, c] = 10 * np.log10(num / den)
                    #C80
                    c80i = np.int(0.08 * fs[s] + 1)
                    num = np.sum(signals[s][c][b][0:c80i] ** 2)
                    den = np.sum(signals[s][c][b][c80i:] ** 2)
                    results[4, b, s, c] = 10 * np.log10(num / den)
    
                    # Transition time.
                    hil = np.abs(sc.hilbert(signals[s][c][b]))
                    before = hil[0:-2]
                    current = hil[1:-1]
                    after = hil[2:]
                    # Determinates the direct sound for window it.
                    minindex = np.where((before > current) * (after > current))[0]
                    minindex = minindex[minindex + 1 > np.argmax(hil)]
                    minval = hil[minindex + 1]
                    before = minval[0:-2]
                    current = minval[1:-1]
                    after = minval[2:]
                    minindex2 = np.where((before > current) * (after > current))[0][0]
                    index = minindex[minindex2 + 1] + 1
                    # Cummulative energy of the RIR outliers.
                    if cmbx_envelope.current() == 1:
                        cumenergy = np.cumsum( 10**( (signalsdb[s][c][b][index:len(envelopes[s][c][b])] + envelopes[s][c][b][index:])/10 ))
                    else:
                        medianmf = ndimage.median_filter(signalsdb[s][c][b][0:IRend + 1 + winlen], winlen)
                        medianmf = medianmf[int(winlen / 2):int(len(medianmf) - winlen / 2)]
                        cumenergy = np.cumsum( 10**( (signalsdb[s][c][b][index:len(medianmf)] + medianmf[index:])/10 ))
                    cumenergy = cumenergy / np.max(cumenergy)
                    Ttindex = np.where(cumenergy >= 0.99)[0][0] + index
                    results[5, b, s, c] = Ttindex / fs[s]
                    # EDTt
                    results[6, b, s, c] = -60 / (slope(envelopes[s][c][b][0:Ttindex + 1]) * fs[s])
                    # CTtt
                    num = np.sum(signals[s][c][b][0:Ttindex] ** 2)
                    den = np.sum(signals[s][c][b][Ttindex:] ** 2)
                    results[7, b, s, c] = 10 * np.log10(num / den)
                    # IACC and IACC early
                    if c == 1:
                        iacc = np.correlate(signals[s][0][b], signals[s][1][b]) / np.sqrt(
                                    np.sum(signals[s][0][b] ** 2) * np.sum(signals[s][1][b] ** 2))
                        results[8, b, s, 0] = np.max(np.abs(iacc))
                        results[8, b, s, 1] = results[8, b, s, 0]
    
                        iacc_early = np.correlate(signals[s][0][b][0:c80i], signals[s][1][b][0:c80i]) / np.sqrt(
                                    np.sum(signals[s][0][b][0:c80i] ** 2) * np.sum(signals[s][1][b][0:c80i] ** 2))
                        results[9, b, s, 0] = np.max(np.abs(iacc_early))
                        results[9, b, s, 1] = results[9, b, s, 0]
                        
                        Ttchannelaverage = int( fs[s]*(results[5, b, s, 1]+results[5, b, s, 0])/2)
                        iacct = np.correlate(signals[s][0][b][0:Ttchannelaverage], signals[s][1][b][0:Ttchannelaverage]) / np.sqrt(np.sum(signals[s][0][b][0:Ttchannelaverage] ** 2) * np.sum(signals[s][1][b][0:Ttchannelaverage] ** 2))
                        results[10, b, s, 0] = np.max(np.abs(iacct))
                        results[10, b, s, 1] = results[10, b, s, 0]
                        
                    # Charge the progress bar.
                    bar_step += 100 / total_steps
                    progress.set(bar_step)
                    APC.update()
        # Averages
        averages = np.zeros([len(rowsaverage), len(bands) + 1, len(parameters) - 3])
        averages[0, :, :] = np.transpose(np.sum(results[0:-3, :, :, :], axis=(2, 3)) / sum(channels))
        averages[1, :, :] = np.transpose(np.max(results[0:-3, :, :, :], axis=(2, 3)))
        minresul = results[:, :, :, 0] + 0
        minindex = (results[:, :, :, 0] > results[:, :, :, 1]) * (results[:, :, :, 1] != 0)
        minresul[minindex] = results[minindex, 1]
        averages[2, :, :] = np.transpose(np.min(minresul[0:-3, :], axis=2))
        averages[3, :, :] = np.transpose(np.std(results[0:-3, :, :, 0], axis=2))
        # Upload data to tables
        for i in range(len(rowsaverage) + 1):
            for j in range(len(bands_labels[-1]) + 1):
                if j <= len(bands):
                    if i == 0 and j != 0:
                        table1[i][j + 1]["text"] = strbands[j - 1]
                        table1[i][j + 1].grid()
                    elif i != 0:
                        table1[i][j + 1]["text"] = str(np.round(averages[i - 1, j, 0], 2))
                        table1[i][j + 1].grid()
                else:
                    table1[i][j + 1].grid_remove()
    
        for i in range(len(parameters) + 1):
            for j in range(len(bands_labels[-1]) + 1):
                if j <= len(bands):
                    if i == 0 and j != 0:
                        table2[i][j + 1]["text"] = strbands[j - 1]
                        table2[i][j + 1].grid()
                    elif i != 0:
                        table2[i][j + 1]["text"] = str(np.round(results[i - 1, j, 0, 0], 2))
                        table2[i][j + 1].grid()
                else:
                    table2[i][j + 1].grid_remove()
        #Update widgets and plots
        cmbx_bandplt["values"] = ("global",) + strbands
        names = ()
        for ir in range(0, len(lstbx_IRs.get(0, "end"))):
            names += tuple([os.path.basename(lstbx_IRs.get(ir))])
        cmbx_IRplt["values"] = names
        cmbx_bandplt.current(0)
        cmbx_IRplt.current(0)
        cmbx_channels["values"] = tuple(range(1, channels[0] + 1))
        cmbx_channels.current(0)
        table1[0][0].current(0)
        refresh_graphtable1("")
        refresh_graph2()
    except:
        messagebox.showerror("Error", message="Calculation fail. Please check data.")
    progress.set(0)
    btn_calculate.grid()
    pbar.grid_remove()
    changestate(1)
    return


def func_exptable():
    # Export the table to .Csv
    if len(cmbx_IRplt["values"]) != 0:
        filename = filedialog.asksaveasfilename(title="Save current table", filetypes=[("CSV File", "*.csv")],
                                                defaultextension=[("CSV File", "*.csv")])
        if not filename: return
        with open(filename, mode='w') as file:
            writer = csv.writer(file, delimiter=',')
            bandsrow = [''] + list(cmbx_bandplt["values"][:])
            writer.writerow(bandsrow)
            if btn_show["text"] == "Show average":
                res = results[:, :, cmbx_IRplt.current(), cmbx_channels.current()]
                for k in range(len(parameters)):
                    row = [parameters[k]] + list(res[k])
                    writer.writerow(row)
                row = ["RIR:"] + [cmbx_IRplt.get()] + ["Channel:"] + [cmbx_channels.get()]
                writer.writerow(row)
            else:
                ave = averages[:, :, table1[0][0].current()]
                for k in range(len(rowsaverage)):
                    row = [rowsaverage[k]] + list(ave[k])
                    writer.writerow(row)
                row = ["Parameter:"] + [table1[0][0].get()] + ["RIRs:"] + [cmbx_IRplt["values"]]
                writer.writerow(row)
    return


def changestate(onoff):
    #Function that disabled settings widgets to avoid using they during calculations
    if onoff == 1:
        btn_clearall["state"] = NORMAL
        btn_clearselected["state"] = NORMAL
        btn_load["state"] = NORMAL
        btn_sweep["state"] = NORMAL
        cmbx_bandfilter["state"] = "readonly"
        cmbx_chunk["state"] = "readonly"
        cmbx_envelope["state"] = "readonly"
        cmbx_maxband["state"] = "readonly"
        cmbx_minband["state"] = "readonly"
    else:
        btn_clearall["state"] = DISABLED
        btn_clearselected["state"] = DISABLED
        btn_load["state"] = DISABLED
        btn_sweep["state"] = DISABLED
        cmbx_bandfilter["state"] = DISABLED
        cmbx_chunk["state"] = DISABLED
        cmbx_envelope["state"] = DISABLED
        cmbx_maxband["state"] = DISABLED
        cmbx_minband["state"] = DISABLED
    return


def slope(data): 
    # Given a set of samples as input, this function calculate the linear regression slope.
    # Unit: data units per sample
    n = len(data) / 1
    x = np.arange(0, n) / 1
    m = (n * np.sum(x * data) - np.sum(x) * np.sum(data)) / (n * np.sum(x ** 2) - np.sum(x) ** 2)
    return m


def getsamplesbetween(data, Ls, Li):
    # This function returns the data samples <Ls and >Li, with Ls and Li the levels relative to the maximum
    supe = np.where(data >= np.max(data) + Ls)[0][-1]
    infe = np.where(data < np.max(data) + Li)[0]
    infe = infe[infe > supe][0]
    samples = data[supe + 1:infe]
    return samples


def mediamovil(data, k):
    # returns Media Moving Filter of data, with window lenght k. 
    # output lenght results: input lenght - window lenght +1.
    w = np.ones([k]) / k
    mmf = sc.fftconvolve(w,data)
    mmf = mmf[k-1:len(mmf) - (k - 1)]
    
    return mmf


def metodo_propio(data, fs, k):
    # Roge method for truncation time of IR.
    mmf = mediamovil(data, k) #media moving filter
    mmf = np.concatenate([mmf, np.ones([fs]) * mmf[-1]]) # append 1s with end sample value, to improve detection
    indexmax = np.argmax(mmf)
    levelmax = np.max(mmf)
    M = (mmf[-1] - levelmax) / (len(mmf) - indexmax) # slope of interpolation line between max and end sample
    B = (levelmax - M * indexmax) #intercept of interpolation line
    cut = np.argmax(M * range(len(mmf))[indexmax:] + B - mmf[indexmax:]) + indexmax  # maximum distance between mmf and interpolation line
    cut = np.int(cut)
    mmf = np.concatenate([mmf[0:cut], np.ones([fs]) * mmf[cut]]) # repeat the process to improve detection 
    M = (mmf[-1] - levelmax) / (len(mmf) - indexmax)
    B = (levelmax - M * indexmax)
    cut = np.argmax(M * range(len(mmf))[indexmax:] + B - mmf[indexmax:]) + indexmax
    cut = np.int(cut)
    return cut, mmf


def lundeby(val2, fs):
    t1 = np.arange(len(val2))/fs
    n = int(np.floor(len(val2)/fs/0.01))   # number of samples on 1 interval.
    v = int(np.floor(len(val2)/n))         # number of intervals.
    val_fil = mediamovil(val2,v)
    t_fil = np.arange(len(val_fil))/fs
    noise = val2[int(0.9*len(val2)):]
    noise = 10*np.log10(np.mean(10**(noise/10)))       # Find the noise in dB from the last 10% of signal

    # Linear regression from max dB to noise level (plus 5-10 dB)
    init = val_fil[np.abs(val_fil - max(val_fil)).argmin()]
    init_sample = np.where(val_fil == init)[0][0]
    end = val_fil[np.where(val_fil < (noise+5))[0][0]]
    end_sample = np.where(val_fil < (noise+5))[0][0]
    slope = (val_fil[end_sample] - val_fil[init_sample]) / (t_fil[end_sample] - t_fil[init_sample])
    intercept = end - slope*t_fil[end_sample]

    # find the preliminary crosspoint.
    xpoint = (noise-intercept)/slope
    tpoint = int(np.abs(t1 - (xpoint)).argmin())
    # new interval time (between 3-10 (p=6) intervals per 10 dB decay)
    p = 6
    max_ite = 5           # Number of max iterations (5 must be enough)
    ite = 1
    error = 1
    # Iterative process.
    while (error > 0.0001) and (ite <= max_ite):
        delta = abs(10/slope)                     # time to get -10 dB decay
        int_time = delta/p                        # interval time
        n = int(np.floor(delta*fs/p))             # number of samples in one interval
        v = int(np.floor(len(t1)/n))
        if v < 2:                               # Case when v is less than 2
            v = 2
        val_fil = mediamovil(val_fil, v)
        t_fil = np.linspace(0, len(val2) / fs, num=len(val_fil))
        delta_t = t_fil[np.abs(t_fil - xpoint - delta/2).argmin()]                    # Where occurs crosspoint + 5 dB.
        noise = val_fil[np.where(t_fil == delta_t)[0][0]::]
        if len(noise) < round(0.1*len(t1)):
            noise = val_fil[len(val_fil) - int(0.1*len(val_fil))::]
        noise = 10*np.log10(np.mean(10**(noise/10)))  # Find the noise in dB from the last 10% of signal

        # Linear regression from 0 dB to noise level (plus 5-10 dB)
        init = val_fil[np.abs(val_fil - max(val_fil)).argmin()]
        init_sample = np.where(val_fil == init)[0][0]
        end = val_fil[np.where(val_fil < (noise+5))[0][0]]
        end_sample = np.where(val_fil < (noise+5))[0][0]
        slope = (val_fil[end_sample] - val_fil[init_sample]) / (t_fil[end_sample] - t_fil[init_sample])
        intercept = end - slope * t_fil[end_sample]

        # find the new crosspoint and error.
        error = abs(xpoint-(noise-intercept)/slope)/xpoint
        xpoint = (noise - intercept) / slope
        ite = ite + 1
    xpointt = int(np.abs(t1 - (xpoint)).argmin())
    mmf = np.concatenate([val_fil[0:xpointt], np.ones([fs]) * val_fil[xpointt]])
    return xpointt, mmf


def pepino(val,fs,k):
    # Creation of initial grids
    t1 = np.linspace(0, len(val) / fs, num=len(val))
    ite = 0
    max_ite = 12
    N = 20
    T = max(t1)
    x3s=0
    R_obs=0
    R = T/N
    error = 1
    # Iterative process
    while (error > 0.0001) and (ite <= max_ite):
        x3 = np.arange(x3s-2*R_obs,x3s+2*R_obs,R)
        if ite==0:
            x3 = np.arange(0,max(t1),R)
        x2 = max(val)
        y_k = np.zeros((N,len(t1)))
        e = np.zeros(N)
        x3s = 0
        for i in range(N):
            if i==0:
                y_k[i][:]=x2
            else:
                init = val[np.abs(val - x2).argmin()]
                end = val[np.abs(t1 - x3[i]).argmin()]
                init_sample = np.where(val == init)[0][0]
                end_sample = np.abs(t1 - x3[i]).argmin()
                x = np.arange(init_sample, end_sample + 1) / fs
                y = val[init_sample:end_sample + 1]
                slope, intercept = stats.linregress(x, y)[0:2]
                x1 = slope
                y_k[i][:] = x1 * t1 + x2
                y_k[i][end_sample::] = x1 * x3[i] + x2
            e[i] = np.square(np.subtract(y_k[i], val)).mean()
        error = abs(x3s - (x3[np.where(e == min(e))[0][0]])) / x3s
        x3s = x3[np.where(e == min(e))[0][0]]
        ite = ite + 1
        R_obs = R
        R = (4*R)/N
        y_f = y_k[np.where(e == min(e))[0][0]]
        cut_db = y_f[np.abs(t1 - x3s).argmin()] - 5
        cut = np.abs(y_f - cut_db).argmin()
    mmf = np.concatenate([val[0:cut], np.ones([fs]) * val[cut]])
    return cut, mmf


def func_load():
    # load multiple RIRs, storage sample rate, number of channels, and audio data, and upload RIRs filename
    name = filedialog.askopenfilenames(title="Load RIRs files", filetypes=[("WAV Audio", "*.wav")])
    try:
        for k in range(len(name)):
            [sgn, frec] = sf.read(name[k])
            if len(np.shape(sgn)) == 2:
                list.append(RIRchannels, np.shape(sgn)[1])
            elif len(np.shape(sgn)) == 1:
                list.append(RIRchannels, 1)
                sgn.resize([len(sgn), 1])
            else:
                messagebox.showerror("Error", message="Data loaded has to many channels. Maximum channels is 2.")
                return
            lstbx_IRs.insert("end", name[k])
            list.append(RIRs, sgn)
            list.append(samplerate, frec)
    except:
        messagebox.showerror("Error", message="Load fail.")
    return


def func_sweep():
    # load sine sweep measurment of RIRs and deconvolve using the inverse filter through FFT method.
    # storage sample rate, number of channels, and upload filename
    SSname = filedialog.askopenfilenames(title="Load Sweep files", filetypes=[("WAV Audio", "*.wav")])
    IFname = filedialog.askopenfilename(title="Load Inverse filter file", filetypes=[("WAV Audio", "*.wav")])
    try:
        if len(IFname) != 0 and len(SSname) != 0:
            [ifilt, fs_ifilt] = sf.read(IFname)
            if len(np.shape(ifilt)) == 2:
                messagebox.showerror("Error", message="Inverse filter must have 1 channel.")
                return
            for k in range(len(SSname)):
                [ss, fs_ss] = sf.read(SSname[k])
                if fs_ss != fs_ifilt:
                    messagebox.showerror("Error", message="Sample rate missmatch")
                    return
                if len(np.shape(ss)) == 2:
                    sgn = np.zeros([np.abs(np.shape(ss)[0] - len(ifilt)) + 1, np.shape(ss)[1]])
                    for c in range(np.shape(ss)[1]):
                        sgn[:, c] = sc.fftconvolve(ifilt, ss[:, c], mode="valid")
                        sgn[:, c] = sgn[:, c] / np.max(np.abs(sgn[:, c]))
                    list.append(RIRchannels, np.shape(sgn)[1])
                elif len(np.shape(ss)) == 1:
                    sgn = sc.fftconvolve(ifilt, ss, mode="valid")
                    sgn = sgn / np.max(np.abs(sgn))
                    list.append(RIRchannels, 1)
                    sgn.resize([len(sgn), 1])
                else:
                    messagebox.showerror("Error", message="Data loaded has to many channels. Maximum channels is 2.")
                    return
                lstbx_IRs.insert("end", SSname[k])
                list.append(RIRs, sgn)
                list.append(samplerate, fs_ss)
    except:
        messagebox.showerror("Error", message="Load fail.")
    return


def func_clearall():
    #clear all the RIRs loaded.
    list.clear(RIRs)
    list.clear(samplerate)
    list.clear(RIRchannels)
    lstbx_IRs.delete(0, "end")
    return


def func_clearselected():
    # clear the selected RIRs loaded
    for k in lstbx_IRs.curselection():
        list.pop(RIRs, lstbx_IRs.curselection()[0])
        list.pop(samplerate, lstbx_IRs.curselection()[0])
        list.pop(RIRchannels, lstbx_IRs.curselection()[0])
        lstbx_IRs.delete(lstbx_IRs.curselection()[0])
    return


def func_expplot():
    # Export the figure that is currently plotted.
    if len(cmbx_IRplt["values"]) != 0:
        filename = filedialog.asksaveasfile(title="Save current plot", filetypes=[("PNG Image", "*.png")],
                                            defaultextension=[("PNG Image", "*.png")])
        if not filename: return
        if btn_show["text"] == "Show average":
            figgraph2.savefig(filename.name, facecolor="w")
        else:
            figgraph1.savefig(filename.name, facecolor="w")
    return


def func_show():
    # switch the plot and table corresponding to RIR view or RIRs average.
    if btn_show["text"] == "Show average":
        btn_show["text"] = "Show RIRs"
        frame_pltsettings.grid_remove()
        get_widz2.grid_remove()
        frame_table2.grid_remove()
        get_widz1.grid()
        frame_table1.grid()
    else:
        btn_show["text"] = "Show average"
        frame_pltsettings.grid()
        get_widz2.grid()
        frame_table2.grid()
        get_widz1.grid_remove()
        frame_table1.grid_remove()
    return


def func_bandfilter(event):
    # Update widgets when user change the band filter
    if cmbx_bandfilter.current() == 0:
        cmbx_minband["values"] = bands_labels[0]
        cmbx_maxband["values"] = bands_labels[0]
        cmbx_minband.current(0)
        cmbx_maxband.current(len(bands_labels[0]) - 1)
    elif cmbx_bandfilter.current() == 1:
        cmbx_minband["values"] = bands_labels[1]
        cmbx_maxband["values"] = bands_labels[1]
        cmbx_minband.current(0)
        cmbx_maxband.current(len(bands_labels[1]) - 1)
    return


def func_limitband(event):
    # Update the posible lower and higher bands, when user change they.
    if cmbx_bandfilter.current() == 0:
        cmbx_maxband["values"] = bands_labels[0][cmbx_minband.current():]
        cmbx_minband["values"] = bands_labels[0][0:cmbx_maxband.current() + cmbx_minband.current() + 1]
    elif cmbx_bandfilter.current() == 1:
        cmbx_maxband["values"] = bands_labels[1][cmbx_minband.current():]
        cmbx_minband["values"] = bands_labels[1][0:cmbx_maxband.current() + cmbx_minband.current() + 1]
    return


def refresh_graphtable1(event):
    #Update plot and table averages, when user change current parameter.
    if len(cmbx_IRplt["values"]) != 0:
        aver = averages[0, :, table1[0][0].current()]
        maxi = averages[1, :, table1[0][0].current()]
        mini = averages[2, :, table1[0][0].current()]
        desvi = averages[3, :, table1[0][0].current()]
        f = range(len(cmbx_bandplt["values"]))
        graph1.cla()
        graph1.plot(aver[0], "ok")
        #graph1.errorbar(0, aver[0], yerr=np.abs(np.array([[mini[0]], [maxi[0]]]) - aver[0]), fmt='--o')
        graph1.errorbar(0, aver[0], yerr=desvi[0]*2, fmt='--o')
        #graph1.fill_between(f[1:], mini[1:], maxi[1:], facecolor='b', alpha=0.3)
        graph1.fill_between(f[1:],aver[1:]-desvi[1:]*2,aver[1:]+desvi[1:]*2, facecolor='b', alpha=0.3)
        graph1.plot(f[1:], aver[1:], "r")
        graph1.set_xlabel("Frequency [Hz]")
        graph1.set_xticks(f)
        graph1.set_xticklabels(cmbx_bandplt["values"])
        graph1.grid()
        if table1[0][0].current() == 0:
            graph1.set_ylabel("EDT [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current() == 1:
            graph1.set_ylabel("T20 [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current() == 2:
            graph1.set_ylabel("T30 [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current() == 3:
            graph1.set_ylabel("C50 [dB]")
        elif table1[0][0].current() == 4:
            graph1.set_ylabel("C80 [dB]")
        elif table1[0][0].current() == 5:
            graph1.set_ylabel("Tt [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current() == 6:
            graph1.set_ylabel("EDTt [s]")
            graph1.set_ylim(ymin=0)
        canv1.draw()
        for i in range(len(rowsaverage)):
            for j in range(len(cmbx_bandplt["values"])):
                table1[i + 1][j + 1]["text"] = str(np.round(averages[i, j, table1[0][0].current()], 2))

    return


def func_channels(event):
    # Callback function for update RIR plot and table, when user change current channel
    if len(cmbx_IRplt["values"]) != 0:
        refresh_table2()
        refresh_graph2()
    return


def func_bandplot(event):
    # Callback function for update plot when user change band
    if len(cmbx_IRplt["values"]) != 0:
        refresh_graph2()
    return


def func_IRplot(event):
    # Callback function for update plot and table when user change RIR
    if len(cmbx_IRplt["values"]) != 0:
        cmbx_channels["values"] = tuple(range(1, channels[cmbx_IRplt.current()] + 1))
        cmbx_channels.current(0)
        refresh_table2()
        refresh_graph2()
    return


def refresh_table2():
    # Update RIR table when is required for other functions
    for i in range(len(parameters)):
        for j in range(len(cmbx_bandplt["values"])):
            table2[i + 1][j + 1]["text"] = str(
                np.round(results[i, j, cmbx_IRplt.current(), cmbx_channels.current()], 2))
    return


def refresh_graph2():
    # Update RIR plot when is required for other functions
    sgn = signalsdb[cmbx_IRplt.current()][cmbx_channels.current()][cmbx_bandplt.current()]
    env = envelopes[cmbx_IRplt.current()][cmbx_channels.current()][cmbx_bandplt.current()]
    frecs = fs[cmbx_IRplt.current()]
    graph2.cla()
    if len(sgn) < 2 * len(env):
        end = len(sgn)
    else:
        end = len(env) * 2
    graph2.plot(np.arange(0, end) / frecs, sgn[0:end])
    graph2.plot(np.arange(0, len(env)) / frecs, env)
    graph2.set_xlim(xmax=end / frecs)
    graph2.set_xlabel("Time [s]")
    graph2.set_ylabel("Level [dBFS]")
    graph2.grid()
    canv2.draw()
    return

# Allocation for RIRs data.
RIRs = []
samplerate = []
RIRchannels = []
# Center frequencies of band filters and they corresponding label .
bands_centers = [(31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000), (
                25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
                2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000)]
bands_labels = [("31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"), (
                "25", "31.5", "40", "50", "63", "80", "100", "125", "160", "200", "250","315", "400", "500", "630",
                "800", "1k","1.25k", "1.6k", "2k", "2.5k", "3.15k", "4k", "5k", "6.3k", "8k", "10k", "12.5k", "16k",
                "20k")]
#Parameters  
parameters = ("EDT", "T20", "T30", "C50", "C80", "Tt", "EDTt","CTt", "IACC", "IACC early", "IACCt")
rowsaverage = ("Average", "Max", "Min", "Sigma")

# Main frame of GUI
APC = Tk()
APC.title("Acoustic Parameters Calculator")
APC.rowconfigure(0, weight=1)
APC.columnconfigure(1, weight=1)
APC["padx"] = "2"
APC["pady"] = "2"
APC.state("zoomed")

# Figure of average plot
figgraph1 = Figure(constrained_layout=True)
graph1 = figgraph1.add_subplot(111)
graph1.set_xlabel("Frecuency [Hz]")
graph1.set_ylabel("")
r, g, b = APC.winfo_rgb(APC["bg"])
figgraph1.set_facecolor([r / 65536, g / 65536, b / 65536])
# place in GUI
canv1 = FigureCanvasTkAgg(figgraph1, master=APC)
canv1.draw()
get_widz1 = canv1.get_tk_widget()
get_widz1.grid(row=0, column=1, padx=4, pady=4, sticky="ewns")
get_widz1.grid_remove()

# Figure of RIR plot
figgraph2 = Figure(constrained_layout=True)
graph2 = figgraph2.add_subplot(111)
graph2.set_xlabel("Time [s]")
graph2.set_ylabel("Level [dBFS]")
graph2.set_ylim([-100, 0])
r, g, b = APC.winfo_rgb(APC["bg"])
figgraph2.set_facecolor([r / 65536, g / 65536, b / 65536])
# place in GUI
canv2 = FigureCanvasTkAgg(figgraph2, master=APC)
canv2.draw()
get_widz2 = canv2.get_tk_widget()
get_widz2.grid(row=0, column=1, padx=4, pady=4, sticky="ewns")


##Frame in which menu is located
frame_menu = Frame(APC)
frame_menu.grid(row=0, column=0, sticky="ewns", padx=2, pady=2)
frame_menu.rowconfigure(3, weight=1)

#Progress bar to show while calculating
progress = DoubleVar()
pbar = Progressbar(frame_menu, variable=progress)
pbar.grid(row=4, column=0, columnspan=2, sticky="ewns", padx=2, pady=2)

#Buttons:

#Load RIRs button
btn_load = Button(frame_menu, text="Load RIRs", command=func_load)
btn_load.grid(row=0, column=0, padx=2, pady=2, sticky="ewns")

# Load sinesweep button
btn_sweep = Button(frame_menu, text="Load Sweep", command=func_sweep)
btn_sweep.grid(row=0, column=1, padx=2, pady=2, sticky="ewns")

#Clear all RIRs button
btn_clearall = Button(frame_menu, text="Clear all", command=func_clearall)
btn_clearall.grid(row=1, column=0, padx=2, pady=2, sticky="ewns")

#Clear selected RIRs button
btn_clearselected = Button(frame_menu, text="Clear selected", command=func_clearselected)
btn_clearselected.grid(row=1, column=1, padx=2, pady=2, sticky="ewns")

#Calculate button
btn_calculate = Button(frame_menu, text="Calculate", font=("Arial", 12, "bold"), command=func_calculate)
btn_calculate.grid(row=4, column=0, columnspan=2, sticky="ewns", padx=2, pady=2)

# Show average/Show RIRs button
btn_show = Button(frame_menu, text="Show average", command=func_show)
btn_show.grid(row=5, column=0, columnspan=2, sticky="ewns", padx=2, pady=2)

# Export current plot button
btn_expplot = Button(frame_menu, text="Export plot", command=func_expplot)
btn_expplot.grid(row=6, column=0, sticky="ewns", padx=2, pady=2)

# Export current table button
btn_exptable = Button(frame_menu, text="Export table", command=func_exptable)
btn_exptable.grid(row=6, column=1, sticky="ewns", padx=2, pady=2)


##Frame inside menu, in which the RIRs list and the scrollbars are located 
frame_IRs = Frame(frame_menu)
frame_IRs.grid(row=3, column=0, columnspan=2, sticky="ewns", padx=2, pady=2)
frame_IRs.columnconfigure(0, weight=1)
frame_IRs.rowconfigure(0, weight=1)
#Vertical scrollbar
scrollbarV = Scrollbar(frame_IRs, orient="vertical")
scrollbarV.grid(row=0, column=1, sticky="ewns")
#Horizontal scrollbar
scrollbarH = Scrollbar(frame_IRs, orient="horizontal")
scrollbarH.grid(row=1, column=0, sticky="ewns")
#Listbox of RIRs
lstbx_IRs = Listbox(frame_IRs, selectmode="extended")
lstbx_IRs.grid(row=0, column=0, sticky="ewns")
lstbx_IRs.config(yscrollcommand=scrollbarV.set, xscrollcommand=scrollbarH.set)
scrollbarV.config(command=lstbx_IRs.yview)
scrollbarH.config(command=lstbx_IRs.xview)


##Frame inside menu, in which the settings are located
frame_settings = Frame(frame_menu)
frame_settings.grid(row=2, column=0, columnspan=2, sticky="ewns", padx=2, pady=2)

# Settings main title
lbl_settings = Label(frame_settings, text="Settings", anchor="center")
lbl_settings.grid(row=0, column=0, columnspan=4, sticky="ewns")

# Band filter selection
lbl_bandfilter = Label(frame_settings, text="Band filter:", anchor="e")
lbl_bandfilter.grid(row=1, column=0, sticky="ewns")
cmbx_bandfilter = Combobox(frame_settings, values=("Octave", "Third"), state="readonly", width=7)
cmbx_bandfilter.grid(row=1, column=1, sticky="ewns", pady=1)
cmbx_bandfilter.bind("<<ComboboxSelected>>", func_bandfilter)

# Minimum band selection
lbl_minband = Label(frame_settings, text="Min band:", anchor="e")
lbl_minband.grid(row=2, column=0, sticky="ewns")
cmbx_minband = Combobox(frame_settings, state="readonly", width=7)
cmbx_minband.grid(row=2, column=1, sticky="ewns", pady=1)
cmbx_minband.bind("<<ComboboxSelected>>", func_limitband)

# Maximum band selection
lbl_maxband = Label(frame_settings, text="Max band:", anchor="e")
lbl_maxband.grid(row=3, column=0, sticky="ewns")
cmbx_maxband = Combobox(frame_settings, state="readonly", width=7)
cmbx_maxband.grid(row=3, column=1, sticky="ewns", pady=1)
cmbx_maxband.bind("<<ComboboxSelected>>", func_limitband)

# Smooth method selection
lbl_envelope = Label(frame_settings, text="Smooth:", anchor="e")
lbl_envelope.grid(row=1, column=2, sticky="ewns")
cmbx_envelope = Combobox(frame_settings, values=("Schoreder", "Median MF", "Mean MF"), state="readonly", width=10)
cmbx_envelope.current(0)
cmbx_envelope.grid(row=1, column=3, sticky="ewns", pady=1)

#Truncation method selection
lbl_chunk = Label(frame_settings, text="IR chunk:", anchor="e")
lbl_chunk.grid(row=2, column=2, sticky="ewns")
cmbx_chunk = Combobox(frame_settings, values=("Lundeby", "Pepino", "Roge"), state="readonly", width=10)
cmbx_chunk.current(0)
cmbx_chunk.grid(row=2, column=3, sticky="ewns", pady=1)


##Frame inside menu, in which the plot settings of calculated parameters are located
frame_pltsettings = Frame(frame_menu)
frame_pltsettings.grid(row=7, column=0, columnspan=2)

# Measurement selection to visualize the calculated parameters and the plot
lbl_IRplt = Label(frame_pltsettings, text="RIR:")
lbl_IRplt.grid(row=0, column=0, sticky="ens", pady=2)
cmbx_IRplt = Combobox(frame_pltsettings, state="readonly")
cmbx_IRplt.grid(row=0, column=1, columnspan=3, sticky="wns", pady=2, padx=2)
cmbx_IRplt.bind("<<ComboboxSelected>>", func_IRplot)

# Band selection to visualize the plot
lbl_bandplt = Label(frame_pltsettings, text="Band:")
lbl_bandplt.grid(row=1, column=0, sticky="ens", pady=2)
cmbx_bandplt = Combobox(frame_pltsettings, state="readonly", width=6)
cmbx_bandplt.grid(row=1, column=1, sticky="wns", pady=2, padx=2)
cmbx_bandplt.bind("<<ComboboxSelected>>", func_bandplot)

# Channel selection to visualize the calculated parameters and the plot
lbl_channels = Label(frame_pltsettings, text="Channel:")
lbl_channels.grid(row=1, column=2, sticky="ens", pady=2)
cmbx_channels = Combobox(frame_pltsettings, state="readonly", width=6)
cmbx_channels.grid(row=1, column=3, sticky="wns", pady=2, padx=2)
cmbx_channels.bind("<<ComboboxSelected>>", func_channels)


#Frame in which the table of averages is located
frame_table1 = Frame(APC)
frame_table1.grid(row=1, column=0, columnspan=2, padx=4, pady=4, sticky="ewns")
frame_table1.grid_remove()

#Table of averages
for i in range(len(rowsaverage) + 1):
    for j in range(len(bands_labels[-1]) + 2):
        if i == 0 and j > 1: # First row. Bands labels
            table1[i] = table1[i] + [
                Label(frame_table1, text=bands_labels[1][j - 2], bg="white", relief="solid", borderwidth=1)]
            table1[i][j].grid(row=i, column=j, sticky="ewns")
            frame_table1.columnconfigure(j, weight=1)
        elif i == 0 and j == 1: # "global" label cell
            table1[i] = table1[i] + [Label(frame_table1, text="global", bg="white", relief="solid", borderwidth=1)]
            table1[i][j].grid(row=i, column=j, sticky="ewns")
            frame_table1.columnconfigure(j, weight=1)
        elif i != 0 and j == 0: # First column. Rows headers with averages labels
            table1 = table1 + [
                [Label(frame_table1, text=rowsaverage[i - 1], bg="white", relief="solid", borderwidth=1)]]
            table1[i][j].grid(row=i, column=j, sticky="ewns")
        elif i != 0 and j != 0: # Empty cells in which calculated averages data will be placed
            table1[i] = table1[i] + [Label(frame_table1, text="", bg="white", relief="solid", borderwidth=1)]
            table1[i][j].grid(row=i, column=j, sticky="ewns")
        elif i == 0 and j == 0: # Upper left cell. Contains parameter selection
            table1 = [[Combobox(frame_table1, values=parameters[0:-3], width=5, state="readonly")]]
            table1[i][j].current(0)
            table1[i][j].bind("<<ComboboxSelected>>", refresh_graphtable1)
            table1[i][j].grid(row=i, column=j, sticky="ewns")


# Frame in which the table of calculated parameters of each measurement is located
frame_table2 = Frame(APC)
frame_table2.grid(row=1, column=0, columnspan=2, padx=4, pady=4, sticky="ewns")

# Table of calculated parameters
for i in range(len(parameters) + 1):
    for j in range(len(bands_labels[-1]) + 2):
        if i == 0 and j > 1: #First row. Bands labels
            table2[i] = table2[i] + [
                Label(frame_table2, text=bands_labels[1][j - 2], bg="white", relief="solid", borderwidth=1)]
            table2[i][j].grid(row=i, column=j, sticky="ewns")
            frame_table2.columnconfigure(j, weight=1)
        elif i == 0 and j == 1: #"global" label cell
            table2[i] = table2[i] + [Label(frame_table2, text="global", bg="white", relief="solid", borderwidth=1)]
            table2[i][j].grid(row=i, column=j, sticky="ewns")
            frame_table2.columnconfigure(j, weight=1)
        elif i != 0 and j == 0: #First column. Rows headers with parameters labels
            table2 = table2 + [[Label(frame_table2, text=parameters[i - 1], bg="white", relief="solid", borderwidth=1)]]
            table2[i][j].grid(row=i, column=j, sticky="ewns")
        elif i != 0 and j != 0: # Empty cells in which calculated parameters data will be placed
            table2[i] = table2[i] + [Label(frame_table2, text="", bg="white", relief="solid", borderwidth=1)]
            table2[i][j].grid(row=i, column=j, sticky="ewns")
        elif i == 0 and j == 0: #Upper left cell. Not used 
            table2 = [[""]]

APC.mainloop()
