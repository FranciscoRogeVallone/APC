from tkinter import * 
from tkinter import filedialog
from tkinter.ttk import Progressbar , Combobox , Treeview
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import soundfile as sf
from scipy import ndimage
import scipy.signal as sc
import numpy as np
import matplotlib.pyplot as plt


def func_calculate():
    if len(RIRs)==0 or cmbx_bandfilter.current()==-1 :
        return
    
    global results
    global signals
    global envelopes
    global averages
    
    strbands = bands_labels[cmbx_bandfilter.current()][cmbx_minband.current():cmbx_minband.current()+cmbx_maxband.current()+1]
    bands = bands_centers[cmbx_bandfilter.current()][cmbx_minband.current():cmbx_minband.current()+cmbx_maxband.current()+1]
    results = np.zeros([len(parameters),len(bands)+1,len(RIRs)])
    signals = []
    envelopes = []
    for s in range(len(RIRs)):
        
        if cmbx_bandfilter.current()==0:
            freclim = (2**(1/2))
        elif cmbx_bandfilter.current()==1:
            freclim = (2**(1/6))
        
        for b in range(len(bands)+1):
            
            if b==0:
                finf = (2/fs[s]) * bands[0]/freclim
                fsup = (2/fs[s]) * bands[-1]*freclim
                winlen = int(1000*np.log10(fs[s]/ (bands[0]/freclim)))
                if fsup>1:
                    fsup = (fs[s]-1)/fs[s]
                sos = sc.butter(8,[finf,fsup],btype="band",output="sos")
                if reversedcheck.get()==1:
                    sgn = sc.sosfiltfilt(sos,np.flip(RIRs[s]))
                    list.append(signals,[np.flip(sgn)])
                else:
                    sgn = sc.sosfiltfilt(sos,RIRs[s])
                    list.append(signals,[sgn])
            else:
                finf = (2/fs[s]) * bands[b-1]/freclim
                fsup = (2/fs[s]) * bands[b-1]*freclim
                winlen = int(1000*np.log10(fs[s]/ (bands[b-1]/freclim)))                
                if fsup>1:
                    fsup = (fs[s]-1)/fs[s]
                sos = sc.butter(8,[finf,fsup],btype="band",output="sos")
                if reversedcheck.get()==1:
                    sgn = sc.sosfiltfilt(sos,np.flip(RIRs[s]))
                    list.append(signals[s],np.flip(sgn))
                else:
                    sgn = sc.sosfiltfilt(sos,RIRs[s])
                    list.append(signals[s],sgn)
            
            signals[s][b] = np.abs(signals[s][b])
            IRstart = np.where(signals[s][b][0:np.argmax(signals[s][b])]<=np.max(signals[s][b])/10)[0]
            if len(IRstart)==0:
                IRstart=0
            else:
                IRstart=IRstart[-1]
            signals[s][b] = signals[s][b][IRstart:]
            signals[s][b][signals[s][b]==0]=0.00000001
            signals[s][b] = 20*np.log10(signals[s][b])
            
            if cmbx_chunk.current()==0:
                print(lundeby)
            elif cmbx_chunk.current()==1:
                print(pepino)
            elif cmbx_chunk.current()==2:
                [IRend,env] = metodo_propio(signals[s][b],fs[s],winlen)
            
            if cmbx_envelope.current()==0:
                h2 = 10**(signals[s][b][0:IRend]/10)
                env =10*np.log10( np.flip(np.cumsum(np.flip(h2)))/np.sum(h2) )
            elif cmbx_envelope.current()==1:
                if IRend+1+winlen<len(signals[s][b]):
                    env = ndimage.median_filter(signals[s][b][0:IRend+1+winlen],winlen)
                else:
                    env = ndimage.median_filter(signals[s][b],winlen)
                env = env[int(winlen/2):int(len(env)-winlen/2)]
            elif cmbx_envelope.current()==2:
                if cmbx_chunk.current()==2:
                    env = env[0:IRend+1]
                else:
                    if IRend+1+winlen<len(signals[s][b]):
                        env = mediamovil(signals[s][b][0:IRend+winlen],winlen)
                    else:
                        env = mediamovil(signals[s][b],winlen)          
            if b==0:
                list.append(envelopes,[env])
            else:
                list.append(envelopes[s],env)
            
            
            edtsamples = getsamplesbetween(envelopes[s][b],0,-10)
            results[0,b,s] = -60/(slope(edtsamples)*fs[s])
            
            t20samples = getsamplesbetween(envelopes[s][b],-5,-25)
            results[1,b,s] = -60/(slope(t20samples)*fs[s])
            
            t30samples = getsamplesbetween(envelopes[s][b],-5,-35)
            results[2,b,s] = -60/(slope(t30samples)*fs[s])
            
            c50i = np.int(0.05*fs[s]+1)
            num = np.sum(10**(signals[s][b][0:c50i]/10))
            den = np.sum(10**(signals[s][b][c50i:]/10))        
            results[3,b,s] = 10*np.log10(num/den)
             
            c80i = np.int(0.08*fs[s]+1)
            num = np.sum(10**(signals[s][b][0:c80i]/10))
            den = np.sum(10**(signals[s][b][c80i:]/10))        
            results[4,b,s] = 10*np.log10(num/den)

                
    averages = np.zeros([len(rowsaverage),len(bands)+1,len(parameters)-1])
    averages[0,:,:] = np.transpose(np.mean(results[0:-1,:,:],axis=2))
    averages[1,:,:] = np.transpose(np.max(results[0:-1,:,:],axis=2))
    averages[2,:,:] = np.transpose(np.min(results[0:-1,:,:],axis=2))
    averages[3,:,:] = np.transpose(np.var(results[0:-1,:,:],axis=2))


    for i in range(len(rowsaverage)+1):
        for j in range(len(bands_labels[-1])+1):
            if j<=len(bands):
                if i==0 and j!=0:
                    table1[i][j+1]["text"] = strbands[j-1]
                    table1[i][j+1].grid()
                elif i!=0:
                    table1[i][j+1]["text"] = str(np.round(averages[i-1,j,0],2)) 
                    table1[i][j+1].grid()
            else:
                table1[i][j+1].grid_remove()
    
    for i in range(len(parameters)+1):
        for j in range(len(bands_labels[-1])+1):
            if j<=len(bands):
                if i==0 and j!=0:
                    table2[i][j+1]["text"] = strbands[j-1]
                    table2[i][j+1].grid()
                elif i!=0:
                    table2[i][j+1]["text"] = str(np.round(results[i-1,j,0],2)) 
                    table2[i][j+1].grid()
            else:
                table2[i][j+1].grid_remove()
    
    
    cmbx_bandplt["values"] = ("global",)+strbands
    cmbx_IRplt["values"] = lstbx_IRs.get(0,"end")
    cmbx_bandplt.current(0)
    cmbx_IRplt.current(0)
    table1[0][0].current(0)
    refresh_graphtable1("")
    refresh_graph2("")
    return

def slope(data):
    n = len(data)/1
    x = np.arange(0,n)/1
    m = ( n*np.sum(x*data)-np.sum(x)*np.sum(data) ) / ( n*np.sum(x**2)-np.sum(x)**2 )
    return m

def getsamplesbetween(data,Ls,Li):
    supe = np.where(data>=np.max(data)+Ls)[0][-1]
    infe = np.where(data<np.max(data)+Li)[0]
    infe = infe[infe>supe][0]
    samples = data[supe+1:infe]
    return samples

def mediamovil(data,k):
    w = np.concatenate([np.zeros([len(data)-1]),np.ones([k])/k])
    dataz = np.concatenate([data,np.zeros([k-1])])
    fftdataz = np.fft.rfft(dataz)
    fftw = np.fft.rfft(w)
    mmf = np.fft.irfft(fftdataz*fftw)
    mmf = mmf[0:len(data)-k+1]
    return mmf

def metodo_propio(data,fs,k):
    mmf = mediamovil(data, k)
    mmf = np.concatenate([mmf,np.ones([fs])*mmf[-1]])
    indexmax = np.argmax(mmf)
    levelmax = np.max(mmf)
    M = (mmf[-1]-levelmax)/(len(mmf)-indexmax)
    B = (levelmax-M*indexmax)
    cut = np.argmax(M*range(len(mmf))[indexmax:]+B-mmf[indexmax:])+indexmax
    cut = np.int(cut)
    print(np.max(mmf)-mmf[cut])
    return cut , mmf

def func_load():
    name = filedialog.askopenfilenames(title="Load RIRs files",filetypes=[("WAV Audio",".wav")])
    for k in range(len(name)):
        lstbx_IRs.insert("end",name[k])
        [sgn,frec]=sf.read(name[k])
        list.append(RIRs,sgn)
        list.append(fs,frec)
    return

def func_clearall():
    list.clear(RIRs)
    list.clear(fs)
    lstbx_IRs.delete(0,"end")
    return

def func_clearselected():
    for k in lstbx_IRs.curselection():
        list.pop(RIRs,lstbx_IRs.curselection()[0])
        list.pop(fs,lstbx_IRs.curselection()[0])
        lstbx_IRs.delete(lstbx_IRs.curselection()[0])
    return

def func_show():
    if btn_show["text"]=="Show average":
        btn_show["text"]="Show RIRs"
        lbl_bandplt.grid_remove()
        lbl_IRplt.grid_remove()
        cmbx_IRplt.grid_remove()
        cmbx_bandplt.grid_remove()
        get_widz2.grid_remove()
        frame_table2.grid_remove()
        get_widz1.grid()
        frame_table1.grid()
    else:
        btn_show["text"]="Show average"
        lbl_bandplt.grid()
        lbl_IRplt.grid()
        cmbx_IRplt.grid()
        cmbx_bandplt.grid()
        get_widz2.grid()
        frame_table2.grid()
        get_widz1.grid_remove()
        frame_table1.grid_remove()
    return

def func_bandfilter(event):
    if cmbx_bandfilter.current()==0:
        cmbx_minband["values"] = bands_labels[0]
        cmbx_maxband["values"] = bands_labels[0]
        cmbx_minband.current(0)
        cmbx_maxband.current(len(bands_labels[0])-1)
    elif cmbx_bandfilter.current()==1:
        cmbx_minband["values"] = bands_labels[1]
        cmbx_maxband["values"] = bands_labels[1]
        cmbx_minband.current(0)
        cmbx_maxband.current(len(bands_labels[1])-1)
    return

def func_limitband(event):
    if cmbx_bandfilter.current()==0:
        cmbx_maxband["values"] = bands_labels[0][cmbx_minband.current():]
        cmbx_minband["values"] = bands_labels[0][0:cmbx_maxband.current()+cmbx_minband.current()+1]
    elif cmbx_bandfilter.current()==1:
        cmbx_maxband["values"] = bands_labels[1][cmbx_minband.current():]
        cmbx_minband["values"] = bands_labels[1][0:cmbx_maxband.current()+cmbx_minband.current()+1]
    return

def refresh_graphtable1(event):
    if len(cmbx_IRplt["values"])!=0:
        aver =  averages[0,:,table1[0][0].current()]
        maxi = averages[1,:,table1[0][0].current()]
        mini = averages[2,:,table1[0][0].current()]
        f = range(len(cmbx_bandplt["values"]))
        graph1.cla()
        graph1.plot(aver[0],"ok")
        graph1.errorbar(0,aver[0],yerr=np.abs(np.array([[mini[0]],[maxi[0]]])-aver[0]),fmt='--o')
        graph1.fill_between(f[1:],mini[1:],maxi[1:],facecolor='b',alpha=0.3)
        graph1.plot(f[1:],aver[1:],"r")
        graph1.set_xlabel("Frecuency [Hz]")
        graph1.set_xticks(f)
        graph1.set_xticklabels(cmbx_bandplt["values"])
        graph1.grid()
        if table1[0][0].current()==0:
            graph1.set_ylabel("EDT [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current()==1:
            graph1.set_ylabel("T20 [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current()==2:
            graph1.set_ylabel("T30 [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current()==3:
            graph1.set_ylabel("C50 [dB]")
        elif table1[0][0].current()==4:
            graph1.set_ylabel("C80 [dB]")
        elif table1[0][0].current()==5:
            graph1.set_ylabel("Tt [s]")
            graph1.set_ylim(ymin=0)
        elif table1[0][0].current()==6:          
            graph1.set_ylabel("EDTt [s]")
            graph1.set_ylim(ymin=0)
        canv1.draw()
        for i in range(len(rowsaverage)):
            for j in range(len(cmbx_bandplt["values"])):
                table1[i+1][j+1]["text"] = str(np.round(averages[i,j,table1[0][0].current()],2)) 
    
    return

def refresh_table2(event):
    if len(cmbx_IRplt["values"])!=0:
        refresh_graph2("")
        for i in range(len(parameters)):
            for j in range(len(cmbx_bandplt["values"])):
                table2[i+1][j+1]["text"] = str(np.round(results[i,j,cmbx_IRplt.current()],2)) 
    return

def refresh_graph2(event):
    if len(cmbx_IRplt["values"])!=0:
        sgn = signals[cmbx_IRplt.current()][cmbx_bandplt.current()]
        env = envelopes[cmbx_IRplt.current()][cmbx_bandplt.current()]
        frecs = fs[cmbx_IRplt.current()]
        graph2.cla()
        graph2.plot(np.arange(0,len(sgn))/frecs,sgn)
        graph2.plot(np.arange(0,len(env))/frecs,env)
        graph2.set_xlim([0,len(sgn)/frecs])
        graph2.set_xlabel("Time [s]")
        graph2.set_ylabel("Level [dBFS]")
        graph2.grid()
        canv2.draw()
    return

RIRs=[]
fs=[]
bands_centers=[(31.5,63,125,250,500,1000,2000,4000,8000,16000),(25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000)]
bands_labels=[("31.5","63","125","250","500","1k","2k","4k","8k","16k"),("25","31.5","40","50","63","80","100","125","160","200","250","315","400","500","630","800","1k","1.25k","1.6k","2k","2.5k","3.15k","4k","5k","6.3k","8k","10k","12.5k","16k","20k")]
parameters=("EDT", "T20", "T30", "C50", "C80","Tt", "EDTt", "IACC")
rowsaverage=("Average","Max","Min","Sigma")

APC = Tk()
APC.title("Acoustic Parameters Calculator")
APC.rowconfigure(0,weight=1)
APC.columnconfigure(1,weight=1)
APC["padx"]="2"
APC["pady"]="2"
APC.geometry("1000x600")

figgraph1 = Figure(constrained_layout=True)
graph1 = figgraph1.add_subplot(111)
graph1.set_xlabel("Frecuency [Hz]")
graph1.set_ylabel("")
r , g , b = APC.winfo_rgb(APC["bg"])
figgraph1.set_facecolor([r/65536,g/65536,b/65536])

canv1 = FigureCanvasTkAgg(figgraph1, master = APC)
canv1.draw()
get_widz1 = canv1.get_tk_widget()
get_widz1.grid(row=0, column=1, padx=4, pady=4, sticky="ewns")
get_widz1.grid_remove()

figgraph2 = Figure(constrained_layout=True)
graph2 = figgraph2.add_subplot(111)
graph2.set_xlabel("Time [s]")
graph2.set_ylabel("Level [dBFS]")
graph2.set_ylim([-100,0])
r , g , b = APC.winfo_rgb(APC["bg"])
figgraph2.set_facecolor([r/65536,g/65536,b/65536])

canv2 = FigureCanvasTkAgg(figgraph2, master = APC)
canv2.draw()
get_widz2 = canv2.get_tk_widget()
get_widz2.grid(row=0, column=1, padx=4, pady=4, sticky="ewns")


frame_menu = Frame(APC)
frame_menu.grid(row=0, column=0, sticky="ewns", padx=2, pady=2)
frame_menu.rowconfigure(4,weight=1)


btn_load = Button(frame_menu, text="Load RIRs", command=func_load)
btn_load.grid(row=0, column=0, padx=2, pady=2, sticky="ewns")

btn_clearall = Button(frame_menu, text="Clear all", command=func_clearall)
btn_clearall.grid(row=0, column=1, padx=2, pady=2, sticky="ewns")

btn_clearselected = Button(frame_menu, text="Clear selected", command=func_clearselected)
btn_clearselected.grid(row=0, column=2, padx=2, pady=2, sticky="ewns")

btn_calculate = Button(frame_menu, text="Calculate", font=("Arial",12,"bold"),command=func_calculate)
btn_calculate.grid(row=5, column=0,columnspan=3, sticky="ewns", padx=2, pady=2)

btn_show = Button(frame_menu, text="Show average", command=func_show)
btn_show.grid(row=6, column=0, sticky="ewns", padx=2, pady=2)

btn_expplot = Button(frame_menu, text="Export plot")
btn_expplot.grid(row=6, column=1, sticky="ewns", padx=2, pady=2)

btn_exptable = Button(frame_menu, text="Export table")
btn_exptable.grid(row=6, column=2, sticky="ewns", padx=2, pady=2)

lbl_IRplt = Label(frame_menu,text="RIR plot:")
lbl_IRplt.grid(row=7, column=0, sticky="ens",pady=2)

cmbx_IRplt = Combobox(frame_menu)
cmbx_IRplt.grid(row=7, column=1, columnspan=2, sticky="wns", pady=2, padx=2)
cmbx_IRplt.bind("<<ComboboxSelected>>",refresh_table2)

lbl_bandplt = Label(frame_menu,text="Band:")
lbl_bandplt.grid(row=8, column=0, sticky="ens", pady=2)

cmbx_bandplt = Combobox(frame_menu)
cmbx_bandplt.grid(row=8, column=1, columnspan=2, sticky="wns", pady=2, padx=2)
cmbx_bandplt.bind("<<ComboboxSelected>>",refresh_graph2)

frame_IRs = Frame(frame_menu)
frame_IRs.grid(row=4, column=0, columnspan=3, sticky="ewns", padx=2, pady=2)
frame_IRs.columnconfigure(0, weight=1)
frame_IRs.rowconfigure(0, weight=1)
scrollbarV=Scrollbar(frame_IRs,orient="vertical")
scrollbarV.grid(row=0,column=1,sticky="ewns")
scrollbarH=Scrollbar(frame_IRs,orient="horizontal")
scrollbarH.grid(row=1,column=0,sticky="ewns")
lstbx_IRs = Listbox(frame_IRs, selectmode="extended")
lstbx_IRs.grid(row=0, column=0, sticky="ewns")
lstbx_IRs.config(yscrollcommand=scrollbarV.set,xscrollcommand=scrollbarH.set)
scrollbarV.config(command=lstbx_IRs.yview)
scrollbarH.config(command=lstbx_IRs.xview)

frame_settings = Frame(frame_menu)
frame_settings.grid(row=3, column=0, columnspan=3, sticky="ewns", padx=2, pady=2)

lbl_settings = Label(frame_settings, text="Settings", anchor="center")
lbl_settings.grid(row=0, column=0,columnspan=4, sticky="ewns")

lbl_bandfilter = Label(frame_settings, text="Band filter:", anchor="e")
lbl_bandfilter.grid(row=1, column=0, sticky="ewns")

cmbx_bandfilter = Combobox(frame_settings, values=("Octave","Third"), state="readonly", width=7)
cmbx_bandfilter.grid(row=1, column=1, sticky="ewns", pady=1)
cmbx_bandfilter.bind("<<ComboboxSelected>>",func_bandfilter)

lbl_minband = Label(frame_settings, text="Min band:", anchor="e")
lbl_minband.grid(row=2, column=0, sticky="ewns")

cmbx_minband = Combobox(frame_settings, state="readonly", width=7)
cmbx_minband.grid(row=2, column=1, sticky="ewns", pady=1)
cmbx_minband.bind("<<ComboboxSelected>>",func_limitband)

lbl_maxband = Label(frame_settings, text="Max band:", anchor="e")
lbl_maxband.grid(row=3, column=0, sticky="ewns")

cmbx_maxband = Combobox(frame_settings, state="readonly", width=7)
cmbx_maxband.grid(row=3, column=1, sticky="ewns", pady=1)
cmbx_maxband.bind("<<ComboboxSelected>>",func_limitband)

lbl_envelope = Label(frame_settings, text="Envelope:", anchor="e")
lbl_envelope.grid(row=1, column=2, sticky="ewns")

cmbx_envelope = Combobox(frame_settings, values=("Schoreder","Median MF","Mean MF"), state="readonly", width=10)
cmbx_envelope.current(0)
cmbx_envelope.grid(row=1, column=3, sticky="ewns", pady=1)

lbl_chunk = Label(frame_settings, text="IR chunk:", anchor="e")
lbl_chunk.grid(row=2, column=2, sticky="ewns")

cmbx_chunk = Combobox(frame_settings, values=("Lundeby", "Pepino", "K0"), state="readonly", width=10)
cmbx_chunk.current(0)
cmbx_chunk.grid(row=2, column=3, sticky="ewns", pady=1)

reversedcheck = IntVar()
chkbtn_reverse = Checkbutton(frame_settings, text="Reversed Filtering", variable=reversedcheck)
chkbtn_reverse.grid(row=3, column=2, columnspan=2, sticky="ewns")

frame_table1 = Frame(APC)
frame_table1.grid(row=1,column=0,columnspan=2,padx=4,pady=4,sticky="ewns")
frame_table1.grid_remove()
for i in range(len(rowsaverage)+1):
    for j in range(len(bands_labels[-1])+2):
        if i==0 and j>1:
            table1[i] = table1[i] + [Label(frame_table1, text=bands_labels[1][j-2], bg="white", relief="solid",borderwidth=1)]
            table1[i][j].grid(row=i,column=j,sticky="ewns")
            frame_table1.columnconfigure(j, weight=1)
        elif i==0 and j==1:
            table1[i] = table1[i] + [Label(frame_table1, text="global", bg="white", relief="solid",borderwidth=1)]
            table1[i][j].grid(row=i,column=j,sticky="ewns")
            frame_table1.columnconfigure(j, weight=1)
        
        elif i!=0 and j==0:
            table1 = table1 + [[Label(frame_table1,text=rowsaverage[i-1], bg="white", relief="solid",borderwidth=1)]]
            table1[i][j].grid(row=i,column=j,sticky="ewns")
        elif i!=0 and j!=0:
            table1[i] = table1[i] + [Label(frame_table1,text="", bg="white", relief="solid",borderwidth=1)]
            table1[i][j].grid(row=i,column=j,sticky="ewns")
        elif i==0 and j==0:
            table1=[[Combobox(frame_table1, values=parameters[0:-1], width=5, state="readonly")]]
            table1[i][j].current(0)
            table1[i][j].bind("<<ComboboxSelected>>",refresh_graphtable1)
            table1[i][j].grid(row=i, column=j, sticky="ewns")

frame_table2 = Frame(APC)
frame_table2.grid(row=1,column=0,columnspan=2,padx=4,pady=4,sticky="ewns")
for i in range(len(parameters)+1):
    for j in range(len(bands_labels[-1])+2):
        if i==0 and j>1:
            table2[i] = table2[i] + [Label(frame_table2, text=bands_labels[1][j-2], bg="white", relief="solid",borderwidth=1)]
            table2[i][j].grid(row=i,column=j,sticky="ewns")
            frame_table2.columnconfigure(j, weight=1)
        elif i==0 and j==1:
            table2[i] = table2[i] + [Label(frame_table2, text="global", bg="white", relief="solid",borderwidth=1)]
            table2[i][j].grid(row=i,column=j,sticky="ewns")
            frame_table2.columnconfigure(j, weight=1)
        elif i!=0 and j==0:
            table2 = table2 + [[Label(frame_table2,text=parameters[i-1], bg="white", relief="solid",borderwidth=1)]]
            table2[i][j].grid(row=i,column=j,sticky="ewns")
        elif i!=0 and j!=0:
            table2[i] = table2[i] + [Label(frame_table2,text="", bg="white", relief="solid",borderwidth=1)]
            table2[i][j].grid(row=i,column=j,sticky="ewns")
        elif i==0 and j==0:
            table2=[[""]]

APC.mainloop()