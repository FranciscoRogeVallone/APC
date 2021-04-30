import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt

fs=192000
tr=2

s=np.random.randn(fs*3)*(10**(-60*np.arange(fs*3)/fs/tr/20))
ir=np.random.randn(fs*4)/1000
ir[fs*1:]=ir[fs*1:]+s
ir=ir/np.max(np.abs(ir))

sf.write(r"C:\Users\Gabriel\Desktop\ir_synthetic_2seg.wav",ir,fs)
