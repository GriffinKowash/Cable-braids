import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable braids\models\test_23-06-05\3D_Braid_model\2'
EFIELD_NAME = 'Field_Probe_1.dat'
CURRENT_NAME = 'Bulk_Current_1.dat'

smoothing = False
plot_all = True

efield_data = np.loadtxt(SOURCE_DIR + '\\' + EFIELD_NAME)
current_data = np.loadtxt(SOURCE_DIR + '\\' + CURRENT_NAME)

et, ex, ey, ez = efield_data.T
it, i = current_data.T

e = np.sqrt(ex**2 + ey**2 + ez**2)

ztrans = e / (i + 1e-12)

dt = et[1] - et[0]
freq = np.fft.rfftfreq(e.size) / dt
ztrans_fft = np.abs(np.fft.rfft(ztrans)) / it.size
i_fft = np.abs(np.fft.rfft(i)) / it.size
e_fft = np.abs(np.fft.rfft(e)) / it.size

#plt.plot(et, e)
#plt.show()

output = np.array(list(zip(freq, ztrans_fft)))
np.savetxt(SOURCE_DIR + '\\' + 'Transfer_Impedance_fft.dat', output)

sns.set()

if plot_all:
    fig, ax = plt.subplots(2, 2)
    
    ax[0, 0].semilogx(it, np.abs(ztrans))
    ax[0, 0].set_yscale('log')
    ax[0, 0].set_ylim(1e-2, 1e15)
    ax[0, 0].set_title('Transfer impedance magnitude time series')
    ax[0, 0].set_xlabel('Time (s)')
    ax[0, 0].set_ylabel('Transfer impedance magnitude (Ω/m)')
    
    if smoothing:
        n = 20
        ztrans_fft = np.convolve(ztrans_fft, np.ones(n), mode='valid') / n
        start = int(n / 2)
        end = freq.size - start + 1
        ax[0, 1].semilogx(freq[start:end], ztrans_fft)
        ax[0, 1].set_title(f'Transfer impedance fft (smoothed, n={n})')
    else:
        ax[0, 1].semilogx(freq, ztrans_fft)
        ax[0, 1].set_title('Transfer impedance fft')
    ax[0, 1].set_xlim(1e6, 5e8)
    ax[0, 1].set_xlabel('Frequency (Hz)')
    ax[0, 1].set_ylabel('Transfer impedance (Ω/m)')
    
    ax[1, 0].semilogx(freq, i_fft)
    ax[1, 0].set_title('Bulk current fft')
    ax[1, 0].set_xlim(1e6, 5e8)
    ax[1, 0].set_ylim(0, 50)
    ax[1, 0].set_xlabel('Frequency (Hz)')
    ax[1, 0].set_ylabel('Current (A)')
    
    ax[1, 1].semilogx(freq, e_fft)
    ax[1, 1].set_title('E-field fft')
    ax[1, 1].set_xlim(1e6, 5e8)
    ax[1, 1].set_xlabel('Frequency (Hz)')
    ax[1, 1].set_ylabel('Electric field (V/m)')
    ax[1, 1].set_ylim(0, 600)
    #ax.set_yscale('log')
    #ax.semilogx(freq, i_fft)
    #ax.semilogx(freq, ztrans_fft)
    #ax.set_ylim((-10, 2e9))
    fig.suptitle('New braid')
    
else:
    fig, ax = plt.subplots(1)
    if smoothing:
        n = 20
        ztrans_fft = np.convolve(ztrans_fft, np.ones(n), mode='valid') / n
        start = int(n / 2)
        end = freq.size - start + 1
        ax.semilogx(freq[start:end], ztrans_fft)
        ax.set_title(f'Transfer impedance fft (smoothed, n={n})')
    else:
        ax.semilogx(freq, ztrans_fft)
        ax.set_title('Transfer impedance fft')
        
    ax.set_xlim((1e6, 5e8))
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Transfer impedance (Ω/m)')

fig.suptitle('Old braid')
fig.show()