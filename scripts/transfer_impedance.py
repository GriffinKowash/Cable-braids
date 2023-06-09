import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable braids\models\test_23-06-08_1MHz_sine\new_braid_sine_test\1'
EFIELD_NAME = 'Field_Probe_1.dat'
CURRENT_NAME = 'Bulk_Current_1.dat'

smoothing = False
plot_all = True

efield_data = np.loadtxt(SOURCE_DIR + '\\' + EFIELD_NAME)
current_data = np.loadtxt(SOURCE_DIR + '\\' + CURRENT_NAME)

et, ex, ey, ez = efield_data.T
it, i = current_data.T

e = np.sqrt(ex**2 + ey**2 + ez**2)

limit = min(e.size, i.size)
#e = e[:limit]
#i = i[:limit]

dt = et[1] - et[0]
freq = np.fft.rfftfreq(e.size) / dt
i_fft = np.abs(np.fft.rfft(i)) / it.size
e_fft = np.abs(np.fft.rfft(e)) / it.size

ztrans_fft = e_fft / i_fft

#plt.plot(et, e)
#plt.show()

output = np.array(list(zip(freq, ztrans_fft)))
np.savetxt(SOURCE_DIR + '\\' + 'Transfer_Impedance_fft.dat', output)

sns.set()

if plot_all:
    fig, ax = plt.subplots(2, 2)
    
    ax[0, 0].plot(it, i)
    #ax[0, 0].set_yscale('log')
    #ax[0, 0].set_ylim(1e-2, 1e4)
    ax[0, 0].set_title('Bulk current timeseries')
    ax[0, 0].set_xlabel('Time (s)')
    ax[0, 0].set_ylabel('Current (A)')
    
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
    #ax[1, 0].set_ylim(-2, 50)
    ax[1, 0].set_xlabel('Frequency (Hz)')
    ax[1, 0].set_ylabel('Current (A)')
    
    ax[1, 1].semilogx(freq, e_fft)
    ax[1, 1].set_title('E-field fft')
    ax[1, 1].set_xlim(1e6, 5e8)
    ax[1, 1].set_xlabel('Frequency (Hz)')
    ax[1, 1].set_ylabel('Electric field (V/m)')
    #ax[1, 1].set_ylim(-2, 5)
    #ax.set_yscale('log')
    #ax.semilogx(freq, i_fft)
    #ax.semilogx(freq, ztrans_fft)
    #ax.set_ylim((-10, 2e9))
    fig.suptitle('New braid')
    
else:
    fig, ax = plt.subplots(1)
    if smoothing:
        n = 50
        ztrans_fft = np.convolve(ztrans_fft, np.ones(n), mode='valid') / n
        start = int(n / 2)
        end = freq.size - start + 1
        ax.semilogx(freq[start:end], ztrans_fft)
        ax.set_title(f'Transfer impedance fft (smoothed, n={n})')
    else:
        ax.plot(et, ex, label='E_x')
        ax.plot(et, ey, label='E_y')
        ax.plot(et, ez, label='E_z')
        ax.set_title('Electric field components time series')
        ax.legend()
        
    #ax.set_xlim((1e6, 5e8))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('E-field (V/m)')

fig.suptitle('New braid - 10MHz sine wave current')
fig.show()