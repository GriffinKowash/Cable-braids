import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\Yazaki\original_100MHz'
SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\new14b_thinwire\1GHz'
#SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\Vance_k85_surface\10GHz'
#SOURCE_DIR2 = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\Vance_k85_surface\10GHz_large-domain_PEC'
#SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\old\Yazaki_validation\Yazaki_cable_EMA3D\original_div-16'
#SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\old\Yazaki_validation\Yazaki_cable_EMA3D_disco\original_div-16'
#SOURCE_DIR2 = r'C:\Users\griffin.kowash\Documents\Projects\Cable braids\models\Griffin_cable_EMA3D_disco\1MHz'
SOURCE_DIR2 = None

EFIELD_NAME = 'Field_Probe.dat'
EFIELD_NAME2 = 'Field_Probe.dat'
CURRENT_NAME = 'Full_Cable.dat'
CURRENT_NAME2 = 'Full_Cable.dat'
TARGET_FREQ = 1e9
PERIOD = 1 / TARGET_FREQ

cumulative = True  #not supported for multiple source files
smoothing = False
plot_all = True
interpolate = True

cutoff = None#5e-6

efield_data = np.loadtxt(SOURCE_DIR + '\\' + EFIELD_NAME)
current_data = np.loadtxt(SOURCE_DIR + '\\' + CURRENT_NAME)

if SOURCE_DIR2 != None:
    efield_data2 = np.loadtxt(SOURCE_DIR2 + '\\' + EFIELD_NAME2)
    current_data2 = np.loadtxt(SOURCE_DIR2 + '\\' + CURRENT_NAME2)
    
if cumulative:
    if cutoff != None:
        arg = np.argmin(np.abs(efield_data[:, 0] - cutoff))
        efield_data = efield_data[:arg, :]
        current_data = current_data[:arg, :]
        
else:
    arg = np.argmin(np.abs(efield_data[:, 0] - (efield_data[-1, 0] - PERIOD)))
    efield_data = efield_data[arg:, :]
    current_data = current_data[arg:, :]
    

et, ex, ey, ez = efield_data.T
it, i = current_data.T

e = np.sqrt(ex**2 + ey**2 + ez**2)


dt = 3.4e-13 #et[1] - et[0]
efreq = np.fft.rfftfreq(efield_data.shape[0]) / dt
ifreq = np.fft.rfftfreq(current_data.shape[0]) / dt
i_fft = np.abs(np.fft.rfft(i, norm='forward')) * 2  #sketchy factor of 2, double check (apparently it's legit)
e_fft = np.abs(np.fft.rfft(ez, norm='forward')) * 2

if e_fft.size == i_fft.size:
    ztrans_fft = e_fft / i_fft
    target_index = np.argmin(np.abs(efreq - TARGET_FREQ))
    
    if interpolate:
        if efreq[target_index] < TARGET_FREQ:
            zi = e_fft[target_index] / i_fft[target_index]
            zf = e_fft[target_index + 1] / i_fft[target_index + 1]
            fi = efreq[target_index]
            ff = efreq[target_index + 1]
            dzdf = (zf - zi) / (ff - fi)
            z0 = zi + (TARGET_FREQ - fi) * dzdf
        
        elif efreq[target_index] > TARGET_FREQ:
            zi = e_fft[target_index - 1] / i_fft[target_index - 1]
            zf = e_fft[target_index] / i_fft[target_index]
            fi = efreq[target_index - 1]
            ff = efreq[target_index]
            dzdf = (zf - zi) / (ff - fi)
            z0 = zi + (TARGET_FREQ - fi) * dzdf
        
        else:
            z0 = e_fft[target_index] / i_fft[target_index]
    
    else:
        z0 = ztrans_fft[target_index]

else:
    ztrans_fft = None
    e_index = np.argmin(np.abs(efreq - TARGET_FREQ))
    i_index = np.argmin(np.abs(ifreq - TARGET_FREQ))
    
    if interpolate:
        if efreq[e_index] < TARGET_FREQ:
            zi = e_fft[e_index] / i_fft[i_index]
            zf = e_fft[e_index + 1] / i_fft[i_index + 1]
            fi = efreq[e_index]
            ff = efreq[e_index + 1]
            dzdf = (zf - zi) / (ff - fi)
            z0 = zi + (TARGET_FREQ - fi) * dzdf
        
        elif efreq[e_index] > TARGET_FREQ:
            zi = e_fft[e_index - 1] / i_fft[i_index - 1]
            zf = e_fft[e_index] / i_fft[i_index]
            fi = efreq[e_index - 1]
            ff = efreq[e_index]
            dzdf = (zf - zi) / (ff - fi)
            z0 = zi + (TARGET_FREQ - fi) * dzdf
        
        else:
            z0 = e_fft[e_index] / i_fft[i_index]
    
    else:
        z0 = e_fft[e_index] / i_fft[i_index]


#plt.plot(et, e)
#plt.show()

#output = np.array(list(zip(freq, ztrans_fft)))
#np.savetxt(SOURCE_DIR + '\\' + 'Transfer_Impedance_fft.dat', output)

if SOURCE_DIR2 != None:
    if cutoff != None:
        arg2 = np.argmin(np.abs(efield_data2[:, 0] - cutoff))
        efield_data2 = efield_data2[:arg2, :]
        current_data2 = current_data2[:arg2, :]

    et2, ex2, ey2, ez2 = efield_data2.T
    it2, i2 = current_data2.T

    e2 = np.sqrt(ex2**2 + ey2**2 + ez2**2)


    dt2 = et2[1] - et2[0]
    freq2 = np.fft.rfftfreq(e2.size) / dt2
    i_fft2 = np.abs(np.fft.rfft(i2, norm='forward')) * 2
    e_fft2 = np.abs(np.fft.rfft(ez2, norm='forward')) * 2

    ztrans_fft2 = e_fft2 / i_fft2

if SOURCE_DIR2 != None:
    target_index2 = np.argmin(np.abs(freq2 - TARGET_FREQ))
    z = ztrans_fft2[target_index2]
    error = abs(z - z0) / z0
    print(f'Transfer impedance at target frequency:\n\t{round(z0, 4)} Ω/m (reference)\n\t{round(z, 4)} Ω/m (test)\nPercent error: {round(error * 100, 2)}')

else:
    print(f'Transfer impedance at target frequency:\n\t{round(z0, 4)} Ω/m')  


sns.set()

if plot_all:
    fig, ax = plt.subplots(2, 2)
    
    ax[0, 0].plot(it, i)
    if SOURCE_DIR2 != None:
        ax[0, 0].plot(it2, i2)
    #ax[0, 0].set_yscale('log')
    #ax[0, 0].set_ylim(1e-2, 1e4)
    ax[0, 0].set_title('Bulk current timeseries')
    ax[0, 0].set_xlabel('Time (s)')
    ax[0, 0].set_ylabel('Current (A)')
    
    if ztrans_fft is not None:
        if smoothing:
            n = 20
            ztrans_fft = np.convolve(ztrans_fft, np.ones(n), mode='valid') / n
            start = int(n / 2)
            end = efreq.size - start + 1
            ax[0, 1].semilogx(efreq[start:end], ztrans_fft, label='reference')
            ax[0, 1].set_title(f'Transfer impedance fft (smoothed, n={n})')
        else:
            ax[0, 1].semilogx(efreq, ztrans_fft, label='reference')
            if SOURCE_DIR2 != None:
                ax[0, 1].semilogx(freq2, ztrans_fft2, label='test')
                ax[0, 1].set_title(f'Transfer impedance fft (ref={round(z0,4)}, test={round(z,4)}; {round(100 * error, 2)}% error at {round(TARGET_FREQ / 1e6)}MHz)')
            else:
                ax[0, 1].set_title('Transfer impedance fft')
        ax[0, 1].set_xlim(TARGET_FREQ * 0.5, TARGET_FREQ * 2)
        ax[0, 1].set_ylim(-1e-3, 1)
        ax[0, 1].set_xlabel('Frequency (Hz)')
        ax[0, 1].set_ylabel('Transfer impedance (Ω/m)')
    
    #ax[1, 0].semilogx(freq, i_fft)
    ax[1, 0].plot(et, ez)
    if SOURCE_DIR2 != None:
        #ax[1, 0].semilogx(freq2, i_fft2)
        ax[1, 0].plot(et2, ez2)
    #ax[1, 0].set_title('Bulk current fft')
    ax[1, 0].set_title('E-field (z) time series')
    #ax[1, 0].set_xlim(TARGET_FREQ * 0.5, TARGET_FREQ * 2)
    #ax[1, 0].set_ylim(-2, 50)
#    ax[1, 0].set_xlabel('Frequency (Hz)')
#    ax[1, 0].set_ylabel('Current (A)')
    ax[1, 0].set_xlabel('Time (s)')
    ax[1, 0].set_ylabel('Electric field (V/m)')
    
    ax[1, 1].semilogx(efreq, e_fft)
    if SOURCE_DIR2 != None:
        ax[1, 1].semilogx(freq2, e_fft2)
    ax[1, 1].set_title('E-field fft')
    ax[1, 1].set_xlim(TARGET_FREQ * 0.5, TARGET_FREQ * 2)
    ax[1, 1].set_xlabel('Frequency (Hz)')
    ax[1, 1].set_ylabel('Electric field (V/m)')
    #ax[1, 1].set_ylim(-2, 5)
    #ax.set_yscale('log')
    #ax.semilogx(freq, i_fft)
    #ax.semilogx(freq, ztrans_fft)
    #ax.set_ylim((-10, 2e9))
    #fig.suptitle('New braid')
    
    if SOURCE_DIR2 != None:
        ax[0,1].legend()

else:
    fig, ax = plt.subplots(1)
    if smoothing:
        n = 50
        ztrans_fft = np.convolve(ztrans_fft, np.ones(n), mode='valid') / n
        start = int(n / 2)
        end = efreq.size - start + 1
        ax.semilogx(efreq[start:end], ztrans_fft)
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

fig.suptitle('Braid transfer impedance comparison - 10MHz sine wave current')
fig.show()