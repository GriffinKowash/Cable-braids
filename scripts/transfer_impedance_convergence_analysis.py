import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\Vance_k85_surface\10GHz_long_finemesh'
SOURCE_DIR2 = None#r'C:\Users\griffin.kowash\Documents\Projects\Cable braids\models\old\Yazaki_validation\Yazaki_cable_EMA3D\original_div-16'
EFIELD_NAME = 'Field_Probe.dat'
EFIELD_NAME2 = 'Field_Probe.dat'
CURRENT_NAME = 'Full_Cable.dat'
CURRENT_NAME2 = 'Full_Cable.dat'
TARGET_FREQ = 1e10
CUMULATIVE = True

smoothing = False
plot_all = True


efield_data = np.loadtxt(SOURCE_DIR + '\\' + EFIELD_NAME)
current_data = np.loadtxt(SOURCE_DIR + '\\' + CURRENT_NAME)

if SOURCE_DIR2 != None:
    efield_data2 = np.loadtxt(SOURCE_DIR2 + '\\' + EFIELD_NAME2)
    current_data2 = np.loadtxt(SOURCE_DIR2 + '\\' + CURRENT_NAME2)
    
step = 1 / TARGET_FREQ
end = 5e-9
cutoffs = np.arange(step, end, step)
zts  = []
zts2 = []

for cutoff in cutoffs:
    if CUMULATIVE:
        arg = np.argmin(np.abs(efield_data[:, 0] - cutoff))
        efield_slice = efield_data[:arg, :]
        current_slice = current_data[:arg, :]
        
    else:
        arg_start = np.argmin(np.abs(efield_data[:, 0] - (cutoff - step)))
        arg_end = np.argmin(np.abs(efield_data[:, 0] - cutoff))
        efield_slice = efield_data[arg_start:arg_end, :]
        current_slice = current_data[arg_start:arg_end, :]
    
    et, ex, ey, ez = efield_slice.T
    it, i = current_slice.T
    
    e = np.sqrt(ex**2 + ey**2 + ez**2)
    
    
    dt = et[1] - et[0]
    freq = np.fft.rfftfreq(e.size) / dt
    i_fft = np.abs(np.fft.rfft(i, norm='forward'))
    e_fft = np.abs(np.fft.rfft(ez, norm='forward'))
    ztrans_fft = e_fft / i_fft
    
    target_index = np.argmin(np.abs(freq - TARGET_FREQ))
    z0 = ztrans_fft[target_index]
    zts.append(z0)

    if SOURCE_DIR2 != None:
        if CUMULATIVE:
            arg2 = np.argmin(np.abs(efield_data2[:, 0] - cutoff))
            efield_slice2 = efield_data2[:arg2, :]
            current_slice2 = current_data2[:arg2, :]
            
        else:
            arg_start2 = np.argmin(np.abs(efield_data2[:, 0] - (cutoff - step)))
            arg_end2 = np.argmin(np.abs(efield_data2[:, 0] - cutoff))
            efield_data2 = efield_data2[arg_start2:arg_end2, :]
            current_data2 = current_data2[arg_start2:arg_end2, :]
    
        et2, ex2, ey2, ez2 = efield_data2.T
        it2, i2 = current_data2.T
    
        e2 = np.sqrt(ex2**2 + ey2**2 + ez2**2)
    
        dt2 = et2[1] - et2[0]
        freq2 = np.fft.rfftfreq(e2.size) / dt2
        i_fft2 = np.abs(np.fft.rfft(i2, norm='forward'))
        e_fft2 = np.abs(np.fft.rfft(ez2, norm='forward'))
    
        ztrans_fft2 = e_fft2 / i_fft2
    

    
    if SOURCE_DIR2 != None:
        target_index2 = np.argmin(np.abs(freq2 - TARGET_FREQ))
        z = ztrans_fft2[target_index2]
        error = abs(z - z0) / z0
        #print(f'Transfer impedance at target frequency:\n\t{round(z0, 4)} Ω/m (reference)\n\t{round(z, 4)} Ω/m (test)\nPercent error: {round(error * 100, 2)}')
        zts2.append(z)

    

plt.plot(cutoffs, zts, color='C0', label='reference')
plt.scatter(cutoffs, zts, color='C0')

if SOURCE_DIR2 != None:
    plt.plot(cutoffs, zts2, color='C1', label='test')
    plt.scatter(cutoffs, zts2, color='C1')

plt.xlabel('End time (s)')
plt.ylabel('Transfer impedance (Ω/m)')
plt.title('Transfer impedance versus simulation end time')
plt.legend()

plt.show()