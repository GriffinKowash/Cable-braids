import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Transfer_Impedance_Analytical_Model import calc_braid_impedance


SOURCE_DIR = r'C:\Users\griffin.kowash\Documents\Projects\Cable_braids\models\new14b_thinwire'
EFIELD_NAME = 'Field_Probe.dat'
CURRENT_NAME = 'Full_Cable.dat'
#FREQS = {'1MHz': 1e6, '2500kHz': 2.5e6, '10MHz': 1e7, '20MHz': 2e7, '30MHz': 3e7, '50MHz': 5e7, '100MHz': 1e8}
FREQS = {'10MHz': 1e7, '20MHz': 2e7, '30MHz': 3e7, '50MHz': 5e7, '100MHz': 1e8, '300MHz': 3e8, '1GHz': 1e9}
#FREQS = {'1GHz': 1e9, '10GHz': 1e10}
#FREQS = {'100MHz': 1e8}

#KOWASH = None
KOWASH = {1e6: 0.0734, 3e6: 0.0744, 1e7: 0.0201, 2e7: 0.0189, 3e7: 0.0212, 5e7: 0.0284, 1e8: 0.0503, 3e8: 0.1458, 1e9: 0.4823}
KOWASH_SURFACE = {5e7: 0.0380, 1e8: 0.0842, 3e8: 0.2450, 1e9: 0.6986}
#KOWASH_STD = {1e10: 4.62}
#KOWASH = {1e6: 0.6311, 2.5e6: 0.14, 1e7: 0.5526, 2e7: 0.98, 3e7: 1.6, 5e7: 2.8, 1e8: 5.3}
#VANCE_CODY = {1e4: 0.11, 1e5: 0.11, 1e6: 0.11, 2.5e6: 0.14, 1e7: 0.48, 2e7: 0.98, 3e7: 1.6, 5e7: 2.8, 1e8: 5.3}
#VANCE = {1e4: 0.02803, 1e5: 0.02803, 1e6: 0.02403, 2.5e6: 0.01309, 1e7: 0.01807, 2e7: 0.03128, 3e7: 0.04661, 5e7: 0.07794, 1e8: 0.15585}
YAZAKI = {1e4: 0.1, 1e5: 0.1, 1e6: 0.1, 1e7: 0.25, 2e7: 0.50, 3e7: 0.75, 5e7: 1.25, 1e8: 2.5}
CUTOFF = None
MODEL = 'Vance'

def calc_zt(e_data, i_data, freq, cumulative=True):
    if not cumulative:
        period = 1 / freq
        tstart = e_data[-1, 0] - period
        istart = np.argmin(np.abs(e_data[:, 0] - tstart))
        e_data = e_data[istart:, :]
        i_data = i_data[istart:, :]
    
    et, ex, ey, ez = e_data.T
    it, i = i_data.T
    e = np.sqrt(ex**2 + ey**2 + ez**2)
    
    edt = np.mean(et[1:] - et[:-1])
    idt = np.mean(it[1:] - it[:-1])
    
    e_freqs = np.fft.rfftfreq(et.size) / edt
    i_freqs = np.fft.rfftfreq(it.size) / idt
    
    e_fft = np.abs(np.fft.rfft(ez, norm='forward')) * 2
    i_fft = np.abs(np.fft.rfft(i, norm='forward')) * 2  #sketchy factor of 2, double check
    
    e_index = np.argmin(np.abs(e_freqs - freq))
    i_index = np.argmin(np.abs(i_freqs - freq))
    
    zt = e_fft[e_index] / i_fft[i_index]
    
    return zt


def load_results(path, e_name='Field_Probe.dat', i_name='Full_Cable.dat'):
    e_data = np.loadtxt(path + '\\' + e_name)
    i_data = np.loadtxt(path + '\\' + i_name)
    
    return e_data, i_data


def trim_results(data, cutoff):
    if cutoff != None:
        arg = np.argmin(np.abs(data[:, 0] - cutoff))
        data = data[:arg, :]
        
    return data


def calc_error(zts, val):
    error = {}
    for freq, zt in zts.items():
        error[freq] = (zt - val[freq]) / val[freq]
        
    return error


def print_zts(zts, val=None):
    if val != None:
        errors = calc_error(zts, val)
        for freq, zt in zts.items():
            print(f'{round(freq / 1e6, 1)} MHz: {round(zt, 4)} Ω/m ({round(errors[freq] * 100, 1)}% error)')
    else:
        for freq, zt in zts.items():
            print(f'{round(freq / 1e6, 1)} MHz: {round(zt, 4)} Ω/m')


def plot(zts, val=None, ref=None):
    sns.set()
    fig, ax = plt.subplots(1)
    fig.suptitle('Cable braid transfer impedance, simulated')
    
    if val != None:
        ax.plot(val.keys(), val.values(), label=MODEL, color='C0')
        fig.suptitle('Cable braid transfer impedance, simulated versus analytical')
        
    if ref != None:
        #ax.plot(ref.keys(), ref.values(), linestyle='--', color='C2')
        ax.scatter(ref.keys(), ref.values(), label='Reference', color='C2')
        fig.suptitle('Cable braid transfer impedance, simulated versus analytical')
    
    #ax.plot(zts.keys(), zts.values(), linestyle='--', color='C1')
    ax.plot(zts.keys(), zts.values(), linestyle='-.', color='C1')
    ax.scatter(zts.keys(), zts.values(), label='Test', color='C1')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Transfer impedance (Ω/m)')
    ax.legend()

    fig.tight_layout()
    fig.show()
    
if KOWASH == None:
    KOWASH = {}
    for name, freq in FREQS.items():
        path = SOURCE_DIR + '\\' + name
        e_data, i_data = load_results(path)
        e_data, i_data = trim_results(e_data, CUTOFF), trim_results(i_data, CUTOFF)    
        
        zt = calc_zt(e_data, i_data, freq, cumulative=True)
        KOWASH[freq] = zt

dataf, ztv = calc_braid_impedance(MODEL)
val = dict(zip(dataf, ztv))

plot(KOWASH_SURFACE, val=val, ref=KOWASH)
print_zts(KOWASH, val=val)