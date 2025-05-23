#!/usr/bin/env python3
import numpy as np
from scipy import signal
import math, os, sys, time
import matplotlib.pyplot as plt


##### PLEASE READ THE FOLLOWING INSTRUCTIONs BEFORE RUNNING SCRIPT #####
####                                                                ####
####  The Format for Running This Script:                           ####
####  python IR_total_KW.py INPUT_FILE DELTA_T WINDOW OUTPUT_FILE   ####
####                                                                ####
####  The values need to input manually when runing this script     ####
####                                                                #### 
####  (1) INPUT_FILE_NAME: The Total_Dipole_Moment_*.Diople file    ####
####           (NOTE: do NOT need to re-split the Dipole file)      ####
####                                                                #### 
####  (2) DELTA_T: The Time_step set in simulation, in unit of fs   ####
####                                                                ####
####  (3) WINDOW: The Name of the Window Function                   ####
####                                                                ####  
####  (4) OUTPUT_FILE_NAME: The Name of the Output File.            ####
####           (NOTE: do NOT need to type > sign!)                  ####
####                                                                ####
#############################  Let's Try It! ###########################






#### The values need to input manually when running this script ####
#path = input("\nPlease Enter the Directory Contained the Dipole Moment File\n")
fname = sys.argv[1]                    # The name of the input file
delta_t = float(sys.argv[2]) * 1.0e-15 # The time step in unit of femtoseconds
window = sys.argv[3]                   # The name of the window function
V = float(sys.argv[4])                 # Unit cell volume in A^3
T = float(sys.argv[5])		       # K
fout = sys.argv[6]                     # The name of the output file



#### The constants will be used in this script ####
EPSILON0=8.854187817*1.0e-12       # Vacuum permittivity.  C^2 N^-1 m^-2
kB = 1.38064852*1.0e-23            # m^-1 N^-1 K-1
beta = 1.0/(kB * T)                # m^-1 N^-1
hbar = 1.0545718*1.0e-34           # m^2 kg / s = N * s * m
betahbar = beta*hbar               # s
#V = np.load('box.npy')
V *=  1.0e-30     #  m^3
c = 299792458     # m/s
PI = np.pi
# 1 e = 1.602*1.0e-19 C
# change unit to C*m for M(0)
unit_basic = 1.602176565*1.0e-19*1.0e-10;       # C m
# change unit to s for dM(0)/dt
unit = unit_basic                               # C m s^-1
unit2 = unit * unit                             # C^2 m^2 s^-2
unit_all = 2.0*PI*beta/3.0/c/V*unit2            # m^-3 N^-1 * s^-1 * C^2

# 1F = A^2 * s^4 * kg^-1 * m^-2
pref = unit_all/(4*PI*EPSILON0);                # m^-1 * s^-1, s^-1 taken by dt.
pref /= 100     # change m^-1 to cm^-1


#### Functions will used in this script ####

def read_data(fname,dt):
  f=np.loadtxt(fname).T
  autocorr=f[1]
  crosscorr=f[2]
  return autocorr,crosscorr

def choose_window(data, kind='string',delta_t = 1e-15):
    ii = int(1e-12 / delta_t) * 0.5
    if kind == 'Gaussian':
        sigma = 2 * math.sqrt(2 * math.log(2))
        window = signal.gaussian(len(data)*2, std=ii/sigma, sym=False)[len(data):]
    elif kind == 'BH':
        window = signal.blackmanharris(len(data), sym=False)
    elif kind == 'Hamming':
        window = signal.hamming(len(data), sym=False)
    elif kind == 'Hann':
        window = signal.hann(len(data), sym=False)
    return window

def zero_padding(sample_data):
    '''
      A series of Zeros will be padded to the end of the dipole moment array 
    (before FFT performed), in order to obtain a array with the length which
    is the "next power of two" of numbers.
    #### Next power of two is calculated as: 2**np.ceil(log2(x))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    '''
    N = 2**int(math.log(len(sample_data)*2-1, 2))
    return N

def calc_FFT(sig, window,delta_t):
    '''
    This function is for calculating the "intensity" of the ACF at each 
    frequency by using the discrete fast Fourier transform.
    
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    '''
    # A series of number of zeros will be padded to the end of the DACF \
    # array before FFT.
    sig-=np.mean(sig[-10:])
    N = zero_padding(sig)
	
    yfft = np.fft.fft(sig, N, axis=0) # / len(sig)
    l=yfft.shape[0]
    wavenumber = np.fft.fftfreq(l, delta_t*c*100)[0:int(l/2)]
    realsig=np.real(yfft[0:int(l/2)])
    return (realsig-np.mean(realsig[-100:]))*wavenumber**2, wavenumber

def calc_FFT2(sig, window, delta_t):
    '''
    This function is for calculating the "intensity" of the ACF at each 
    frequency by using the discrete fast Fourier transform.
    
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    '''
    # A series of number of zeros will be padded to the end of the DACF \
    # array before FFT.
    LL = 5000
    c = 299792458 
    wavenumber = np.array(range(LL))
    yfft = np.array(range(LL)) * 0.0
    tt = np.array(range(len(sig))) * delta_t
    for ii in range(LL):
        yfft[ii] = np.sum(2 * np.cos(wavenumber[ii] * 2*np.pi / 0.01 * c * tt) * sig[ii])

    return yfft, wavenumber

######## Save The Results to A TEXT File ########
def save_results(fout, wavenumber, intensity_auto, intensity_cross):
    title = ("#Wavenumber", "IR Intensity (autocorr)", "IR Intensity (cross)", "cm^-1", "a.u.","a.u.")
    np.save(fout,np.array([wavenumber, intensity_auto, intensity_cross]))

######## Plot The Spectrum by Using Matplotlib module ########
def visualization(wavenumber, intensity, delta_t, label):
    plt.plot(wavenumber, intensity, color='black',  linewidth=1.5, label=label)
    plt.axis([0, 6000,
             -1.1*np.min(intensity), 1.1*np.max(intensity)])
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
    plt.subplots_adjust(hspace = 0.5)
    plt.legend()
    plt.show()

######## The main program ########
def main(fname, delta_t, window, fout):
    autocorr,crosscorr = read_data(fname,delta_t) 
    yfft_auto,wavenumber = calc_FFT(autocorr, window, delta_t)
    yfft_cross,wavenumber = calc_FFT(crosscorr, window, delta_t)
    yfft_auto *= delta_t 
    yfft_cross *= delta_t 
    intensity_a = pref * yfft_auto #/ 100
    visualization(wavenumber,intensity_a,delta_t,'Intra')
    intensity_c = pref * yfft_cross #/ 100
    visualization(wavenumber,intensity_c,delta_t,'Inter')
    save_results(fout, wavenumber, 2 * intensity_a, 2 * intensity_c)

if __name__ == '__main__':
    main(fname, delta_t, window, fout)
