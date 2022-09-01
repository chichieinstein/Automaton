import matlab.engine as engine
import numpy as np
import scipy.io as sio
from SSCA_Srivatsan_numba_version import *
import os
import shutil

class Analyzer:

    '''
    This class performs the following actions:

    1. Creates a signal by calling the appropriate atomic signal generator. 
    Radhika's script producing baseband non-resampled signals is called by the class.

    2. Channelises the signal by calling the channogram function.

    3. It also computes and stores the SSCA. 

    4. It plots the SSCA output in K-Q domain, as well as the one dimensional plots vs cycle frequency
    '''

    def __init__(self, snr_list, channel_list, signal_name, N, Np):

        folder_path = '/home/srajagopal/ssca_analysis' + signal_name
        if os.path.exists(folder_path):
            shutil.rmtree(folder_path, ignore_errors=True)

        os.makedirs(folder_path)
        eng = engine.start_matlab()

        self.snr_list = snr_list
        self.channel_list = channel_list
        self.signal_name = signal_name

        s = eng.genpath('/home/srajagopal/signalGenerator') # Srivatsan has the channelizer inside signalGenerator, so one command sufficed.
        eng.addpath(s, nargout=0)
        # If they are different folders, repeat these two commands with the path to the channelizer.

        arr, sample_rate = eng.python_integrator(signal_name, nargout=2)
        Obj_1 = CycloGram(100.0, 200.0, 8192, 128, False) # Computes the non conjugate scd
        Obj_2 = CycloGram(100.0, 200.0, 8192, 128, True) # Computes the conjugate scd

        for snr in snr_list:

            noisy_iq = eng.noise_adder(arr, snr)

            for Nch in channel_list:
            
                [s, f, t] = eng.channogram(noisy_iq, Nch, sample_rate, 16, 1, nargout=3)
            
                s = np.asarray(s)

                for ind, channel in enumerate(s):

                    ssca_output = Obj_1.cyclogram(channel, avg=True)

                    conj_ssca_output = Obj_2.cyclogram(channel, avg=True)

                    sio.savemat(folder_path + '/' + str(snr) + ' ' + str(ind) + '_non_conjugate_ssca.mat', {'arr' : ssca_output, 'downsampled_rate' : 2 * sample_rate / Nch})
                    sio.savemat(folder_path + '/' + str(snr) + ' ' + str(ind) + '_conjugate_ssca.mat', {'arr' : conj_ssca_output, 'downsampled_rate' : 2 * sample_rate / Nch})

    

    
    



