# This is a Python script for performing the SSCA of a given discrete time signal.

# References : 

# 1. "On the implementation of the SSCA for cyclic spectrum estimation by April, in Defence Research Establishment Ottawa, Feb. 1994 ".

# 2. "Computationally Efficient Algorithms for Cyclic Spectral Analysis by Roberts, Brown and Loomis, in IEEE SP Magazine, Apr. 1991".

from xmlrpc.client import FastParser
import numpy as np 

import pyfftw as ftw

import multiprocessing as mp

import pdb

import math

from collections import defaultdict

import scipy.signal as sig

from numba import njit, prange


@njit(parallel=True)
def create_matrix(record, N, Np):

	'''

		The create_matrix function takes the data record x[n] and constructs an N x Np matrix as follows. Note that for this to work

		the length of the record must be N+Np

			x[1] x[2] ....  x[Np]

			x[2] x[3] ....  x[Np+1]

			x[3] x[4] ....  x[Np+2]

	  				.

	  				.

	  				.

			x[N] x[N+1] ... x[Np+N-1] 


	'''
	output = np.zeros((N, Np), dtype='complex128')

	for j in prange(N) :
		output[j] = record[j : j + Np]

	return output



def tiled_window(window, N, Np):

	'''

		The tiled_window function takes a data tapering window 

				a[1] a[2] .... a[K]

		where K = N or Np and forms an N x Np matrix as follows. 

				If K = Np, then it forms

					a[1] a[2] .... a[Np]

					a[1] a[2] .... a[Np]

						.

						.

						.

					a[1] a[2] .... a[Np]

				Otherwise if K = N, then we form

					a[1]  a[1]  a[1] ....  a[1]

					a[2]  a[2]  a[2] ....  a[2]

					a[3]  a[3]  a[3] ....  a[3]

						.
			
						.
			
						.

					a[N]  a[N]  a[N] ....  a[N]

	'''

	
	if len(window) == N :

		return np.tile(window, Np).reshape(Np, N).transpose().astype('complex128')

	elif len(window) == Np :

		return np.tile(window, N).reshape(N, Np).astype('complex128')


@njit(parallel=True)
def form_conj(record, N, Np, conj):

	'''
		
		Form

			x*[1] x*[1] .... x*[1]

			x*[2] x*[2] .... x*[2]

			x*[3] x*[3] .... x*[3]

					.
					.
					.
			x*[N] x*[N] .... x*[N]

	'''

	inter = np.zeros(N, dtype='complex128')

	for j in prange(N):

		if not conj :

			inter[j] = np.conj(record[j + (Np // 2)])

		else :

			inter[j] = record[j + (Np // 2)]

	output = np.zeros((Np, N), dtype='complex128')

	for j in prange(Np) :

		output[j] = inter

	return output.transpose()

@njit(parallel=True)
def create_exponential(exp_args, N, Np):

	out = np.zeros((N, Np), dtype='complex128')

	for n in prange(N):

		out[n] = np.exp(-1j * 2 * np.pi * n * exp_args)

	return out

@njit(parallel=True)
def create_freqs(K, Q, N, Np, cycle_bool):

	freqs = np.zeros((N, Np), dtype='float64')

	if not cycle_bool:

		for i in prange(N):

			for j in prange(Np):

				freqs[i][j] = 0.5 * (Q[j] - K[i])

	else:

		for i in prange(N):

			for j in prange(Np):

				freqs[i][j] = Q[j] + K[i]


	return freqs




def SSCA(data_matrix, conj_mat, exp_mat, window_1, window_2):

	'''
		
		The main function of the script that implements the SSCA. 

		The variable data_matrix is the matrix formed from the record that was discussed before.

		Variable conj_matrix is the matrix of conjugates of the record. 

		Variables window_1 and window_2 are tiled data tapering windows.

	'''

	thread_count = mp.cpu_count()

	N, Np = data_matrix.shape

	inter_matrix = np.multiply(data_matrix, window_1)

	out_matrix = np.zeros((N, Np), dtype='complex128')

	# This step performs the Np-point FFT of each row of the product of the window and data.

	fft_obj_inter = ftw.FFTW(inter_matrix, out_matrix, axes=(1,), direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE',), threads=thread_count)

	fft_obj_inter()

	# We want centered ffts. So shift each row.

	out_matrix = np.fft.fftshift(out_matrix, axes=(1,))

	# Multiply everything.

	inter_matrix_2 = np.multiply(np.multiply(np.multiply(out_matrix, exp_mat), conj_mat), window_2)

	out_matrix_2 = np.zeros((N, Np), dtype='complex128')

	# Now perform the N-point FFT of each coloumn of the result.

	fft_obj_inter_2 = ftw.FFTW(inter_matrix_2, out_matrix_2, axes=(0, ), direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE',), threads=thread_count)

	fft_obj_inter_2()

	# Center the result in each row and return.

	out_matrix_2 = np.fft.fftshift(out_matrix_2, axes=(0,))

	return out_matrix_2

@njit(parallel=True)
def apply_ROI(corr_density, freqs, cycles, N, Np, sample_rate, data_center, freq_lo, freq_hi, conj):

	# sample_rate * freqs[i][j] + data_center, sample_rate * cycles[i][j]

	mask = np.zeros((N, Np),dtype='float64')
	center = 0.5 * (freq_lo + freq_hi)
	bandwidth = freq_hi - freq_lo

	for i in prange(N):

		for j in prange(Np):

			if not conj:
				abs_freq = sample_rate * freqs[i][j] + data_center
				abs_cycle_freq = sample_rate * cycles[i][j]
				if ((abs(abs_cycle_freq) <= bandwidth) and (abs(abs_freq-center) < 0.5 * (bandwidth - abs(abs_cycle_freq)))):
					mask[i][j] = 1.0

				elif ((abs(abs_cycle_freq) <= bandwidth) and (abs(abs_freq-center) == 0.5 * (bandwidth - abs(abs_cycle_freq)))):
					mask[i][j] = 0.5

				elif ((abs(abs_cycle_freq) == bandwidth) and (abs(abs_freq-center) == bandwidth)):
					mask[i][j] = 0.25

				else:
					mask[i][j] = 0.0

			else:
				abs_freq = sample_rate * freqs[i][j]
				abs_cycle_freq = sample_rate * cycles[i][j] + 2 * data_center
				if ((abs(abs_cycle_freq - 2 * center) <= bandwidth) and (abs(abs_freq) < 0.5 * (bandwidth - abs(abs_cycle_freq - 2 * center)))):
					mask[i][j] = 1.0

				elif ((abs(abs_cycle_freq - 2 * center) <= bandwidth) and (abs(abs_freq) == 0.5 * (bandwidth - abs(abs_cycle_freq - 2 * center)))):
					mask[i][j] = 0.5

				elif ((abs(abs_cycle_freq - 2 * center) == bandwidth) and (abs(abs_freq) == bandwidth)):
					mask[i][j] = 0.25

				else:
					mask[i][j] = 0.0

	return np.multiply(mask, corr_density)

def psd(x_vec, win, nfft):
	
#     windows_dict = {'bartlett' : np.bartlett, 'blackman' : np.blackman, 'hamming' : np.hamming, 'hanning' : np.hanning, 'kaiser' : lambda x : np.kaiser(x, beta)}

	ncpu = mp.cpu_count()

#     # pdb.set_trace()

	if len(x_vec) < nfft :
		
		in0 = np.pad(x_vec, (0, nfft - len(x_vec)), 'constant')

	else :
		
		in0 = x_vec[:nfft]
	
	fft_vec = np.zeros(nfft, dtype='complex128')

	fft_obj = ftw.FFTW(in0, fft_vec, direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE',), threads=ncpu, planning_timelimit=None)

	fft_obj()

	spec = np.fft.fftshift(np.multiply(fft_vec, np.conjugate(fft_vec)))

	window = win / np.sum(win)

	int_out = np.zeros(nfft, dtype='complex128')
#     # int_win = np.zeros(nfft, dtype='complex64')

	win_fft = np.fft.fft(window, nfft)

	spec = np.multiply(fft_obj(spec, int_out), win_fft)

	spec = np.fft.ifft(spec, nfft) / nfft

#     # spec = np.fft.fftshift(spec)

	return np.abs(spec)


@njit(parallel=True)
def create_psd_matrix(psd, conj, freqs, cycles, N, Np):
	output = np.zeros((N, Np), dtype='float64')
	nfft = len(psd)
	for i in prange(N):
		for j in prange(Np):
			freq_bin_1 = min(nfft // 2 + round(nfft * (freqs[i][j] + 0.5 * cycles[i][j])), nfft-1)
			freq_bin_2 = min(nfft // 2 + (round(nfft * (freqs[i][j] - 0.5 * cycles[i][j])) if not conj else round(nfft * (-freqs[i][j]+0.5*cycles[i][j]))), nfft-1)
			output[i][j] = (psd[freq_bin_1]*psd[freq_bin_2])**0.5
	return output 


# @njit(parallel=True)
# def coherence(signal, win1, win2, win3, N, Np, len_avg, psd_oversample, sample_rate, data_center, kaiser_beta1=10, kaiser_beta2=20, kaiser_beta3=10, conj=False, bandlimits=None, welch=False):

# 	windows_dict = {'bartlett' : np.bartlett, 'blackman' : np.blackman, 'hamming' : np.hamming, 'hanning' : np.hanning, 'kaiser1' : (lambda x : np.kaiser(x, kaiser_beta1)), 

# 	'kaiser2' : (lambda x : np.kaiser(x, kaiser_beta2)), 'kaiser' : (lambda x : np.kaiser(x, kaiser_beta3))}

# 	corr_density = scd_ssca(signal, win1, win2, N, Np, sample_rate, data_center, kaiser_beta1, kaiser_beta2, conj, bandlimits)

# 	power_spectrum = psd(signal, win3, len_avg, psd_oversample, kaiser_beta3) if not welch else np.fft.fftshift(sig.welch(signal, nperseg=1024, return_onesided=False)[1])

# 	coherence = np.zeros((N, Np), dtype='complex128')

# 	K = np.fft.fftshift(np.fft.fftfreq(N))

# 	Q = np.fft.fftshift(np.fft.fftfreq(Np))

# 	fun = lambda x : math.floor(x) if x > 0 else math.ceil(x)

# 	nfft = psd_oversample if not welch else len(power_spectrum)

# 	# pdb.set_trace()

# 	for indk in prange(N):

# 		for indq in prange(Q):

# 			alpha = K[indk] + Q[indq] 

# 			freq = 0.5 * (Q[indq] - K[indk]) 

# 			# freq_bin_1 =  psd_oversample // 2 + fun(psd_oversample * (freq + alpha * 0.5))

# 			# freq_bin_2 = min(psd_oversample // 2 + (fun(psd_oversample * (freq - alpha * 0.5)) if not conj else fun(psd_oversample * (-freq + alpha * 0.5))), psd_oversample-1)

# 			freq_bin_1 =  min(nfft // 2 + round(nfft * (freq + alpha * 0.5)), nfft-1)
# 			freq_bin_2 = min(nfft // 2 + (round(nfft * (freq - alpha * 0.5)) if not conj else round(nfft * (-freq + alpha * 0.5))), nfft-1)
# 			denom = (power_spectrum[freq_bin_1] * power_spectrum[freq_bin_2]) ** 0.5
# 			coherence[indk][indq] = (0.0 + 1j * 0.0) if (math.isclose(denom, 0.0)) else (corr_density[indk][indq] / denom)

# 	return coherence




# def extract_psd(output):

# 	N, Np = output.shape

# 	psd_dict = defaultdict(float)

# 	K = np.fft.fftshift(np.fft.fftfreq(N))

# 	Q = np.fft.fftshift(np.fft.fftfreq(Np))

# 	for indk, k in enumerate(K) :

# 		for indq, q in enumerate(Q):

# 			alpha = k + q

# 			freq = 0.5 * (q - k)

# 			if math.isclose(alpha, 0.0) :

# 				psd_dict[freq] = abs(output[indk][indq])

# 	return psd_dict
@njit(parallel=True)
def create_alpha_estimates(scd_result, cycles, N, Np):
	output = np.zeros((N, Np, 2), dtype='float64')
	for i in prange(N):
		for j in prange(Np):
			output[i][j][0] = scd_result[i][j]
			output[i][j][1] = cycles[i][j]
	return output


def sorted_alpha_estimates(output, cycles, N, Np):
	
	output = output.reshape((N * Np, 2))

	ind_arr = np.argsort(output[:, 1], axis=0, kind='stable')

	output = output[ind_arr]

	return output

def pick_relevant_values(output, freq_lo, freq_hi):
	
	bool_arr_1 = output[:, 1] > freq_lo
	red_arr_1 = np.compress(bool_arr_1, output, axis=0)
	bool_arr_2 = red_arr_1[:, 1] < freq_hi
	red_arr_2 = np.compress(bool_arr_2, red_arr_1, axis=0)
	
	return red_arr_2

def collect_alphas(output):
	
	return np.array(sorted(list(set(output[:, 1]))))

def pick_max_values(output, cycle_array):
	
	max_dict = {cycle : -1.0 for cycle in cycle_array}
	
	for out, freq in output:
	
		if (freq > cycle_array[0] and freq < cycle_array[-1] and max_dict[freq] < out):
			max_dict[freq] = out

	return max_dict

def sum_values(output, cycle_array):

	sum_dict = {cycle : 0.0 for cycle in cycle_array}

	for out, freq in output:

		if (freq > cycle_array[0] and freq < cycle_array[-1]):
			sum_dict[freq] += np.abs(out)

	return sum_dict


# def max_val_array(output, cycle_array):

# 	max_arr = np.zeros(len(cycle_array), dtype='float64')

# 	count = 0

# 	for out, freq in output:
		
# 		if (freq > cycle_array[0] and freq < cycle_array[-1] and max_arr[count] < out):
# 			max_arr[count] = out
# 		count += 1

# 	return max_arr

class CycloGram:
	
	def __init__(self, beta_1, beta_2, N, Np, conj):
		
		self.N = N
		self.Np = Np
		self.beta_1 = beta_1
		self.beta_2 = beta_2
		self.conj = conj
		K = np.fft.fftshift(np.fft.fftfreq(N))
		Q = np.fft.fftshift(np.fft.fftfreq(Np))
		freqs = create_freqs(K, Q, N, Np, False)
		cycles = create_freqs(K, Q, N, Np, True)
		self.freqs = freqs
		self.cycles = cycles
		
	def cyclogram(self, signal_1, signal_2=None, avg=False):
		if signal_2 is None:
			signal_2 = signal_1
		window_1 = np.kaiser(self.Np, self.beta_1)
		window_2 = np.kaiser(self.N, self.beta_2)
		norm_1 = (np.sum(np.square(window_1)))**0.5
		norm_2 = np.sum(window_2)
		window_1 = window_1 / norm_1
		window_2 = window_2 / norm_2
		win_1 = tiled_window(window_1, self.N, self.Np)
		win_2 = tiled_window(window_2, self.N, self.Np)
		exp_args = np.fft.fftshift(np.fft.fftfreq(self.Np))
		exp_mat = create_exponential(exp_args, self.N, self.Np)
		ssca = np.zeros((self.N, self.Np), dtype='float64')
		if avg :
			ncol = (2 * len(signal_1)) // (self.N + self.Np)
			for ind in prange(ncol - 2):
				new_sig = signal_1[ind * (self.N + self.Np) // 2 :]
				new_sig_2 = signal_2[ind * (self.N + self.Np) // 2 :]
				data_matrix = create_matrix(new_sig, self.N, self.Np)
				conj_mat = form_conj(new_sig_2, self.N, self.Np, self.conj)
				ssca += np.abs(SSCA(data_matrix, conj_mat, exp_mat, win_1, win_2))
			ssca = ssca / ncol
		else:
			data_matrix = create_matrix(signal_1, self.N, self.Np)
			conj_mat = form_conj(signal_2, self.N, self.Np, self.conj)
			ssca = np.abs(SSCA(data_matrix, conj_mat, exp_mat, win_1, win_2))
		return ssca
		
		
