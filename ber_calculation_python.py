# -*- coding: utf-8 -*-
"""BER_calculation_surya.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1izVLUJz_-0iXUvcFXs3K3RZPAHOA0S-E
"""

from google.colab import drive

drive.mount('/content/gdrive/',force_remount=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import h5py

file_path = '/content/gdrive/MyDrive/mimo_data_surya_final.mat'

with h5py.File(file_path, 'r') as file:
  # data=file['tx_data']
  print(file.keys())

  tx_data_init=file['tx_data'][:]
  rx_data_init=file['rx_data'][:]
  snr_values_init=file['snr_values'][:]

tx_data=tx_data_init.T
rx_data=rx_data_init.T
snr_values=snr_values_init.T

# 7000 iterations are there so calculate the BER for each iteration (with 400 bits each)

ber_values=[]
snr_levels=list(np.unique(snr_values))
errors=0
total_bits=0

for i in range(len(tx_data)): #looping through each iteration

  for j in range(tx_data.shape[1]): # looping through each bit in each iteration
    if tx_data[i][j]!=rx_data[i][j]:
      errors+=1
    total_bits+=1
  if (i+1001)%1000 ==0: # this is to reset the errors and bits for each SNR

    print(total_bits)
    ber_values.append(errors/total_bits)

    errors=0
    total_bits=0
print("BER Values for each SNR: " + str(ber_values))

plt.figure(figsize=(10, 6))
plt.plot(snr_levels, ber_values, 'bo-', linewidth=2, markersize=8)


plt.xlabel('SNR (dB)', fontsize=14)
plt.ylabel('Bit Error Rate (BER)', fontsize=14)
plt.title('BER vs. SNR for 4x4 MIMO System with BPSK and Rayleigh Fading', fontsize=16)


plt.grid(True, which='both', linestyle='--', linewidth=0.7)


plt.yscale('log')


plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig('\content\gdrive\My Drive\BER_vs_SNR.jpeg', format='jpeg', dpi=300)
plt.savefig('\content\gdrive\My Drive\BER_vs_SNR.eps', format='eps', dpi=300)

plt.show()

