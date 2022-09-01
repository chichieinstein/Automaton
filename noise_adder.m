function iq = noise_adder(input, snr)
iq = awgn(input, snr, 'measured');
end