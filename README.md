## Introduction

This code evaluates two Digital Cancellation methods for a Full Duplex system, in the Frequency and Time domain respectivelly. To that end, two
WiFi frames will be generated, passed through a TGn channel Model A, and AWGN. Impairments will be incorporated at the receiver such as Frequency Offset or IQ Imbalanced, in order to emulate a real radio behavior.

Since no WLAN System Toolbox was available at the moment of implementing this code, the work in [1] was taken as a baseline and their C++ code was addapted to Matlab accordingly.

## Simulation configuration

The executable of the project is *FD_full_duplex_simulation_refactor*, which loads a prestored 802.11a frame with default configurations (Modulation, encoding rate, etcetera), generates a desired and interference signal and adds the pertinent RF impairments to it. The Interference emulates the signal the device transmitts and slips through its own receiver (Full-Duplex mode).

The parameters can be configured in the main script prior execution:

- **SNR**: Signal-to-Noise Ratio for the desired signal in dB
- **Prate**: Power Ratio between the desired signal and the interference in linear scale.
- **offset**: Time Offset between the Desired and the Interfering signal.
- **can_meth**: Cancellation method: 'freq' or 'time'
- **FO**: Impairment at RX: Frequency offset
- **IQImb_gain**: Impairment at RX: IQ Imbal - Gain missmatch
- **IQImb_phase**: Impairment at RX: IQ Imbal - Phase Missmatch

## References

[1] Bastian Bloessl, Michele Segata, Christoph Sommer, and Falko Dressler. 2013. An IEEE 802.11a/g/p OFDM receiver for GNU radio. In Proceedings of the second workshop on Software radio implementation forum (SRIF '13). ACM, New York, NY, USA, 9-16. Code: https://github.com/bastibl/gr-ieee802-11.git