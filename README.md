# 3D-tagless-object-reconstruction

The hardware system consists of several commercial off-the-shelf (COTS) radio-frequency identification (RFID) readers, and a multitude of passive RFID tags. The transmission between the reader and tags are described as follows:

The reader actively transmits modulated continuous-wave signals through passive antennas in the frequency range of 902-928 MHz. Baseband signal is modulated onto the carrier wave according to the Generation-2 (Gen2) Electronic Product Code (EPC) protocol. 50 channels, each with 500 kHz bandwidth, are successively used for transmission, where frequency hopping spread spectrum (FHSS) is adopted and the channel sequence is pre-determined by the random number generator in the reader.
The passive tags do not contain any battery, as they are able to harvest energy from the incoming reader's transmitted signal to power up the tag chip. The tag chip enables the tag to generate its own baseband signal (including tag ID and other command-specific information compliant with the Gen2 EPC protocol), and modulate its baseband signal onto the incoming reader's transmitted signal by changing its antenna's reflection coefficient. This process is called backscattering.
The reader receives the tag backscattering signal and decodes tag ID, as well as other RF information from the tags. Most current RFID readers support received signal strength indicator (RSSI) and phase information, which are the two most important RF parameters for this project.


This project aims to reconstruct the reflectivity image of the capture volume of interest, based on the received RF parameters from backscattering signal from RFID tags. Conventional matched filtering (MF) is used as the inverse solution for this project, with a novel calibration method to eliminate background clutter and proposed postprocessing techniques to cancel out noisy channels due to background dynamics and channels affected by line-of-sight (LoS) blockage by large target objects.

The experiment is carried out as follows:
1. We first run the hardware system to collect the data when no object is in the capture volume.
2. We then run the hardware system to collect the data when the object(s) is(are) in the capture volume.
3. Then, data processing is carried out, where we reconstruct the reflectivity image in the capture volume using MF.

The first two steps were carried out beforehand, with collected data stored in the folder "DataCollected" (committed to this repository), and the MATLAB scripts and functions committed into this repository achieve the third step above.
If you download the repository as a package, you would like to run "ReconstructScript.m", and you will get the reconstructed images of the object at specified locations. You might need to change lines 21 and 22 of "ReconstructScript.m" according to the directory of your computer after downloading the repository.
