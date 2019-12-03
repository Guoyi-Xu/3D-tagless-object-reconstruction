# 3D-tagless-object-reconstruction

This project aims to reconstruct the object from the signals received by RFID reader antennas, which are backscattered from RFID tags.
Both the tags and receiver antennas are placed around the capture volume, which is around 3.6 m by 3.6 m by 3 m.

The experiment is carried out as follows:
1. We first run the hardware system to collect the data when no object is in the capture volume.
2. We then run the hardware system to collect the data when the object(s) is(are) in the capture volume.
3. Then, post-processing of the collected data is carried out, in which we reconstruct the object location and size using our algorithm.

The first two steps were carried out beforehand, with collected data stored in the folder "DataCollected" (committed to this repository), and the MATLAB scripts and functions committed into this repository achieve the third step above.
If you download the repository as a package, you would like to run "ReconstructScript.m", and you will get the reconstructed images of the object at specified locations. You might need to change lines 21 and 22 of "ReconstructScript.m" according to the directory of your computer after downloading the repository.
