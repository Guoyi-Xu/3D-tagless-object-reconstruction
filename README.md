# 3D-tagless-object-reconstruction

This project aims to reconstruct the object from the signals received by RFID reader antennas, which are backscattered from RFID tags.
Both the tags and receiver antennas are placed around the capture volume, which is around 3.6 m by 3.6 m by 3 m.

The experiment is carried out as follows:
We first run the hardware system to collect the data when no object is in the capture volume.
We then run the hardware system to collect the data when the object(s) is(are) in the capture volume.
Then, post-processing of the collected data is carried out, in which we reconstruct the object location and size using our algorithm.
