These files are a list of vectors that most efficiently bin a sphere. It's based on the HEALpix discretisation, but since I only use basically one function from this library (pix2vec) I have elected to reduce the setup requirements by just giving the 2nd, 4th, 8th and 16th order discretisations. 

The formula for the number of pixels is, N_pix = 12 n^2 where n is the order of the discretisation which needs to be a power of two. I only go up to 8 here because that's easily enough to get a nice smooth distribution. 

This tends to be the fastest and cleanest way to do the binning of the full directional distribution in 3D.
