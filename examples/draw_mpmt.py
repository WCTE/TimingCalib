"""
Draw a single mPMT in 3D


"""

import k3d
import numpy as np

from TimingCalib.MPMT import MPMT
from TimingCalib.PMT import PMT
from TimingCalib.LED import LED

my_mpmt = MPMT(None, 'M2', {}, {})

origins = []
vectors = []

for pmt in my_mpmt.pmts:
    location, direction = pmt.get_z_orientation('true')

    lv = [direction[i]*100 for i in range(len(direction))]
    origins.append(location)
    vectors.append(lv)

plt_vectors = k3d.vectors(origins=origins, vectors=vectors)

plot = k3d.plot()
plot += plt_vectors
plot.display()