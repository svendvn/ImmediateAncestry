import tmhmm
import numpy
import os
import operator
import functools
import matrix_generation
import operator
import matplotlib
import matplotlib.pyplot as plt
from math import log
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection

# -----------------------------------------------------

color_lookup = {
                  'gneg': (1., 1., 1.),
                'gpos25': (.6, .6, .6),
                'gpos50': (.4, .4, .4),
                'gpos75': (.2, .2, .2),
               'gpos100': (0., 0., 0.),
                  'acen': (.8, .4, .4),
                  'gvar': (.8, .8, .8),
                 'stalk': (.9, .9, .9),
               }

height = 0.9
spacing = 0.9

def ideograms(fn):
    last_chrom = None
    fin = open(fn)
    fin.readline()
    xranges, colors = [], []
    ymin = 0

    for line in fin:
        chrom, start, stop, label, stain = line.strip().split('\t')
        start = int(start)
        stop = int(stop)
        width = stop - start
        if chrom == last_chrom or (last_chrom is None):
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            last_chrom = chrom
            continue

        ymin += height + spacing
        yrange = (ymin, height)
        yield xranges, yrange, colors, last_chrom
        xranges, colors = [], []
        xranges.append((start, width))
        colors.append(color_lookup[stain])
        last_chrom = chrom

    # last one
    ymin += height + spacing
    yrange = (ymin, height)
    yield xranges, yrange, colors, last_chrom

fig = plt.figure()
ax = fig.add_subplot(111)
d = {}
yticks = []
yticklabels = []

# ideogram.txt downloaded from UCSC's table browser
for xranges, yrange, colors, label in ideograms('ideogram.txt'):
    coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
    ax.add_collection(coll)
    center = yrange[0] + yrange[1]/2.
    yticks.append(center)
    yticklabels.append(label)
    d[label] = xranges

ax.axis('tight')
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)
ax.set_xticks([])
plt.show()


