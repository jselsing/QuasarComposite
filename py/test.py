


import matplotlib.pyplot as pl
import numpy as np
import seaborn as sns; sns.set_style('ticks')
print(sns.color_palette())
cmap = sns.cubehelix_palette(n_colors=3, start=2, rot=1.2, gamma=1.0, hue=0.8, light=0.80, dark=0.20, reverse=True, as_cmap=False)

x = np.arange(0, 100, 0.1)

pl.plot(x, 2*x, color = cmap[0])
pl.plot(x, 3*x, color = cmap[1])
pl.plot(x, 4*x, color = cmap[2])
pl.show()