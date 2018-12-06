import numpy as np
import os
import matplotlib.pyplot as plt

def plot_slices(outputdir, var, log):

    runname = 'small'

    data = np.load(os.path.join(os.path.dirname(__file__), outputdir, runname + '.npz'))
    nx, ny, nz = data[var].shape

    data_var = np.log(data[var]) if log else data[var]


    fig, axes = plt.subplots(2,2)

    plt.subplot(221)
    plt.pcolor(data_var[:,0,:].T)
    plt.title('y = 0')
    plt.colorbar()

    plt.subplot(222)
    plt.pcolor(data_var[:,ny//4,:].T)
    plt.title('y = ' + str(ny//4))
    plt.colorbar()

    plt.subplot(223)
    plt.pcolor(data_var[:,ny//2,:].T)
    plt.title('y = ' + str(ny//2))
    plt.colorbar()
    
    plt.subplot(224)
    plt.pcolor(data_var[:,-1,:].T)
    plt.title('y = ' + str(ny))
    plt.colorbar()

    fig.suptitle(outputdir + ': ' + var)
    plt.tight_layout()
    plt.draw()

plot_slices('small', 'facies', False)
plot_slices('small', 'azim', False)
plot_slices('small', 'dip', False)
plot_slices('small', 'ha', False)
plt.show()
