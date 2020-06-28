import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider#, Button, RadioButtons

if __name__ == "__main__":

    data = np.load('/home/evan/Documents/Work/JackmanLab/EPSCs.200615_3.npz')

    fig, ax = plt.subplots()
    E = ax.imshow(data['arr_4'][:,0,:], interpolation = None, cmap = cm.coolwarm)
    fig.colorbar(E, shrink = .3, label = 'Error')

    axcolor = 'white'
    axT_D = plt.axes([0.1, 0.1, 0.5, 0.05], facecolor=axcolor)
    T_D_0 = data['arr_2'][0]
    delta_T_D = np.diff(data['arr_2'])[0]
    sT_D = Slider(axT_D, 'T_D', min(data['arr_2']), max(data['arr_2']), valinit=T_D_0, valstep=delta_T_D)

    ax.set_xlabel("k_max")
    ax.set_ylabel("K_D")
    # ax.set_xlim(min(data['arr_3']), max(data['arr_3']))
    # ax.set_ylim(min(data['arr_1']), max(data['arr_1']))
    # ax.set_xticks(data['arr_3'])
    # ax.set_yticks(data['arr_1'])

    def update(val):
        T_D = sT_D.val
        E.set_data(data['arr_4'][:,int(T_D)-1,:])
        fig.canvas.draw()

    sT_D.on_changed(update)

    plt.show()

    data.close()
