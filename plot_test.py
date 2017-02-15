import numpy as np
import matplotlib.pyplot as plt

analy = np.load('./results/Analytical_correct_all.npz')
calc_ = np.load('./results/CalcPotential4_correct_all.npz')

for key_val in analy.keys():
    plt.close('all')
    plt.subplot(121)
    plt.imshow(analy[key_val].reshape(180, 180), cmap=plt.cm.PRGn)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(calc_[key_val], cmap=plt.cm.PRGn)
    plt.colorbar()
    plt.savefig(key_val + '.png')
