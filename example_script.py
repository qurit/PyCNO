# %% Example Usage
import matplotlib.pyplot as plt
swept_parameters = ['Tumor1Volume', 'Rden_Tumor1']
swept_values = [(Tumor1Volume,Rden_Tumor1)
                for Tumor1Volume in np.linspace(0.1, 0.5, 24)
                for Rden_Tumor1 in np.arange(0.1, 0.5, 0.1)
                ]

time, TACs = runPBPK(swept_parameters=swept_parameters, swept_values=swept_values)
# plt.plot(time, TACs)
# plt.show()