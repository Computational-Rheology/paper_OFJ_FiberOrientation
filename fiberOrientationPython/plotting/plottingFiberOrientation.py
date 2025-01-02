
#%% Modules
import matplotlib.pyplot as plt

#%% Funcs
def plotting(time, A2):
    fZ = 14
    lW = 2.0
    plt.rc('text', usetex=True) # If you have problems with latex, switch to usetex=False
    plt.rc('font', family='serif')
    plt.rcParams["figure.figsize"] = (10, 7)
    plt.rcParams['axes.linewidth'] = lW
    
    fig, ax = plt.subplots()

    ax.plot(time, A2[:, 0], label="$A_{11}$", color="darkblue")
    ax.plot(time, A2[:, 4], label="$A_{22}$", color="darkred")
    ax.plot(time, A2[:, 8], label="$A_{33}$", color="darkgreen")
    
    plt.ylim([0, 1])
    plt.xlim([min(0, min(time)), max(time)])
    
    plt.yticks(fontsize = fZ)
    plt.xticks(fontsize = fZ)
    
    ax.tick_params(which = 'both', direction='out', length=5, width=lW, labelsize=fZ+2)
    plt.xlabel("Time [s]", fontsize=fZ+4)
    plt.ylabel("$A_{ij}$ components", fontsize=fZ+4)
    plt.grid(visible=True, which='major', color='black', linestyle='-')
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='darkgrey', linestyle='--', alpha=0.5)
    plt.tight_layout()

    plt.show()

    return ax
