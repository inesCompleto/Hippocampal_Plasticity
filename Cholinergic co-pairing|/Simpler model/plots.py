import numpy as np
import matplotlib.pyplot as plt

# ----------  FUNCTIONS  ---------- #

def Read_Column_One_File(file_name):
    with open(file_name, 'r') as data:
        x = []
        for line in data:
            p = line.split()
            x.append(float(p[0]))
            
        x=np.array(x)
            
    return x

def Read_Column_Two_File(file_name):
    with open(file_name, 'r') as data:
        y = []
        for line in data:
            p = line.split()
            y.append(float(p[1]))
        y=np.array(y)

    return y


def plot_EPSC(EPSC_x,files):
    
    EPSC_y=np.zeros((len(files),40))

    a =100 # lower bound
    b =150 # upper bound

    f=0
    while f<len(files):
        EPSC_y[f,:] = Read_Column_Two_File(files[f])
        f+=1

    #re-scale data and then calculate mean and std
    C = np.min(EPSC_y)
    D = np.max(EPSC_y)

    f=0
    while f<len(files):         
        EPSC_y[f,:] = [a+(x-C)*(b-a)/(D-C) for x in EPSC_y[f,:]]
        f+=1

    EPSC_y_mean = np.mean(EPSC_y,axis=0)
    EPSC_y_std = EPSC_y.std(0)

    plt.figure(figsize=(10,7))
    plt.errorbar(EPSC_x, EPSC_y_mean, yerr=EPSC_y_std ,linewidth=2.0, marker ='s', markerfacecolor='none', markeredgewidth=2.0,markersize=12, color='black')
    plt.ylim(50, 170)
    plt.yticks(np.arange(50, 160+1, 25))
    plt.xlim(-1, 45.5)
    plt.xticks(np.arange(0, 41+1, 5)) 
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Normalized EPSC amplitude', fontsize=28)
    plt.xlabel('time (min)', fontsize=28, x=1, y=1)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()




# ----------  PLOT FIGURES  ---------- #

EPSC_x = Read_Column_One_File('EPSC_simpler_v1.txt')
files = ['EPSC_simpler_v1.txt', 'EPSC_simpler_v2.txt', 'EPSC_simpler_v3.txt' ]
plot_EPSC(EPSC_x, files)


