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


def plot_EPSC(EPSC_x,files_5min,files_8min):
    
    EPSC_5y=np.zeros((len(files_5min),45))
    EPSC_8y=np.zeros((len(files_8min),45))

    a =100 # lower bound
    b =150 # upper bound

    f=0
    while f<len(files_8min):
        EPSC_5y[f,:] = Read_Column_Two_File(files_5min[f])
        EPSC_8y[f,:] = Read_Column_Two_File(files_8min[f])
        f+=1

    #re-scale data and then calculate mean and std
    A = np.min(EPSC_5y)
    B = np.max(EPSC_5y)
    C = np.min(EPSC_8y)
    D = np.max(EPSC_8y)

    f=0
    while f<len(files_8min):
        EPSC_5y[f,:] = [a+(x-A)*(b-a)/(B-A) for x in EPSC_5y[f,:]]           
        EPSC_8y[f,:] = [a+(x-C)*(b-a)/(D-C) for x in EPSC_8y[f,:]]
        f+=1


    EPSC_5y_mean = np.mean(EPSC_5y,axis=0)
    EPSC_5y_std = EPSC_5y.std(0)
    EPSC_8y_mean = np.mean(EPSC_8y,axis=0)
    EPSC_8y_std = EPSC_8y.std(0)

    plt.figure(figsize=(10,7))
    plt.errorbar(EPSC_x, EPSC_5y_mean, yerr=EPSC_5y_std ,linewidth=2.0, marker='o', markerfacecolor='none', markeredgewidth=2.0, markersize=12 ,color='#2197d5', label='Disinhibition period of 5 minutes')    
    plt.errorbar(EPSC_x, EPSC_8y_mean, yerr=EPSC_8y_std ,linewidth=2.0, marker ='s', markerfacecolor='none', markeredgewidth=2.0,markersize=12, color='darkorange', label='Disinhibition period of 8 minutes')
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
    plt.legend(fontsize=24,loc=3, frameon=False)

    plt.show()



# ----------  PLOT FIGURES  ---------- #

#FIGURE 4(D)
EPSC_x = Read_Column_One_File('EPSC_8mint_v1.txt')
files_5min = ['EPSC_5mint_DUR3_v1.txt', 'EPSC_5mint_DUR3_v2.txt',  'EPSC_5mint_DUR3_v6.txt', 'EPSC_5mint_DUR3_v8.txt', 'EPSC_5mint_DUR3_v9.txt' ]
files_8min = ['EPSC_8mint_DUR3_v1.txt', 'EPSC_8mint_DUR3_v2.txt', 'EPSC_8mint_DUR3_v6.txt', 'EPSC_8mint_DUR3_v8.txt', 'EPSC_8mint_DUR3_v9.txt' ]


plot_EPSC(EPSC_x,files_5min,files_8min)
