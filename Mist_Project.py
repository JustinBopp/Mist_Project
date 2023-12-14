# %%
import numpy as np
import matplotlib.pyplot as plt
import read_mist_models as rmm
import pandas as pd
import scipy as spy
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic
import argparse

parser = argparse.ArgumentParser(description="Mist Project")
parser.add_argument('file', type=str, metavar='', help='Name of the File  \n  Example: "00100" which represents 1.00 Solar Mass')
parser.add_argument('process',type=str, metavar='', help="\n The process you would like to run \n Options: Interpolation, Integration, Visualization")
parser.add_argument('style',type=str, metavar='', help="\n The type based on Process \n Interpolation: 'linear', 'gpr',‘nearest’, ‘nearest-up’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, or ‘next’. ‘zero’, ‘slinear’, ‘quadratic’ \n Integration: 'trapezoid' , 'simpson', 'midpoint \n Visualization: 'compare','hr','core',radius,'surface' ")
parser.add_argument('data_1', type=str,nargs='?',default='star_age', metavar='', help='\n X value \n  Default: "star_age" ')
parser.add_argument('data_2', type=str,nargs='?',default='log_L', metavar='', help='\n Y value  \n  Default: "log_L"')

args = parser.parse_args()


def MIST_project(file,process, style,data_1,data_2):
    eep = rmm.EEP(file+'M.track.eep')
    initial_mass=eep.minit
    print ('initial mass: ', eep.minit)
    data1=eep.eeps[data_1]
    data2=eep.eeps[data_2]

    if process == 'Interpolation':
        rand = np.random.choice(np.arange(len(data1)),50, replace=False)
        rand= np.sort(rand)
        xnew  = np.linspace(np.min(data1[rand]),np.max(data1[rand]),100 )

        plt.plot(data1[rand], data2[rand], '--' , label = "original")
        
        if style == "linear":
            ynew = linear_interp(data1[rand],data2[rand],xnew)
            plt.plot(xnew,ynew,'.', label = 'linear')
        elif style == "gpr":
            X_scaled = data1[rand].reshape(-1,1)
            kernel = RationalQuadratic()
            gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, n_restarts_optimizer=9)
            gpr.fit(X_scaled,data2[rand])
            Xnew_scaled = xnew.reshape(-1,1)
            f_gpr = gpr.predict(Xnew_scaled, return_std=True)
            ynew_gpr = np.asarray(f_gpr[0])
            plt.plot(xnew, ynew_gpr,'.', label="gpr")
        else:
            f= spy.interpolate.interp1d(data1[rand], data2[rand], kind=style)
            plt.plot(xnew, f(xnew), '.',markersize=4, label=style)
        plt.title(process)
        plt.xlabel(data_1)
        plt.ylabel(data_2)
        plt.legend()
        plt.show(block=True)

    if process == 'Integration':
        star_age = np.log(eep.eeps['star_age'])
        log_L = eep.eeps['log_L']
        phase = eep.eeps['phase']
        ms = np.argwhere(phase == 0)
        post_ms = np.argwhere(phase >= 2)
        if style == 'trapezoid':
            print('Total Energy on Main Sequence:',trap(ms,file))
            print('Total Energy on Post-Main Sequence:',trap(post_ms,file))
        elif style == 'simpson':
            print('Total Energy on Main Sequence:',simp(ms,file))
            print('Total Energy on Post-Main Sequence:',simp(post_ms,file))
        elif style == 'midpoint':
            print('Total Energy on Main Sequence:',midpoint(ms,file))
            print('Total Energy on Post-Main Sequence:',midpoint(post_ms,file))
        plt.plot(star_age[ms], log_L[ms],'-', label="Main Sequence")
        plt.plot(star_age[post_ms], log_L[post_ms],'-', label="Post Main Sequence")  
        plt.title('initial mass: '+ str(initial_mass))
        plt.xlabel(data_1)
        plt.ylabel(data_2)
        plt.legend()
        plt.show()
    if process == 'Visualize':
        if style == 'compare':
            plt.plot(data1,data2,'-')
            plt.title('initial mass: '+ str(initial_mass))
            plt.xlabel(data_1)
            plt.ylabel(data_2)
            plt.show()
        if style == 'hr': 
            eep.plot_HR(color='Black', phases=[0, 6], phasecolor=['Red', 'Blue'])
            plt.show()
        if style == 'core':         
            star_age = eep.eeps['star_age']
            center_h1 = eep.eeps['center_h1']
            center_he4 = eep.eeps['center_he4']
            center_c12 = eep.eeps['center_c12']
            plt.plot(star_age, center_h1, label='H1')
            plt.plot(star_age, center_he4, label='He4')
            plt.plot(star_age, center_c12, label='C12')
            plt.xlabel('Star Age')
            plt.ylabel('Mass Fraction')
            plt.axis([1e7, 1.5e10, 1e-6, 3])
            plt.xscale('log')
            plt.yscale('log')
            leg = plt.legend(loc=3, fontsize=16)
            plt.show()
        if style == 'radius':
            star_age = eep.eeps['star_age']
            radius = eep.eeps['log_R'] 
            plt.plot(star_age, radius )
            plt.xlabel('Star Age (log)')
            plt.ylabel('Radius (log)')
            plt.xscale('log')
            plt.show()    
        if style == 'surface':
            star_age = eep.eeps['star_age']
            surface_h1 = eep.eeps['surface_h1']
            surface_he3 = eep.eeps['surface_he3']
            surface_he4 = eep.eeps['surface_he4']
            surface_c12 = eep.eeps['surface_c12']
            surface_o16 = eep.eeps['surface_o16']
            plt.plot(star_age, surface_h1, label='H1')
            plt.plot(star_age, surface_he3, label='He3')
            plt.plot(star_age, surface_he4, label='He4')
            plt.plot(star_age, surface_c12, label='C12')
            plt.plot(star_age, surface_o16, label='O16')
            plt.xlabel('Star Age')
            plt.ylabel('Mass Fraction')
            plt.xscale('log')
            plt.yscale('log')
            leg = plt.legend(loc=3, fontsize=16)
            plt.show() 

def linear_interp(x,y,xnew):
        ynew=np.zeros(xnew.size)     
        n = x.size  
        for i in np.arange(xnew.size):
            
            ind = max(np.argwhere(x <= xnew[i]))   

            if (ind == n-1):  
                ind=ind-1

            m = (y[ind+1]-y[ind])/(x[ind+1]-x[ind])   
            ynew[i] = m*(xnew[i]-x[ind]) + y[ind]       
        return ynew


def trap(type,file):
    eep = rmm.EEP(file+'M.track.eep')
    star_age = np.log(eep.eeps['star_age'])
    log_L = eep.eeps['log_L']
    phase = eep.eeps['phase']
    ms = np.argwhere(phase == 0)
    post_ms = np.argwhere(phase >= 2)
    x = (star_age[type])
    y = 10**log_L[type] *3.846e26
    integral = 0.5 * sum((y[i] + y[i + 1]) * (x[i + 1] - x[i]) for i in range(len(x) - 1))
    return integral


def simp(type,file):
    eep = rmm.EEP(file+'M.track.eep')
    star_age = np.log(eep.eeps['star_age'])
    log_L = eep.eeps['log_L']
    phase = eep.eeps['phase']
    ms = np.argwhere(phase == 0)
    post_ms = np.argwhere(phase >= 2)
    integral = 0
    x = (star_age[type])
    y = 10**log_L[type] *3.846e26
    for i in range(0, len(x) - 2, 2):
        h = x[i + 2] - x[i]
        integral += (h / 3) * (y[i] + 4 * y[i + 1] + y[i + 2])
    return integral


def midpoint(type,file):
    eep = rmm.EEP(file+'M.track.eep')
    star_age = np.log(eep.eeps['star_age'])
    log_L = eep.eeps['log_L']
    phase = eep.eeps['phase']
    ms = np.argwhere(phase == 0)
    post_ms = np.argwhere(phase >= 2)
    x = (star_age[type])
    y = 10**log_L[type] *3.846e26
    dx = np.diff(x)  
    x_mid = (x[:-1] + x[1:]) / 2
    integral = np.sum(y[:-1] * (x[1:] - x[:-1])) 
    return integral

if __name__ == '__main__':
    MIST_project(args.file, args.process, args.style, args.data_1, args.data_2)


