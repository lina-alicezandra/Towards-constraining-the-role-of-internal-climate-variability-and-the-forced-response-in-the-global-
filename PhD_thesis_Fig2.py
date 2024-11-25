# -*- coding: utf-8 -*-
"""

This is to plot Figure 2 of PhD thesis
"""

import xarray as xr
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


def lowess(x, y, f=3. / 4., iter=3):
    
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    check the websit for reference https://www.weisang.com/en/documentation/loessandlowessalgorithm_en/
    
    # Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
    #
    # License: BSD (3-clause)
    """
    
    n = len(x)
    r = int(math.ceil(f * n))
    # h represent the largest distance from xi to other x in the smoothing range
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)] 
    # for one iteration, the whole x - xi, the distance would be sorted, then only select the furthest distance with within smoothing range f
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    # in this way, for each xi, its distance to other x within the smoothing range is selected, larger than this distance is set to 1 and then to 0.
    w = (1 - w ** 3) ** 3
    # w is the weight matrix, for each column it saves for each xi the weight of each x within the smoothing range.

    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = np.linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]
        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2
    return yest

def filt(carbon):
    train_x = np.arange(400, len(carbon)+400)
    train_y = np.asarray(carbon).ravel()
    f1 = 0.25# 0.25 for 54 years 0.45 for 30 years
    yest = lowess(train_x, train_y, f = f1, iter=3)
    carbon2 = train_y - yest
    
    return yest, carbon2

def carbon_treatment(f_out, carbon):
    date = pd.date_range(start='1959-01', end='2024-1', freq='m')
    df = pd.read_excel(carbon)
    df = df.set_index(date)
    print(df)
    pre, carbon2 = filt(df['de-seasonalized'])
    carbon2 = pd.DataFrame(carbon2, index=date)
    
    annu = carbon2.resample(rule='Y').sum()
    
    date2 = np.arange(1959, 2024)
    
    
    fig, axs =  plt.subplots(5, 1, figsize=(7, 10))
    
    axs[0].plot(date, df['monthly_average'], c='#2c7fb8')
    axs[0].set_ylabel(r'$ppm \cdot month^{-1}$')
    
    axs[0].set_ylim(300, 450)
    axs[0].set_title('(a)')
    axs[0].text(1980, 420.5,  '$CO_2$ concentration, monthly average', color = '#3D3D3D', fontsize=10)
    
    axs[1].plot(date,df['seasonal'], c='#2c7fb8')
    axs[1].set_ylabel(r'$ppm \cdot month^{-1}$')
    axs[1].set_ylim(-5, 7)
    axs[1].set_title('(b)')
    axs[1].text(1980, 4.5, '$CO_2$ concentration, seasonal trend', color = '#3D3D3D', fontsize=10)
    
    axs[2].plot(date, pre, c='#2c7fb8')
    axs[2].set_ylabel(r'$ppm \cdot month^{-1}$')
    axs[2].set_ylim(300, 450)
    axs[2].set_title('(c)')
    axs[2].text(1980, 420,  '$CO_2$ concentration, decadal trend', color = '#3D3D3D', fontsize=10)
    
    axs[3].plot(date, carbon2, c='#2c7fb8')
    axs[3].set_ylabel(r'$ppm \cdot month^{-1}$')
    axs[3].set_ylim(-3, 4)
    axs[3].set_title('(d)')
    axs[3].text(1980, 2.5,  '$CO_2$ concentration, monthly anomalies', color = '#3D3D3D', fontsize=10)
    
    
    axs[4].plot(date2, annu, c='#2c7fb8')
    axs[4].set_ylabel(r'$ppm \cdot year^{-1}$')
    axs[4].set_ylim(-15, 20)
    axs[4].set_title('(e)')
    axs[4].set_xlabel('Year')
    axs[4].text(1975, 12.5,  '$CO_2$ concentration, inter-annual variability', color = '#3D3D3D', fontsize=10)
    
    for a in axs.flat:
        a.label_outer()
    fig.tight_layout(pad=2.0)

    plt.savefig(f_out + "CO2_decomposition.pdf")
    
    plt.show()
    

    
    # car_new = pd.DataFrame({'Date':df.index, 'AGR':agr_dr})
    # car_new.to_csv(f_out + 'Global_Carbon_1959_2017_detr.csv', index=False)
    
if __name__ == '__main__':
    carbon_orig = ".../Book1.xlsx"
    carbon_treatment('/', carbon_orig)
    
    
    
    
    
    

  
