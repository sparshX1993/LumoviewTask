# Author : Sparsh Tiwari.

import csv
from scipy.interpolate import interp1d
import math
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt


def cart2pol(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)

# Data reading.
with open('out_startplatz_cut.txt') as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    NewRotations = []
    MeasurementQualities = []
    MeasurementAngles = []
    MeasurementDistances = []
    MeasurementDistancesX = []
    MeasurementDistancesY = []
    for row in readCSV:
        NewRotation = bool(row[0])
        MeasurementQuality = float(row[1])
        MeasurementAngle = float(row[2])
        if MeasurementQuality > 11 and float(row[3]) < 5000 and float(row[3]) > 500:
            MeasurementDistance = float(row[3])
            [MeasurementDistanceX, MeasurementDistanceY] = pol2cart(MeasurementDistance,
                                                                    -(2 * math.pi * MeasurementAngle / 360))
        else:
            MeasurementDistance = None
            MeasurementDistanceX = None
            MeasurementDistanceY = None
            
        NewRotations.append(NewRotation)
        MeasurementQualities.append(MeasurementQuality)
        MeasurementAngles.append(MeasurementAngle)
        MeasurementDistances.append(MeasurementDistance)
        MeasurementDistancesX.append(MeasurementDistanceX)
        MeasurementDistancesY.append(MeasurementDistanceY)
    
    # Creation of numpy arrays for parameters.
    NP_MeasurementDistances = np.array(MeasurementDistances)
    NP_MeasurementAngles = np.array(MeasurementAngles)
    
    # Sampling of useful data.
    NP_MeasurementAngles = NP_MeasurementAngles[NP_MeasurementDistances != None]
    NP_MeasurementDistances = NP_MeasurementDistances[NP_MeasurementDistances != None]

    # 1D interpolator.
    interpolator = interp1d(NP_MeasurementAngles, NP_MeasurementDistances)
    
    # Creation of numpy arrays for parameters.
    NP_MeasurementDistances = np.array(MeasurementDistances)
    NP_MeasurementAngles = np.array(MeasurementAngles)
    
    # Creation of numpy array.
    NP_MeasurementDistancesX = np.array(MeasurementDistancesX)
    NP_MeasurementDistancesY = np.array(MeasurementDistancesY)
    
    # Updating the missing data points.
    for i in range(0,np.size(NP_MeasurementAngles)):
        if NP_MeasurementDistances[i] == None:
            NP_MeasurementDistances[i] = interpolator(NP_MeasurementAngles[i])
            [NP_MeasurementDistancesX[i], NP_MeasurementDistancesY[i]] = pol2cart(NP_MeasurementDistances[i],
                                                                    -(2 * math.pi * NP_MeasurementAngles[i] / 360))
# Plotting the data.
fig = go.Figure(data= go.Scatterpolar(r=NP_MeasurementDistances,theta=NP_MeasurementAngles,mode='markers'))

fig.update_layout(showlegend=False)
fig.show()

fig = go.Figure(data=go.Scatter(x=NP_MeasurementDistancesX, y=NP_MeasurementDistancesY, mode='markers'))
fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))
fig.show()

# CSV output file.
with open('xdata.csv', 'w') as f:
    np.savetxt(f, NP_MeasurementDistancesX, delimiter=',')
    
with open('ydata.csv', 'w') as f:
    np.savetxt(f, NP_MeasurementDistancesY, delimiter=',')

    
