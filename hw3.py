import scipy.io as sio
import numpy as np
import re
from scipy import signal
from scipy.ndimage import imread
import matplotlib.pyplot as plt
import cmath
import math
import time
from PIL import Image

#=========================  WAVELET FUNCTIONS  ===================================

#define the morlet function that return the real part
def morlet_real(x, y, sig, theta, C1, C2):
    # set variables
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    pie = np.pi
    # one peak morlet wave function with greek letter equal to 4
    exponentOfEInsideBrackets = (pie / (2 * sig)) * ((x * cosTheta) + (y * sinTheta))
    exponentOfEOutsideBrackets = -(x**2 + y**2)/ (2 * sig**2)

    #morlet wave function
    #cmath.rect(r, phi) Return the complex number x with polar coordinates r and phi.
    z = C1 / sig * (cmath.rect(1, exponentOfEInsideBrackets) - C2) * np.exp(exponentOfEOutsideBrackets)
    return z.real


#define the morlet function that return the imaginary part
def morlet_imag(x, y, sig, theta, C1, C2):
    # set variables
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    pie = np.pi
    # one peak morlet wave function with greek letter equal to 4
    exponentOfEInsideBrackets = (pie / (2 * sig)) * ((x * cosTheta) + (y * sinTheta))
    exponentOfEOutsideBrackets = -(x**2 + y**2)/ (2 * sig**2)

    #morlet wave function
    #cmath.rect(r, phi) Return the complex number x with polar coordinates r and phi.
    z = C1 / sig * (cmath.rect(1, exponentOfEInsideBrackets) - C2) * np.exp(exponentOfEOutsideBrackets)
    return z.imag


#finds the constants c2
def find_c2(xymin, xymax, sig, theta):
    numerator = 0
    denominator = 0
    cosine = np.cos
    cosineTheta = np.cos(theta)
    sineTheta = np.sin(theta)
    pie = np.pi
    for x in range(xymin, xymax+1, 1):
        for y in range( xymin, xymax+1, 1):
            numerator = numerator + (cosine((pie / (2 * sig)) * ((x * cosineTheta) + (y * sineTheta))) * np.exp(-(x**2 + y**2)/(2 * sig**2)))
            denominator = denominator + (np.exp(-(x**2 + y**2)/(2 * sig**2)))

    C2 = numerator/denominator
    return C2


#finds the constant c1
def find_c1(xymin, xymax, sig, theta, C2):
    Z = 0
    pie = np.pi
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    cosine = np.cos

    for x in range(xymin, xymax+1, 1):
        for y in range( xymin, xymax+1, 1):
            Z = Z + (1 - 2* C2 * cosine(pie/(2*sig) * ((x * cosTheta) + (y * sinTheta))) + C2**2) * np.exp((-(x**2 + y**2)/sig**2))
    C1 = 1/np.sqrt(Z)

    return C1


#plot the morlet function for the real
def morletMatrix_real(xymin, xymax, sig, theta):

    #find c1 and c2
    C2 = find_c2(xymin, xymax, sig, theta)
    C1 = find_c1(xymin, xymax, sig, theta, C2)

    #define grid over which the function should be plotted
    xx, yy = np.meshgrid(np.linspace(xymin, xymax, 33),np.linspace(xymin, xymax, 33))

    # fill a matrix with the morlet function values
    zz= np.zeros(xx.shape)
    for i in range(yy.shape[0]):
        for j in range(xx.shape[0]):
            zz[i,j] = morlet_real(xx[i,j], yy[i,j], sig, theta, C1, C2)

    return zz

# plot morlet function for imiginary
def morletMatrix_imag(xymin, xymax, sig, theta):
    #determine constants
    C2 = find_c2(xymin, xymax, sig, theta)
    C1 = find_c1(xymin, xymax, sig, theta, C2)

    #define grid over which the function should be plotted
    xx = np.meshgrid(np.linspace(xymin, xymax, xymax-xymin+1))
    yy = np.meshgrid(np.linspace(xymin, xymax, xymax-xymin+1))

    # fill a matrix with the morlet function values
    zz= np.zeros(xx.shape)
    for i in range(yy.shape[0]):
        for j in range(xx.shape[0]):
            zz[i,j] = morlet_imag(xx[i,j], yy[i,j], sig, theta, C1, C2)

    return zz



#========================== LOAD FILE ==============================================

#load mathlab file
imageFile = sio.loadmat('I013.mat')
# get info about symmetry axis for 1 line of symmetry
contents = imageFile['segments']
print ("\n =================================================\n")
print("\nCoordinates of starting point of symmetry are "+ str(contents[0][0][0]) + " and end points are " + str(contents[0][0][1]))
print ("\n =================================================\n")

image= imread("i013.png", True)

# dimensions of image that we are minapulating
column = len(image[0])
imageRow = len(image)

threshold = 100
valuesAboveThreshold = [[0]*column for i in range(imageRow)]# this array will hold all of the wavelet responses above threshold, for each of the 4 angles

print ("\nStarting wavelet computations")
start_time = time.clock()
for i in range(len(valuesAboveThreshold)):
    for j in range(len(valuesAboveThreshold[0])):
        arr = []
        valuesAboveThreshold[i][j] = arr

# # compute wavelets at different angles
kernel = morletMatrix_real(-16, 16, 6, 0)
waveletZero = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(waveletZero):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)

kernel = morletMatrix_real(-16, 16, 6, np.pi/6)
waveletQuaterPie = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(waveletQuaterPie):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)

kernel = morletMatrix_real(-16, 16, 6, 2*np.pi/6)
waveletHalfPie = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(waveletHalfPie):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)

kernel = morletMatrix_real(-16, 16, 6, 3*np.pi/6)
wavelet3Pie = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(wavelet3Pie):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)

kernel = morletMatrix_real(-16, 16, 6, 4*np.pi/6)
wavelet3Pie = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(wavelet3Pie):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)

kernel = morletMatrix_real(-16, 16, 6, 5*np.pi/6)
wavelet3Pie = signal.convolve2d(image, kernel, boundary='symm', mode='same')
for (x,y), value in np.ndenumerate(wavelet3Pie):
    if(value > threshold):
        valuesAboveThreshold[x][y].append(value)
    else:
        valuesAboveThreshold[x][y].append(0)



print ("we now have wavelet responses for the 6 angles that pass the threshold for each pixel stored in an array")
print "Time taken was " + str(time.clock() - start_time), "seconds\n"
print ("\n =================================================\n")
#===================================================================================


symmetryAccumulatorRow = 0
symmetryAccumulatorPhi =0
iterArr = [-(np.pi), 0, np.pi/4]

coordinatesAngles = []

# values above threshold is a 2x2 array containing x,y,tuppleOfSixAngles indicies
for i in range(len(valuesAboveThreshold)):
    for j in range(len(valuesAboveThreshold[0])):
        for k in range(len(valuesAboveThreshold[i][j])):
            # if any 1 angle of the 6 are above the threshold then we add it to the coordinateAngles array (which is a 1d array)
            if valuesAboveThreshold[i][j][k] != 0:
                xyPosition = (i,j, valuesAboveThreshold[i][j])
                coordinatesAngles.append(xyPosition)
                break
# debug
# for i in range(len(valuesAboveThreshold)):
#     for j in range(len(valuesAboveThreshold[0])):
#         print valuesAboveThreshold[i][j]
print("We now have an array with each index containing the 6 wavelet responses\nStartin rho min and max computation")
start_time2 = time.clock()

#==================================== CALCULATE ROW MIN AND ROW MAX =================================
rowMax = -100000010000000
rowMin = 100000000000000
approxDirectionsms = [round(0,7), round(np.pi/6,7), round(2*np.pi/6,7),round(3*np.pi/6,7),round(4*np.pi/6,7),round(5*np.pi/6,7)]
approxDirections = [0, np.pi/6, 2*np.pi/6,3*np.pi/6,4*np.pi/6,5*np.pi/6]
for i in range(len(coordinatesAngles)):
    x1 = coordinatesAngles[i][0]# x coordinate
    y1 = coordinatesAngles[i][1]# y coordinate
    for j in range(i+1,len(coordinatesAngles)):
        x2 = coordinatesAngles[j][0]# x' coordinate
        y2 = coordinatesAngles[j][1]# y' coordinate
        #EXTRACT

        #calculate phi
        den = x2-x1

        if den == 0:
            phi = 0
        else:
            phi = (np.pi/2) + np.arctan((y2-y1)/(x2-x1))
        # check the 6 angles to find the closest match to pi
        myArray = np.array(approxDirections)
        pos = (np.abs(myArray-phi)).argmin()
        phi = pos# replace phi with closest angle that we find

        # calculate c
        cx = (x1 + x2) * 0.5
        cy = (y1 + y2) * 0.5

        # calculate row
        row = (cx * math.cos(phi - np.pi/2)) +(cy* math.sin(phi - np.pi/2))

        if(row < rowMin):
            rowMin = row
        if(row > rowMax):
            rowMax = row
print("Rho min is " + str(rowMin))
print("Rho max is " + str(rowMax))
print "Time taken was " + str((time.clock() - start_time2)/60), "minutes\n"
print ("\n =================================================\n")

#================================================ ALGORITHM ======================================

start_time3 = time.clock()
MSSymmetryMatrix = [[0]*column for i in range(imageRow)]# this will be the same size as the image and will store the MS at whichever pixel has symmetry
print("Starting algorithm")
# for each index I' in array, select I' x and y values
for i in range(len(coordinatesAngles)):
    x1 = coordinatesAngles[i][0]# x coordinate
    y1 = coordinatesAngles[i][1]# y coordinate

    # for each angle after that, I'+1 till the end of the array,
    # get the x and y value for each of the indexes
    for j in range(i+1,len(coordinatesAngles)):
        x2 = coordinatesAngles[j][0]# x' coordinate
        y2 = coordinatesAngles[j][1]# y' coordinate

        #EXTRACT

        #calculate phi
        # chech to see if the denominator will equal to 0, if it is set phi to 0 and skip computation
        den = x2-x1
        if den == 0:
            phi = 0
        else:
            phi = (np.pi/2) + np.arctan((y2-y1)/(x2-x1))# this is the angle that is between the 2 pixels

        # check the 6 angles to find the closest match to pi
        myArray = np.array(approxDirections)
        pos = (np.abs(myArray-phi)).argmin()
        phi = pos# replace phi with closest angle that we find

        # calculate c
        cx = (x1 + x2) * 0.5
        cy = (y1 + y2) * 0.5

        # calculate rho
        row = (cx * math.cos(phi - np.pi/2)) +(cy* math.sin(phi - np.pi/2))

        # discretize rho
        rowIndex = math.floor(((row-rowMin)/(rowMax-rowMin))*99)


        # we are looking for contrast along these three angles. we are looking
        # to see if the 3 responses match with the values we have in order to detect smmetry
        threeAngleResponses = [-np.pi/3, 0, np.pi/3]
        for k in range(len(threeAngleResponses)):

            # this angle may be more or less than the range we want, so we subtract and add pi accordingly to balance it off
            thetai = phi - np.pi/2 + threeAngleResponses[k] # this is x's angle
            if(thetai > np.pi):
                thetai = thetai - np.pi
            if(thetai < 0):
                thetai = thetai + np.pi
            # print thetai

            thetaj = 2*phi - thetai# this will be x' angle
            if(thetaj > np.pi):
                thetaj = thetaj - np.pi
            if(thetaj < 0):
                thetaj = thetaj + np.pi
            # print thetaj
            # thetas will equal to one of these 6 indicies so we choose the index number(0-5) in order to do the next calculation
            for z in range(len(approxDirectionsms)):
                if (round(thetai,7) == approxDirectionsms[z]):
                    indexi = z
                if (round(thetaj,7) == approxDirectionsms[z]):
                    indexj = z
            # check if these two angles are good responses with respect to our measure of symmetry
            numerator = ((valuesAboveThreshold[x1][y1][indexi]) * np.conjugate((valuesAboveThreshold[x2][y2][indexj]))) + (np.conjugate((valuesAboveThreshold[x1][y1][indexi])) * (valuesAboveThreshold[x2][y2][indexj]))
            denominator = (abs(valuesAboveThreshold[x1][y1][indexi]))**2 + (abs(valuesAboveThreshold[x1][y1][indexj]))**2

            if denominator == 0:
                MS = 0
            else:
                complexMs = numerator/denominator
                MS = complexMs.real
                print("got MS " + str(MS) + " and storing it")

            newindx = int(rowIndex)
            newindy = int(phi)
            # debug
            # print "rho index into the symmetry array is "+ str(newindx)
            # print "phi index into the symmetry array is "+ str(newindy)

            #put caluclation at the correct index which will by at x = rho and y = phi
            MSSymmetryMatrix[newindx][newindy] = MS
print "Time taken was " + str((time.clock() - start_time3)/60), "minutes\n"
print("\nDone with Symmetry accumulator")
print ("\n =================================================\n")


#============================= PLOT IMAGE ==================================
print("Beginning image plotting")

img = Image.open("i013.png")
width, height = img.size
pixels = img.load()

for x in range(width):
        for y in range(height):
            if(MSSymmetryMatrix[x][y] != 0):
                pixels[x,y] = (255,0,0)

img.show()


# def draw(ind):
#     phi = angles[ind[0]]
#     rho = rhos[ind[1]]
#     slope = np.tan(phi)
#     b = rho*np.sin(phi-np.pi/2)
#     x = 0
#     while(x<len(image1)):
#         y = int(np.around(slope*x+b))
#         if(y<0 or y>= len(image1[0])):
#             x+=1
#             continue
#         image1[x][y] = (0.5,0.5,0.5)
#         x+=1


