# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 08:49:26 2019

@author: herism
"""

import time
import os
import numpy as np
import rasterio
import pandas as pd
import math as math
import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import *
from glob import glob
from rasterio import Affine as A
from rasterio.warp import calculate_default_transform, reproject, Resampling
import h5py


# set the timer
t1=time.time()

# landsat folder address
address =r'D:\EU_JRC\Layers\LandsatRome\LC08_L1TP_191031_20170823_20170912_01_T1.tar'
# set the coordinate system to USA_Contiguous_Albers_Equal_Area_Conic (https://epsg.io/102003)

# set the emissivity multiplier here
emissivity=0.96
# under this temperature will set as nan (Degree C)
tempLevel=10
# clear sky
# https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.where.html
lstClearSkyQA=[2720, 2724, 2728, 2732]
np.seterr(divide='ignore', invalid='ignore')


lstErrors=[]

def cal_TOA_Sp_Raidance(Band, Rai_Mu_B, Rai_add_B):
    return (Rai_Mu_B*Band)+Rai_add_B

# TOA planetary reflectance with correction for sun elevation (angel)
def cal_TOA_pl_ref(Band, Ref_Mu_B, Ref_add_B, SunElev):
    Toa_withouthSunCorrection=(Ref_Mu_B*Band)+Ref_add_B
    return (Toa_withouthSunCorrection/(math.sin(SunElev)))
# this function calculates brightness temperature    
def calStemperature (TOA_Sp_Ra_band10,K1ConstandBand10,K2ConstandBand10,emissivityAr):
    #d0=np.true_divide(K1_CONSTANT_BAND_10,TOA,where=TOA>0)
    d0=K1ConstandBand10/(TOA_Sp_Ra_band10+1)
    d1=(np.log(d0))
    #BT= (np.true_divide(K2_CONSTANT_BAND_10,d1,where=d1!=0))-273.15
    BT= (K2ConstandBand10/d1)-273.15
    #c0=np.true_divide(BT,1.4388)
    c0=BT/1.4388
    c1=(1 + (0.00115 * c0) * math.log(emissivityAr))
    #ST=np.true_divide(BT,c1)
    ST=BT/c1
    return ST

# this function finds the dataset (level) that is empty
## gets the hdf file and the window to get the data
def findEmpty (f,rowStart,rowEnd,colStart,colEnd):
    maxVal = -1
    lstKey=[]
    for key in f.items():
        lstKey.append (key[0])
    global count1
    count1=0
    while maxVal !=0 & count1<len(lstKey):
        maxVal=(f[lstKey[count1]][rowStart:rowEnd,colStart:colEnd]).max()
        count1+=1
    return lstKey[count1]

# this function gets a mother raster and a child raster
   # we want to write the child raster in the mother raster
   # to make sure that the bounaries fit and we find the overlap
   # we use this function to find the overlapping boundaries
def calculatableBounds2(mainRstBounds,fillingRstBounds):
    # match the top side
    if fillingRstBounds[3]>mainRstBounds[3]:
        endY = mainRstBounds[3]
    else:
        endY = fillingRstBounds[3]
    # match the bottom side  
    if mainRstBounds[1]>fillingRstBounds[1]:
        startY = mainRstBounds[1]
    else:
        startY = fillingRstBounds[1]
    # match the left side
    if mainRstBounds[0]>fillingRstBounds[0]:
        startX = mainRstBounds[0]
    else:
        startX = fillingRstBounds[0]
    # match right side
    if fillingRstBounds[2]>mainRstBounds[2]:
        endX = mainRstBounds[2]
    else:
        endX = fillingRstBounds[2] 
    return (startX,startY,endX,endY)


# this function recieves the calculated boundary and rasters
        # then calculates the start and end rows and cols that will be used to fill the arrays
def sliceFillingRst (calBounds,motherRstBounds,motherRstShape, cellSizeAr):
    # calBounds is (startX,startY,endX,endY)
    startRow = int((motherRstBounds[3] - calBounds[3])/cellSizeAr)
    endRow   = int(motherRstShape[0] - ((calBounds[1] - motherRstBounds[1])/cellSizeAr))
    startCol = int((calBounds[0] - motherRstBounds[0])/cellSizeAr)
    endCol   = int(motherRstShape[1] - ((motherRstBounds[2]-calBounds[2])/cellSizeAr))
    return (startRow,endRow,startCol,endCol)


def getMetaValues (mlt_file):

    with open(mlt_file,'r') as txtFile:
        lines = txtFile.readlines()
    dic={}
    for li in lines:
        if 'K1_CONSTANT_BAND_10 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['k1_constant_band_10']=float(n)
            
        if 'K2_CONSTANT_BAND_10 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['k2_constant_band_10']=float(n)
            
        if 'K1_CONSTANT_BAND_11 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['k1_constant_band_11']=float(n)
            
        if 'K2_CONSTANT_BAND_11 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['k2_constant_band_11']=float(n)
            
        if 'RADIANCE_MULT_BAND_10 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Radiance_mult_band_10']=float(n)
            
        if 'REFLECTANCE_MULT_BAND_4 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Reflectance_mult_band_4']=float(n)
            
        if 'REFLECTANCE_MULT_BAND_5 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Reflectance_mult_band_5']=float(n)
            
        if 'REFLECTANCE_ADD_BAND_4 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Reflectance_add_band_4']=float(n)   
            
        if 'REFLECTANCE_ADD_BAND_5 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Reflectance_add_band_5']=float(n)
            
        if 'RADIANCE_ADD_BAND_10 = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Radiance_add_band_10']=float(n)
            
        if 'SUN_ELEVATION = ' in li:
            locEsign=li.find('=')
            n=li[locEsign+1:-1]
            dic['Sun_Elevation']=float(n)
            
        if 'LANDSAT_PRODUCT_ID = ' in li:
            locEsign=li.find('"')
            p=li[locEsign+1:-2]
            dic['Landsat_Product_Id']=(p)
    return dic

mainRasterAddresstest=r'D:\EU_JRC\Layers\PilotFiles\Temp\lc_30.tif'
with rasterio.open(mainRasterAddresstest) as dst:
    kwdsMainRaster = dst.meta.copy()
    dst_crs=kwdsMainRaster['crs']

#directories=[]
#for d in os.walk(address):
#    directories.append(d[0])
#
##directories=directories[137:]
#
#lstErrorAddresses=[]
#for f in range(1,len(directories)):
#    folderAddress = (directories[f])
#    print ('number ', f, ' out of 455; started working on folder: ',folderAddress)

folderAddress =r'D:\EU_JRC\Layers\LandsatRome\LC08_L1TP_191031_20170823_20170912_01_T1.tar'
# here we get the bands of each image. In one folder, there are more than one scene and more than one band for each
band_4s=(glob(folderAddress+'\\*B4.TIF'))
band_5s=(glob(folderAddress+'\\*B5.TIF'))
band_10s=(glob(folderAddress+'\\*B10.TIF'))
band_MTLs=(glob(folderAddress+'\\*MTL.TXT'))
band_BQAs=(glob(folderAddress+'\\*BQA.TIF'))

# create a list for the arrays of the satellite images
# we also hold the info about each scene
#try:
# the loop goes to each scene and does all calcualtions
for sc in range(0,len(band_4s)):
    #print ('sarted working on image ',band_10s[sc])
    b4 = band_4s[sc]
    b5 = band_5s[sc]
    b10= band_10s[sc]
    MLT_text_file=band_MTLs[sc]
    bQA=band_BQAs[sc]
    
    with rasterio.open(b10,'r') as b10Rst:
        arb10=b10Rst.read(1)
        transform, width, height = calculate_default_transform(
                b10Rst.crs, dst_crs, b10Rst.width, b10Rst.height,
                resolution=(30,30), *b10Rst.bounds)
        sourceCRS=b10Rst.crs
        sourceTransofmr=b10Rst.transform
        kwargs = b10Rst.meta.copy()
        kwargs.update({
            'dtype': 'uint16',
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height})

        
        
    with rasterio.open(b4) as b4Rst:
        arb4=b4Rst.read(1)
    
    with rasterio.open(b5,'r') as b5Rst:
        arb5=b5Rst.read(1)
    
    # find the values of quality assessments
    # https://landsat.usgs.gov/collectionqualityband
    with rasterio.open(bQA,'r') as bQArst:
        arbQA=bQArst.read(1)
    
        
    #plt.imshow(arb10)
    ## here we get all the constant values out of the metadata file
    metaDic=getMetaValues (MLT_text_file)
    # K1_CONSTANT_BAND_10
    K1ConsBand10 = metaDic['k1_constant_band_10']
    # K2_CONSTANT_BAND_10
    K2ConstantBand10 =metaDic['k2_constant_band_10']
    
    # RADIANCE_MULT_BAND_4= float(lines[168][27:-1])
    # RADIANCE_MULT_BAND_5= float(lines[169][27:-1])
    # RADIANCE_MULT_BAND_10
    RadianceMultBand10=metaDic['Radiance_mult_band_10']
    # REFLECTANCE_MULT_BAND_4
    ReflectanceMultBand4 = metaDic['Reflectance_mult_band_4']
    # REFLECTANCE_MULT_BAND_5
    ReflectanceMultBand5 = metaDic['Reflectance_mult_band_5']
    
    # REFLECTANCE_ADD_BAND_4
    ReflectanceAddBand4 = metaDic['Reflectance_add_band_4']
    # REFLECTANCE_ADD_BAND_5
    ReflectanceAddBand5 = metaDic['Reflectance_add_band_5']
    
    # RADIANCE_ADD_BAND_4 = float(lines[179][26:-1])
    # RADIANCE_ADD_BAND_5 = float(lines[180][26:-1])
    # RADIANCE_ADD_BAND_10
    RainanceAddBand10= metaDic['Radiance_add_band_10']
    # SUN_ELEVATION
    SunElevation = metaDic['Sun_Elevation']
    
    LandsatProductID= metaDic['Landsat_Product_Id']
    
    
    # https://www.usgs.gov/land-resources/nli/landsat/using-usgs-landsat-level-1-data-product
    # TOA spectral radiance (Watts/( m2 * srad * μm))
                
    # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.where.html
    # find the clearSky Values in the BQA array
    clearSkyAr = np.isin(arbQA,lstClearSkyQA)
    notClearSkyAr = np.invert(clearSkyAr)
    
    #TOA_Sp_Ra_b10 = RADIANCE_MULT_BAND_10 * arb10 + RADIANCE_ADD_BAND_10
    TOA_Sp_Ra_b10 =cal_TOA_Sp_Raidance(arb10, RadianceMultBand10, RainanceAddBand10)
    # TOA planetary reflectance with correction for sun elevation (angel)
    TOA_pl_ref_b4 = cal_TOA_pl_ref(arb4,ReflectanceMultBand4, ReflectanceAddBand4, SunElevation)
    TOA_pl_ref_b5 = cal_TOA_pl_ref(arb5,ReflectanceMultBand5, ReflectanceAddBand5, SunElevation)
    
    # calculate NDVI
    # Allow division by zero
    np.seterr(divide='ignore', invalid='ignore')
    # Calculate NDVI
    arNDVI = (((TOA_pl_ref_b5.astype(float) - TOA_pl_ref_b4.astype(float)) / (TOA_pl_ref_b5 + TOA_pl_ref_b4))*100).astype('int16')
    
    # calculate surface temperature
    arST=((calStemperature(TOA_Sp_Ra_b10,K1ConsBand10,K2ConstantBand10,emissivity))*100).astype('int16')
    
    # set temperatures under tempeLevel in the mask
    #notClearSkyAr[arST<10]=True
    # set nan for not clear sky in the ndvi array
    arNDVI[notClearSkyAr]=-9999
    # set nan for not clear sky in the ndvi array
    arST[notClearSkyAr]=-9999
    #print (arST.shape)
    # add the arrays to the lists
    stAdr   = folderAddress+'\\'+LandsatProductID+'_ST.tif'  
    ndviAdr = folderAddress+'\\'+LandsatProductID+'_NDVI.tif'      
    
    with rasterio.open(stAdr, 'w', **kwargs) as dst:
        reproject(
            source=arST,
            destination=rasterio.band(dst, 1),
            src_transform=sourceTransofmr,
            src_crs=sourceCRS,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest)
    
    
    # write the arrays in rasters       
    
    with rasterio.open(ndviAdr, 'w', **kwargs) as dst:
        reproject(
            source=arNDVI,
            destination=rasterio.band(dst, 1),
            src_transform=sourceTransofmr,
            src_crs=sourceCRS,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest)    # if there was an error running the scene
#    #except:
#        print ('I got an aerror doing analysis for this address: ', folderAddress)
#        lstErrors.append(folderAddress)
#t2=time.time()
#print ('time passed: ',t2-t1)




