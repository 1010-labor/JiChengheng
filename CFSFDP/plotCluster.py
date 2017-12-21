#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
# File Name: plotCluster.py
# Author: Zhang Jiafeng
# mail: 761369516@qq.com
# Created Time: 2017-05-05
from numpy import *
def loadData():
    dataMat = []; labelMat = [];
    fr = open('result.dat')
    for line in fr.readlines():
        lineArr = line.strip().split()
        dataMat.append([float(lineArr[0]), float(lineArr[1])])
        labelMat.append(int(lineArr[2]))
    return dataMat,labelMat

def show():
    import matplotlib.pyplot as plt
    dataMat,labelMat=loadData()
    dataArr = array(dataMat)
    n = shape(dataArr)[0]
    xcord1 = []; ycord1 = []
    xcord2 = []; ycord2 = []
    xcord3 = []; ycord3 = []
    xcord4 = []; ycord4 = []
    for i in range(n):
        if int(labelMat[i]) == 1:
            xcord1.append(dataArr[i,0]); ycord1.append(dataArr[i,1])
        elif int(labelMat[i]) == 2:
            xcord2.append(dataArr[i,0]); ycord2.append(dataArr[i,1])
        elif int(labelMat[i]) == 3:
            xcord3.append(dataArr[i,0]); ycord3.append(dataArr[i,1])
        elif int(labelMat[i]) == 4:
            xcord4.append(dataArr[i,0]); ycord4.append(dataArr[i,1])
    plt.figure(figsize=(14,12))
    plt.scatter(xcord1, ycord1, s=5, c='red', marker='o')
    plt.scatter(xcord2, ycord2, s=5, c='blue', marker='o')
    plt.scatter(xcord3, ycord3, s=5, c='green', marker='o')
    plt.scatter(xcord4, ycord4, s=5, c='yellow', marker='o')
    plt.show()

show()
