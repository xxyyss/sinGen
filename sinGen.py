from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import scipy.fftpack
import cmath
import sys
import random

NOISE_LEVEL_LIST=[30.0,50.0]
OFFSET_VALUE_LIST=[-1,1,5,10]

def sineGen(offset=0.0):
    index=np.linspace(0,4*math.pi,num=100,endpoint=True)
    print 'index length:',len(index)
    print 'offset',offset
    sinWave=[math.sin(item+offset) for item in index]
    max_amplitude = 32767.0
    audioSamples = [int(point * max_amplitude) for point in sinWave]
    return audioSamples

def squareGen(offset=0.0):
    result=[]
    for index in range(100):
        if ((index+offset)>10 and (index+offset)<20 or (index+offset)>80 and (index+offset)<90):
            result.append(1)
        else:
            result.append(0)
    return result

def sineIn0(offset=0.0,hasNoise=False):
    result=[]
    for index in np.arange(0,100,0.25):
        noise=0
        if (hasNoise):
            noise=random.uniform(-1,1)*32767.0*NOISE_LEVEL/100.0
        if ((index-offset)>=30 and (index-offset)<40):
            amplitude=math.sin(20*math.pi/10*(index+offset-30))*32767.0+noise
        elif ((index-offset)>=60 and (index-offset)<70):
            amplitude=math.sin(20*math.pi/10*(index+offset-60))*32767.0+noise
        else:
            amplitude=noise
        result.append(amplitude)
    return result

def plotData(data,pngFileName,auxData=None):
    """plot numpy array data into png file"""
    #separate different axes
    figure = plt.gcf() # get current figure
    figure.set_size_inches(16, 6)
    plt.xticks(range(0,len(data),2))
    plt.plot(data)
    if (auxData != None):
        reverseAuxData=[-value for value in auxData]
        plt.plot(reverseAuxData)
    plt.ylabel('audio data')
    d = os.path.dirname(pngFileName)
    if not os.path.exists(d):
        os.makedirs(d)
    plt.savefig(pngFileName,dpi=100)
    plt.close()
#    print 'finished saving png file',pngFileName

def gccphatequation13(signalList):
    """calc gcc-PHAT value of two signals in the signalList"""
    BETA_LOW=-20
    BETA_HIGH=20
    BETA_RESOLUTION=1
    signalfftList=[]
    #calculate fft
    num=0
    for signal in signalList:
#        print '\nsingal received from equation:\n',signal
        fft=scipy.fftpack.fft(signal)
#        print 'fft value:',fft
        signalfftList.append(fft)
        fftMag=[abs(value) for value in fft]
        fftPhase=[cmath.phase(value) for value in fft]
 #       plotData(fftPhase,'./fftPhase'+str(num)+'.png')
 #       plotData(fftMag,'./fft'+str(num)+'.png')
        num+=1
#    print 'each fftlist length: %d, fftlist: %s\n' %(len(signalfftList[0]),signalfftList)
    equation13ResultList=[]
    #for each beta calculate the sum of cosine
    for eachBeta in np.arange(BETA_LOW,BETA_HIGH,BETA_RESOLUTION):
#        print 'beta:',eachBeta
        result=equation13(signalfftList,eachBeta,100)
        equation13ResultList.append(result)
#    print 'betaResult length: %f, betaResult: %s\n' %(len(equation13ResultList),equation13ResultList)
    resultArray=np.array(equation13ResultList)
    maxValue=np.amax(resultArray)
    maxValIndex=np.argmax(resultArray)
#    plotData(np.array(resultArray),'./resultArray.png')
    bestBeta=maxValIndex*BETA_RESOLUTION+BETA_LOW
    print 'maxVal : %f, maxValIndex:%d, best beta : %f' %(maxValue,maxValIndex,bestBeta)
    return bestBeta

def equation13(signalList,beta,numSample):
    """implementing the sum of cosine in equation 13"""
    result=0.0
    for frequencyIndex in range(numSample):
        #get the angle
        signal1Angle=mapto2piRange(cmath.phase(signalList[0][frequencyIndex]))
#        print 'signal 1 angle: ' , signal1Angle
        signal2Angle=mapto2piRange(cmath.phase(signalList[1][frequencyIndex]))
#        complexDiv=signalList[0][frequencyIndex]/signalList[1][frequencyIndex]
#        signalAngleDiff=mapto2piRange(cmath.phase(complexDiv))
#        print 'signal 2 angle: ' , signal2Angle
        cosResult=math.cos (signal1Angle - signal2Angle - (2.0 * math.pi *beta * frequencyIndex)/numSample)
#        cosResult=math.cos (signalAngleDiff - (2.0 * math.pi *beta * frequencyIndex)/numSample)
#        print 'cosResult: ', cosResult
        result=result+cosResult
#    print 'final result:',result
    return result

def mapto2piRange(angle=0.0):
    """map the return value of cmath.phase() in [-pi,pi] to [0,2*pi] region"""
    result=angle
    if (angle<0):
        result=2*math.pi-abs(angle)
    # print type(angle)
    # print 'angle: %f, result: %f'%(angle,result)
    return result

def test(NOISE_LEVEL,OFFSET_VALUE):
    outputFolder='./'+str(NOISE_LEVEL)+'%noise/'
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
    sys.stdout=open(outputFolder+'offset'+str(OFFSET_VALUE)+'.log','w+')
    totalCnt=0
    correctCnt=0
    for runTimes in np.arange(0,5):
        totalCnt+=1
        sineWave=sineIn0(0.0,True)
        sineWaveOffset=sineIn0(OFFSET_VALUE,True)
    # print 'sine wave:x',sineWave
    # print 'sine wave offset:',sineWaveOffset
        plotData(np.array(sineWave),outputFolder+str(OFFSET_VALUE)+'graph/'+'mSinWave'+str(runTimes)+'.png')
        plotData(np.array(sineWaveOffset),outputFolder+str(OFFSET_VALUE)+'graph/'+'mSinWave'+str(runTimes)+'_offset.png')
    #test for gccphat
        signalList=[sineWave,sineWaveOffset]
#        print 'signal list sent from test:\n',signalList
        diff=gccphatequation13(signalList)
        #diff should give how late the first signal is compare to the second. if diff<0, then the first signal is ahead of the second
        print 'offset: %d, perceived diff: %d'%(OFFSET_VALUE,diff)
        if (OFFSET_VALUE!=diff):
            print 'WRONG'
        else :
            print 'correct'
            correctCnt+=1
    print 'total cnt: %d, correct cnt: %d, ratio:%f'%(totalCnt,correctCnt,correctCnt/totalCnt)

if __name__ == "__main__":
#    sys.stdout = open('./testGCCPHAT.log', 'w')
#    sineWave=sineGen()
#    sineWaveOffset=sineGen(4*math.pi/100*10)
#    sineWave=squareGen()
#    sineWaveOffset=squareGen(10)
    for NOISE_LEVEL in NOISE_LEVEL_LIST:
        for OFFSET in OFFSET_VALUE_LIST:
            test(NOISE_LEVEL,OFFSET)

