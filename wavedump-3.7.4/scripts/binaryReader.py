#!/usr/bin/python
# specialReader.py

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from binaryReaderUtilities  import peakdet
from binaryReaderUtilities  import getData, getStats

def rotate(l, n):
    return l[-n:] + l[:-n]

def myWavePlot( id, aTitle, times, waves, waveNumb, nev, voffset ) :
	if id < 0 : return
	import matplotlib.pyplot as pylab
        #print "length of times is ", len(times)," length of waves is ", len(waves)
	#print "length of waves[0] is ", len(waves[0])
        pylab.figure(id)
        P=pylab.plot(times,waves[waveNumb+0], linewidth=1.0)
        #P=P+pylab.plot(timesAtThreshold[0],threshold, marker='o')
	for iw in range(1,nev) :
	   wave=waves[waveNumb+iw]
	   for i in range(len(wave)) : wave[i] = wave[i] + voffset*iw
           P=P+pylab.plot(times,wave, linewidth=1.0)
           #P=P+pylab.plot(timesAtThreshold[1],threshold, marker='o')
	#wave=waves[waveNumb+2]
	#for i in range(len(wave)) : wave[i] = wave[i] + 20.0
        #P=P+pylab.plot(times,wave, linewidth=1.0)
        #P=P+pylab.plot(timesAtThreshold[2],threshold, marker='o')

        xl=pylab.xlabel('time (sec)')
        yl=pylab.ylabel('voltage (mv)')
        [xlow,xhgh,ylow,yhgh]=pylab.axis()
        xlow=130.0;xhgh=170.0
        xlow=0.0;xhgh=200.0
        pylab.axis([xlow, xhgh, ylow, yhgh])
	#aTitle=aTitle+" thres=%f mv"%threshold
        #aTitle=" waveforms %d %d %d" % ( waveNumb, waveNumb+1, waveNumb+2)
	bTitle = ": %d events starting at %d" % ( nev, waveNumb ) 
        aTitle=pylab.title(aTitle+bTitle)
        pylab.grid(True)
        #pylab.savefig("timeAtThres.png")
        #pylab.show(P)

def bincenters(binedges) :
	centers=[]
	for i in range( len(binedges)-1 ) :
		centers.append( (binedges[i]+binedges[i+1])/2 )
	return centers
def generateerrorbars(values) :
	err=[]
	for i in range( len(values) ) :
		err.append( math.sqrt( values[i] ) )
	return err

def printTypeOf(objDesc,  object ) :
	print "Type of ", objDesc, "is", type(object) 

run2V={'Run1035':2.7, 'Run1036':2.9, 'Run1037':3.1, 
       'Run1038':3.3, 'Run1039':3.5, 'Run1040':3.9, 
	'Run2001':2.8, 'Run2002':2.8, 'Run2003':2.8 ,
	'CAENData2001':2.8}
f1=open("binaryReader.dat", 'a')

if len(sys.argv) == 3 :
   numbEvents = int (sys.argv[1])
   pathName   = sys.argv[2]
else :
   print "Usage: ",sys.argv[0]," numbEvts pathToData"
   sys.exit()


wavesW0=[]
times=[]
for i in range(1024) : times.append(0.2*i)
wantNevs=numbEvents
print "Reading wave_0"
wavesW0=getData(pathName+"/wave_0.dat", wantNevs)
print "Reading wave_2"
wavesW2=getData(pathName+"/wave_2.dat", wantNevs)
print "Reading wave_4"
wavesW4=getData(pathName+"/wave_4.dat", wantNevs)
print "Reading wave_6"
wavesW6=getData(pathName+"/wave_6.dat", wantNevs)
print "Reading wave_Tr"
wavesTr=getData(pathName+"/TR_0_0.dat", wantNevs)
print "Done Reading wave data"
waveNumb = 0
nev = 10
voffset = 20
"""
myWavePlot( 1, "raw W0", times, wavesW0, waveNumb, nev, voffset )
myWavePlot( 2, "raw W2", times, wavesW2, waveNumb, nev, voffset )
myWavePlot( 3, "raw W4", times, wavesW4, waveNumb, nev, voffset )
myWavePlot( 4, "raw W6",  times, wavesW6, waveNumb, nev, voffset )
myWavePlot( 5, "raw TR0", times, wavesTr, waveNumb, nev, 100.0 )
"""


wavesTrrot=[]
wavesW0rot=[]
wavesW2rot=[]
wavesW4rot=[]
wavesW6rot=[]
thresholdTr=-100.0	# mV
thresholdWs=-4.0	# mV
peakBin=615		# time bin to align to
W0timebins=[]; W0phs=[]
W2timebins=[]; W2phs=[]
W0mW2timebins=[]
W6timebins=[]; W6phs=[]
W0peakCnts=[]
W0peakCharges=[]
for iev in range(len(wavesTr)) :
	if iev%100 == 0 : print "Working on event %d" % iev
	awaveTr = wavesTr[iev]
	maxtabTr, mintabTr = peakdet(awaveTr,-thresholdTr)
	peakLoc = int(mintabTr[0][0])
        binOffset = peakLoc - peakBin
	awaveTrrot = rotate( awaveTr, -binOffset )
	wavesTrrot.append( awaveTrrot )
	awaveW0 = wavesW0[iev]
	awaveW0rot = rotate( awaveW0, -binOffset )
	wavesW0rot.append(awaveW0rot)
	maxtabW0, mintabW0 = peakdet(awaveW0rot,-thresholdWs)
	W0peakCnts.append(len(mintabW0))
	if iev < 10 :
		print "iev,len(mintabW0) = ", iev,len(mintabW0)
        	if len(mintabW0)>0 : print iev, mintabW0[0]
	if len(mintabW0)>0 :
		W0timebins.append(0.200*mintabW0[0][0])
		W0phs.append(-mintabW0[0][1])
	else :
		W0timebins.append(0)
		W0phs.append(0.0)
	awaveW2 = wavesW2[iev]
	awaveW2rot = rotate( awaveW2, -binOffset )
	wavesW2rot.append(awaveW2rot)
	maxtabW2, mintabW2 = peakdet(awaveW2rot,-thresholdWs)
	if len(mintabW0)>0 and len(mintabW2)>0 :
		W0mW2timebins.append(0.200*(mintabW0[0][0]-mintabW2[0][0]))
	else :
		W0mW2timebins.append(-10)
	awaveW4 = wavesW4[iev]
	awaveW4rot = rotate( awaveW4, -binOffset )
	wavesW4rot.append(awaveW4rot)
	awaveW6 = wavesW6[iev]
	awaveW6rot = rotate( awaveW6, -binOffset )
	wavesW6rot.append(awaveW6rot)
	maxtabW6, mintabW6 = peakdet(awaveW6rot,-thresholdWs)
	if len(mintabW6)>0 :
	   W6timebins.append(0.200*mintabW6[0][0]); W6phs.append(-mintabW6[0][1])
	else :
	   W6timebins.append(0); W6phs.append(0.0)
	if len( mintabW0 ) > 0 :
	   W0peakCharge=0
	   ipeak=0
	   ipeakG80 =  0
	   while ipeak < len(mintabW0) :
		if 0.2*mintabW0[ipeak][0] > 80.0 : 
		   ipeakG80 += 1 
		   break
		ipeak += 1
	   #print 'len(mintabW0),ipeak,ipeakG80 = ',len(mintabW0),ipeak,ipeakG80
	   if ipeakG80 > 0 :
	      #print "iplow, iphgh = ",int(mintabW0[ipeak][0])-15, int(mintabW0[ipeak][0])+26
              iplow = max(0, int(mintabW0[ipeak][0])-15)
	      iphgh = min(int(mintabW0[ipeak][0])+26, 1023)
	      for ip in range(iplow, iphgh ) :
	         W0peakCharge += 0.2 * awaveW0rot[ip]	#nsec*mV
	      W0peakCharges.append(W0peakCharge)
	   else :
	      W0peakCharges.append(0)
	else :
	   W0peakCharges.append(0)
"""
	if iev < 100 :
	   print 'ev ipeak mintabW0[ipeak][0], mintabW0[ipeak][0], W0peakCharge'
	   print iev, ipeak, mintabW0[ipeak][0], mintabW0[ipeak][0], W0peakCharge
"""
myWavePlot(  11, 'waves W0 rotated', times, wavesW0rot, waveNumb, nev, voffset )
myWavePlot(  12, 'waves W2 rotated', times, wavesW2rot, waveNumb, nev, voffset )
myWavePlot(  13, 'waves W4 rotated', times, wavesW4rot, waveNumb, nev, voffset )
myWavePlot(  14, 'waves W6 rotated', times, wavesW6rot, waveNumb, nev, voffset )
myWavePlot(  15, 'waves Tr rotated', times, wavesTrrot, waveNumb, nev, 100.0 )




#aInd = np.argmin(wave)
#print 'minimum of wave is ', np.min(wave), ' at index = ', aInd, ' time = ', 0.2*aInd
#minTimes.append( 0.2*aInd )
#waves.append( wave )


plt.figure(101)
h1=np.histogram(W0peakCnts, bins=200)
values = h1[0]
centers=bincenters(h1[1])
errs=generateerrorbars(values)
plt.errorbar(centers, values, yerr=errs, fmt='o')
plt.title("W0 number of peaks"  )
plt.xlabel('number of peaks')
meanNumPksW0=np.mean(W0peakCnts)

plt.figure(102)
h2=np.histogram(W0timebins, bins=200)
values = h2[0]
centers=bincenters(h2[1])
print 'delta time bin (ns)is ', (centers[1]-centers[0])*0.200
errs=generateerrorbars(values)
plt.errorbar(centers, values, yerr=errs, fmt='o')
(W0lenght, W0min, W0max, W0mean, W0std)=getStats( 'Time W0', W0timebins, 50.0, 150.0 )
plt.title("W0 time  est rms=%0.3f (ns)"%W0std)
plt.xlabel('Time (ns)')
getStats( 'Time W6', W6timebins, 50.0, 150.0 )

plt.figure(103)
bins=np.linspace(-15,+15,300)
h3=np.histogram(W0mW2timebins, bins)
values = h3[0]
centers=bincenters(h3[1])
print 'delta time bin (ns)is ', (centers[1]-centers[0])*0.200
errs=generateerrorbars(values)
plt.errorbar(centers, values, yerr=errs, fmt='o')
(W0mW2lenght, W0mW2min, W0mW2max, W0mW2mean, W0mW2std)=getStats( 'Time W0-W2', W0mW2timebins, 50.0, 150.0 )
plt.title("W0-W2 time  est rms=%0.3f (ns)"%W0std)
plt.xlabel('Time (ns)')
getStats( 'Time W6', W6timebins, 50.0, 150.0 )

plt.figure(104)
h4=np.histogram(W6phs, bins=100)
values = h4[0]
centers=bincenters(h4[1])
print 'delta ph bin is ', centers[1]-centers[0]
errs=generateerrorbars(values)
plt.errorbar(centers, values, yerr=errs, fmt='o')
getStats( 'ph W0', W0phs, 0.0, 150.0 )
meanPhW0 = np.mean(W0phs)
npeW0 = ( np.mean(W0phs)/np.std(W0phs) )**2
getStats( 'ph W6', W6phs, 0.0, 150.0 )
npeW6 = ( np.mean(W6phs)/np.std(W6phs) )**2
print " npeW0, npeW6 = ", npeW0, npeW6
plt.title("W0 pulse height est npe=%0.3f"%npeW0)
plt.xlabel('pulse height (mV)')

plt.figure(105)
h5=np.histogram(W0peakCharges, bins=50)
values = h5[0]
centers=bincenters(h5[1])
errs=generateerrorbars(values)
plt.errorbar(centers, values, yerr=errs, fmt='o')
plt.title("W0 peak total Charge "  )
plt.xlabel('Charge (nsec*mV)')
meanChargeW0=np.mean(W0peakCharges)
stdChargeW0=np.std(W0peakCharges)
npeChargeW0= (meanChargeW0/stdChargeW0)**2
print "meanChargeW0, npeChargeW0 = ",meanChargeW0, npeChargeW0


fixedrun2V=2.8
f1.write("# Run nevts Volt mNumPksW0 mPhW0 npeW0 tW0mean tW0std mCharge npeCharge\n")
f1.write("%s %d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n" % ( 
					pathName, wantNevs, fixedrun2V, 
					meanNumPksW0, meanPhW0, 
				        npeW0, W0mean, W0std, 
					meanChargeW0, npeChargeW0  ))

f1.close()
plt.show()

sys.exit()
	
print 'linecount = ',  linecount


wavesW0rot=[]
wavesTrrot=[]
for i in range(len(waves)) :
	awaveW0 = waves[i]
	offset = TrBinOffsets[i]
	awaveW0rot = rotate( awaveW0, -offset )
	wavesW0rot.append( awaveW0rot )
	awaveTr = wavesTr[i]
	awaveTrrot = rotate( awaveTr, -offset )
	wavesTrrot.append( awaveTrrot )

from operator import add
#map(add, list1, list2)

#sumW0 = [0]*1024
#print 'length of sumW0 = ', len(sumW0)
#for i in range(len(wavesW0rot)) :
#	sumW0 = map(add, sumW0, wavesW0rot[i] )

# collect data into 3 groups of 30 triggers
sumW0s = []
aveW0s = []
for grp in range(0,3) :
   sumW0 = [0]*1024
   for i in range( 30*grp, 30*(grp+1) ) :
        sumW0 = map(add, sumW0, wavesW0rot[i] )
   sumW0s.append(sumW0)
   aveW0 = newList = [x / 30 for x in sumW0]
   aveW0s.append(aveW0)



waveNumb = 25
voffset = 10.0 # mV
voffset = 100.0 # mV
nev = 25
print "in main: length of times is ", len(times)," length of waves is ", len(waves)
voffset = 10.0 # mV
#voffset = 50.0 # mV
#myWavePlot( 100, times, waves, waveNumb, nev, voffset )
myWavePlot( 100, times, wavesW0rot, waveNumb, nev, voffset )

voffset = 100.0 # mV
#myWavePlot( 101, times, wavesTr, waveNumb, nev, voffset )
myWavePlot( 101, times, wavesTrrot, waveNumb, nev, voffset )

voffset=10.
#voffset=100.
myWavePlot( 102, times, sumW0s , 0, 3, 30*voffset  )
myWavePlot( 103, times, aveW0s , 0, 3, voffset  )

plt.show()

sys.exit()

text_file = open("binaryReader.dat", "r")
#lines = text_file.readlines()
lines = text_file.read().split(' ')
print lines[0]
print 'length of lines is ',len(lines)
print 'length of line is ',len(lines[0])
print 'lines[-1] = ', lines[-1], "\n = lines[-1]"
wave=[]
waves=[]
for i in range(len(lines)) :
	wave=[]
	aLine=lines[i]
	print 'length of aLine = ',len(aLine)
	for j in range(len(aLine)) :
		print j, aLine[j]
		wave.append( float(aLine[j] ) )
	waves.append(wave)
print wave
print 'length of wave is ', len(wave)
print 'length of waves is ', len(waves)
text_file.close()

sys.exit()

with open('special.dat') as f:
    for line in f:
        #inner_list = [elt.strip() for elt in line.split(' ')]
        # in alternative, if you need to use the file content as numbers
        inner_list = [int(elt.strip()) for elt in line.split(' ')]
        list_of_lists.append(inner_list)
print 'length of list_of_lists is ', len(list_of_lists )
