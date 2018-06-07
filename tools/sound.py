#!/usr/bin/env python

import numpy as np        #import needed modules
import pyaudio     #sudo apt-get install python-pyaudio
import pylab as pl
import wave as wv

PyAudio = pyaudio.PyAudio     #initialize pyaudio

#See https://en.wikipedia.org/wiki/Bit_rate#Audio
#BITRATE = 16000     #number of frames per second/frameset.
BITRATE = 40000     #number of frames per second/frameset.      

#FREQUENCY = 1.5*440     #Hz, waves per second, 261.63=C4-note.
LENGTH = 1     #seconds to play sound

#if FREQUENCY > BITRATE:
#    BITRATE = FREQUENCY*10

NUMBEROFFRAMES = int(BITRATE * LENGTH)
RESTFRAMES = NUMBEROFFRAMES % BITRATE
WAVEDATA = ''    

base = 1000 # Hz for deepest bass

dtm=pl.loadtxt(str( 100 )+'/particles.dat')
xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

y_vx = sorted( zip( ym , vxm) )

#y   = [ y   for y , vx in y_vx]

vx  = [ vx  for y , vx in y_vx]

vx_l= pl.size( vx )

nn = int( np.sqrt( vx_l ) )

vxi = np.zeros( nn )

#integrate on x

for i in range(nn) :
    vxi[i]= np.average( vx[ i : vx_l : nn ] )

vx_min = min(vxi)
vx_max = max(vxi)
vx_amp = vx_max - vx_min 

vxi = ( vxi - vx_min ) / vx_amp   # normalized between 0 and 1

vxi_l= pl.size( vxi )

#def wave(time) :
#    phase = 2 * math.pi * time * FREQUENCY
#    return 0.2 * math.sin( phase ) + 0.5 * math.sin( 2 * phase ) + 0.3 * math.sin( 3 * phase )

chunks = BITRATE / base

#generating wawes
for x in range(NUMBEROFFRAMES):
    time = x / BITRATE

    #    WAVEDATA = WAVEDATA+chr( int(  wave( time ) * 127 + 128 ))    

    # 0 to 1 in each one of the samples
    time_in = ( x % chunks ) / chunks
    wave = vxi[ int( time_in * vxi_l )  ]
    #    wave = vxi[ ( int(time * base) % int(base)   ) * vxi_l ]
    #wave = vxi[ x % vxi_l ]

    WAVEDATA = WAVEDATA+chr( int(  wave * 127 + 128 ))    

for x in range(RESTFRAMES): 
 WAVEDATA = WAVEDATA+chr(128)

p = PyAudio()

FORMAT = p.get_format_from_width(1)
stream = p.open(format = FORMAT, 
                channels = 1, 
                rate = BITRATE, 
                output = True)

stream.write(WAVEDATA)
stream.stop_stream()
stream.close()
p.terminate()

waveFile = wv.open('sound1.wav', 'wb')
waveFile.setnchannels( 1 )
waveFile.setsampwidth(p.get_sample_size(FORMAT))
waveFile.setframerate( BITRATE )
waveFile.writeframes( WAVEDATA.encode() ) # b''.join(WAVEDATA))
waveFile.close()
