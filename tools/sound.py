#!/usr/bin/env python

import math        #import needed modules
import pyaudio     #sudo apt-get install python-pyaudio

PyAudio = pyaudio.PyAudio     #initialize pyaudio

#See https://en.wikipedia.org/wiki/Bit_rate#Audio
#BITRATE = 16000     #number of frames per second/frameset.
BITRATE = 44100     #number of frames per second/frameset.      

FREQUENCY = 1.5*440     #Hz, waves per second, 261.63=C4-note.
LENGTH = 1     #seconds to play sound

if FREQUENCY > BITRATE:
    BITRATE = FREQUENCY*10

NUMBEROFFRAMES = int(BITRATE * LENGTH)
RESTFRAMES = NUMBEROFFRAMES % BITRATE
WAVEDATA = ''    

def wave(time) :
    phase = 2 * math.pi * time * FREQUENCY
    return 0.2 * math.sin( phase ) + 0.5 * math.sin( 2 * phase ) + 0.3 * math.sin( 3 * phase )

#generating wawes
for x in range(NUMBEROFFRAMES):
    time = x / BITRATE
    WAVEDATA = WAVEDATA+chr( int(  wave( time ) * 127 + 128 ))    
    
for x in range(RESTFRAMES): 
 WAVEDATA = WAVEDATA+chr(128)

p = PyAudio()
stream = p.open(format = p.get_format_from_width(1), 
                channels = 1, 
                rate = BITRATE, 
                output = True)

stream.write(WAVEDATA)
stream.stop_stream()
stream.close()
p.terminate()

