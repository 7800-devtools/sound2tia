sound2tia - a WAV/OGG/FLAC to TIA conversion tool, for use with the Atari 7800
and 2600 console sound chip.

Legal Stuff
-----------
   sound2tia is created by Michael Saarna, copyright 2020.
   sound2tia is provided under the GPLv2 license. See the included LICENSE.txt
   file for details.


What Is This?
-------------
sound2tia will take as input a  WAV, OGG, or FLAC sound recording, and it will
figure out the dominant sound frequency in each 60Hz chunk of the sound, and 
provide output TIA register data in formats suitable for batari Basic, 6502 
assembly, or 7800basic.

sound2tia is a command-line utility. Run it without any arguments to see the
usage information.


Limitations
-----------
While sound2tia can create some neat conversions, the FFT analysis approach is 
doomed to mediocrity, due to basic trade-offs with FFT analysis; either you 
have good frequency accuracy over a very large period of time, or damn poor 
frequency accuracy over a relatively short period of time. 

Unfortunately 1/60th of a second sized windows have terrible frequency 
coarseness. The situation is made worse because FFTs have more coarseness of 
frequency in the lower end, and TIA has more coarseness of frequency in the 
higher end.

sound2tia tries to make the best of this situation by using a sliding window 
to cheat the trade-off, and manages some success.  It's mainly useful for 
recordings that consist of mostly pure tones, without a lot of overlapping 
tone layers, and that don't change character super quick. Whole songs aren't 
great, and intelligible human voices are a no-go.

