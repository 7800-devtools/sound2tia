// Provided under the GPL v2 license. See the included LICENSE.txt for details.

#define PROGNAME "Sound2TIA 0.6 RC1, by Michael Saarna, 2020."

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <unistd.h>
#include "tiavals.h"

#include "tia-fft-1.h"		// TIA frequencies in varying phase, window size=1
#include "tia-fft-2.h"		// TIA frequencies in varying phase, window size=2

int FFTDATASIZE[] = { 0, 269, 539 };

int SLICES;
long SAMPLERATE;

int WINDOWSIZE = 2;
int lowpass = -1;
int lowpassfreq;
int lowpassindex = 1;
int highpass = -1; // 1 Hz
int highpassfreq;
int highpassindex = 1;

int N, N_OUT;
long samplecount;
int algorithm = 0;

float *samplebuffer = NULL;
double *fftbuffer = NULL;

double normalizefactor;
double samplemax;

int loadwave(char *filename);
void usage(char *programname);


int main(int argc, char **argv)
{
    fftw_complex *out;
    double *in;
    double *frequency1;
    double *channel1;
    double *volume1;
    fftw_plan p;
    long t;

    int outformat = 0;
    double frequencyscale = 1.0;
    int halfrate = 0;

    int c;
    extern char *optarg;
    extern int optind, optopt;
    int errflg = 0;

    char *infilename;
    long windowindex;

    infilename = NULL;
    windowindex = 0;
    SLICES = 60;

    while ((c = getopt(argc, argv, ":i:to:m:b:s:a:l:h:")) != -1)
    {
	switch (c)
	{
	case 'i':
	    infilename = optarg;
	    break;
	case 'o':
	    outformat = atoi(optarg);
	    if ((outformat > 3) || (outformat < 0))
	    {
		fprintf(stderr, "*** ERR: option -o requires a number from 0 to 3.\n");
		errflg++;
	    }
	    break;
	case 'a':
	    algorithm = atoi(optarg);
	    if ((algorithm > 2) || (algorithm < 0))
	    {
		fprintf(stderr, "*** ERR: option -a requires a number from 0 to 2.\n");
		errflg++;
	    }
	    else if (algorithm > 0)
	    {
		WINDOWSIZE = algorithm;
	    }
            // check if a high pass filter argument was already used, and adjust if necessary
            if( (WINDOWSIZE==1) && (highpassindex>1) ) 
                highpassindex = highpassindex / 2; 
	    break;
	case 's':
	    frequencyscale = atof(optarg);
	    break;
	case ':':		/* -i or -w without operand */
	    fprintf(stderr, "*** ERR: option -%c requires an operand\n", optopt);
	    errflg++;
	    break;
	case 't':
	    halfrate = 1;
	    SLICES = 30;
	    break;
	case 'l':
	    lowpassfreq = atoi(optarg);
	    if (lowpassfreq <1 )
	    {
		fprintf(stderr, "*** ERR: lowpass value isn't valid.\n");
		errflg++;
	    }
            // calculate the index # that corresponds to that frequency
            lowpassindex = (float) lowpassfreq / (SLICES/WINDOWSIZE) ;
            if (lowpassindex == 0)
                lowpassindex = 1;
            break;
	case 'h':
	    highpassfreq = atoi(optarg);
	    if (highpassfreq <1 )
	    {
		fprintf(stderr, "*** ERR: highpass value isn't valid.\n");
		errflg++;
            }
            // calculate the index # that corresponds to that frequency
            highpassindex = (float) highpassfreq / (SLICES/WINDOWSIZE) ;
            if (highpassindex == 0)
                highpassindex = 1;
            break;
	case '?':
	    fprintf(stderr, "*** ERR: unrecognised option \"-%c\"\n", optopt);
	    errflg++;
	    break;
	}
    }
    if (errflg)
    {
	usage(argv[0]);
	return (2);
    }

    if (infilename == NULL)
    {
	fprintf(stderr, "-i wav file parameter required\n\n");
	usage(argv[0]);
	return 1;
    }

    if (loadwave(infilename) != 0)
    {
	fprintf(stderr, "loading wav file failed\n");
	return 1;

    }

    N = (SAMPLERATE * WINDOWSIZE / SLICES);
    N_OUT = ((N << 1) + 1);

    in = (double *) fftw_malloc(sizeof(double) * 2 * N);	// too much
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N_OUT);

    frequency1 = (double *) malloc(sizeof(double) * (samplecount * WINDOWSIZE / N));
    volume1 = (double *) malloc(sizeof(double) * (samplecount * WINDOWSIZE / N));
    channel1 = (double *) malloc(sizeof(double) * (samplecount * WINDOWSIZE / N));

    memset(frequency1, (samplecount * WINDOWSIZE / N), sizeof(double));
    memset(volume1, (samplecount * WINDOWSIZE / N), sizeof(double));
    memset(channel1, (samplecount * WINDOWSIZE / N), sizeof(double));

    if ((in == NULL) || (out == NULL) || (frequency1 == NULL) || (volume1 == NULL))
    {
	fprintf(stderr, "** error allocating double array\n");
	return (1);
    }

    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    if (p == NULL)
    {
	fprintf(stderr, "** error creating plan\n");
	return (1);
    }
    if ((windowindex * N) > (samplecount - N))
    {
	fprintf(stderr, "** Window index %ld is out of sample range\n", windowindex);
	return (1);
    }

    //normalize the entire sample to range from -100 to 100
    samplemax = 0;
    for (t = 0; t < samplecount; t++)
    {
	if (samplebuffer[t] > samplemax)
	    samplemax = samplebuffer[t];
	if (samplebuffer[t] < 0)
	    if ((0.0 - samplebuffer[t]) > samplemax)
		samplemax = (0.0 - samplebuffer[t]);
    }
    normalizefactor = (100.0 / samplemax);

    for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
    {
	for (t = 0; t < N; t++)
	{
	    in[t] = samplebuffer[t + (windowindex * N / WINDOWSIZE)] * normalizefactor;
	}
	fftw_execute(p);	// run the FFT on our window

	// Figure out the maximum 2 frequency values for this window
	float maximumfftval1;
	float maximumfftfreq1;

	maximumfftval1 = 0;
	maximumfftfreq1 = 0;

	float currentfftval;
	for (t = highpassindex; t < N / 2; t++)
	{
	    // we use 1.25 instead of 1 to favor fundamentals. 
	    // Ideally this should be made tunable...
	    out[t][1] = fabs(out[t][1]) / (t * 1.25);

	    currentfftval = 2 * out[t][1];
	    if (currentfftval > maximumfftval1)
	    {
		maximumfftval1 = currentfftval;
		maximumfftfreq1 = t * SLICES / WINDOWSIZE;
	    }
	}

	if (algorithm == 0)
	{
	    frequency1[windowindex] = maximumfftfreq1;
	    volume1[windowindex] = maximumfftval1;
	}
	else			// algorithm > 0
	{

	    if (maximumfftval1 == 0.0)	// if no frequency content was found, handle it here so we don't divide by zero
	    {
		channel1[windowindex] = 0x4;
		frequency1[windowindex] = 0;
		volume1[windowindex] = 0;
	    }
	    else
	    {
		int peakvolume = 0;
		for (t = highpassindex; t < N; t++)
		{
		    if (peakvolume < (int) fabs(in[t]) * 15.0 / 100.0)
			peakvolume = (int) fabs(in[t]) * 15.0 / 100.0;
		}

		long bestscore;
		int bestscoreidx;
		long currentscore;
		int tiaindex;
		int maxhigh;

		maxhigh = FFTDATASIZE[algorithm];	//change size according to fft data-set
		bestscore = 2000000;
		bestscoreidx = 0;

		if (N_OUT < maxhigh)
		    maxhigh = N_OUT;

/*
			// TODO: lowpass isn't working right...
                        long lowpassindex;
                        lowpassindex=(4000*WINDOWSIZE)/SLICES; // 4000hz cutoff
                        if((lowpass==1)&&(maxhigh>lowpassindex))
                                maxhigh=lowpassindex;
*/

		if (algorithm == 1)
		{
		    for (tiaindex = 0; fftlists_1[tiaindex] != NULL; tiaindex++)
		    {
			currentscore = 0;
			for (t = highpassindex; t <= maxhigh; t++)
			    currentscore +=
				fabs((float) fftlists_1[tiaindex][t] - ((out[t][1] * 200.0) / maximumfftval1));
			if (currentscore < bestscore)
			{
			    bestscore = currentscore;
			    bestscoreidx = tiaindex;
			}
		    }
		    assert(bestscore != 2000000);	// even the worst match will score better than this.
		    channel1[windowindex] = tiachannels_1[bestscoreidx];
		    frequency1[windowindex] = tiafreqs_1[bestscoreidx];
		    volume1[windowindex] = peakvolume;
		}
		else
		{
                    // TODO: join the data structures for fft1 and fft2 better, so this code block
                    // doesn't need to essentially be a duplicate of the above one.
		    for (tiaindex = 0; fftlists_2[tiaindex] != NULL; tiaindex++)
		    {
			currentscore = 0;
			for (t = highpassindex; t <= maxhigh; t++)
			    currentscore +=
				fabs((float) fftlists_2[tiaindex][t] - ((out[t][1] * 200.0) / maximumfftval1));
			if (currentscore < bestscore)
			{
			    bestscore = currentscore;
			    bestscoreidx = tiaindex;
			}
		    }
		    assert(bestscore != 2000000);	// even the worst match will score better than this.
		    channel1[windowindex] = tiachannels_2[bestscoreidx];
		    frequency1[windowindex] = tiafreqs_2[bestscoreidx];
		    volume1[windowindex] = peakvolume;

		}
	    }
	}
    }




    float averagevolume = 0;
    // first we calculate the sample average volume...
    samplemax = 0.0;
    for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	averagevolume = averagevolume + volume1[windowindex];
    averagevolume = averagevolume / windowindex;

    // then normalize the volume of the sample range, so the the average volume is now 10 out of 15.
    normalizefactor = (10.00 / averagevolume);

    // we boosted to 10/15, so some values likely exceed TIA's max volume of 15. Let's clip them...
    for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
    {
	volume1[windowindex] = volume1[windowindex] * normalizefactor;
	if (volume1[windowindex] > 15)
	    volume1[windowindex] = 15;
    }

    if (algorithm == 0)
    {
	if (frequencyscale != 1.0)
	{
	    // frequency scaling is specified
	    for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
		frequency1[windowindex] = frequency1[windowindex] * frequencyscale;
	}


	int frequencydelta1, diff1;
	int chan1, freq1;
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    frequencydelta1 = 30000;
	    for (t = 0; TIAVALS[t][0] != 0; t++)
	    {
		diff1 = abs((int) frequency1[windowindex] - (int) TIAVALS[t][0]);
		if (diff1 < frequencydelta1)
		{
		    frequencydelta1 = diff1;
		    chan1 = TIAVALS[t][1];
		    freq1 = TIAVALS[t][2];
		}
	    }
	    channel1[windowindex] = chan1;
	    frequency1[windowindex] = freq1;
	}

    }
    if (outformat == 0)
    {
	printf("\nAUDF0\tAUDC0\tAUDV0\n-----\t-----\t-----\n");
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x\t$%02x\t$%02x\n",
		   (int) frequency1[windowindex], (int) channel1[windowindex], (int) volume1[windowindex]);
	}
    }

    if (outformat == 1)
    {

	printf("\n data AUDCVDATA0\n");
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x", (int) channel1[windowindex] * 16 + (int) volume1[windowindex]);
	    if (windowindex % 8 == 7)
		printf("\n");
	    else if (windowindex == (samplecount / N - 1))
		printf("\n 0\nend\n");
	    else
		printf(", ");
	}
	printf("\n data AUDFDATA0\n");
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x", (int) frequency1[windowindex]);
	    if (windowindex % 8 == 7)
		printf("\n");
	    else if (windowindex == (samplecount / N - 1))
		printf("\nend\n");
	    else
		printf(", ");
	}
    }
    if (outformat == 2)
    {

	printf("\nAUDCVDATA0\n .byte");
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x", (int) channel1[windowindex] * 16 + (int) volume1[windowindex]);
	    if (windowindex % 8 == 7)
		printf("\n .byte");
	    else if (windowindex == (samplecount / N - 1))
		printf("\n");
	    else
		printf(", ");
	}
	printf("\nAUDFDATA0\n .byte");
	for (windowindex = 0; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x", (int) frequency1[windowindex]);
	    if (windowindex % 8 == 7)
		printf("\n .byte");
	    else if (windowindex == (samplecount / N - 1))
		printf("\n");
	    else
		printf(", ");
	}
    }
    if (outformat == 3)
    {

	printf("\n data sfx_samplesound\n");
	if (halfrate == 0)
	    printf(" $10,$10,$00 ; version, priority, frames per chunk\n");
	else
	    printf(" $10,$10,$01 ; version, priority, frames per chunk\n");
	printf(" $%02x,$%02x,$%02x ; first chunk of freq,channel,volume\n",
	       (int) frequency1[0], (int) channel1[0], (int) volume1[0]);
	for (windowindex = 1; windowindex < (samplecount * WINDOWSIZE / N); windowindex = windowindex + 1)
	{
	    printf(" $%02x,$%02x,$%02x\n",
		   (int) frequency1[windowindex], (int) channel1[windowindex], (int) volume1[windowindex]);
	}
	printf(" $00,$00,$00\nend\n");
    }


    // we're done. time to clean-up.
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    free(frequency1);
    free(channel1);
    free(volume1);
    free(samplebuffer);
    return (0);
}


int loadwave(char *filename)
{
    SF_INFO sndInfo;
    SNDFILE *sndFile;
    long numFrames;
    int channels;

    // Open sound file
    sndFile = sf_open(filename, SFM_READ, &sndInfo);
    if (sndFile == NULL)
    {
	fprintf(stderr, "Error reading wav file '%s': %s\n", filename, sf_strerror(sndFile));
	return 1;
    }

    channels = sndInfo.channels;

    SAMPLERATE = sndInfo.samplerate;

    N = (SAMPLERATE / SLICES);

    // Allocate memory
    samplecount = sndInfo.frames;
    if ((samplecount % N) != 0)
	samplecount = samplecount + N - (samplecount % N);

    samplebuffer = calloc(samplecount * 2, sizeof(float));	// too much
    if (samplebuffer == NULL)
    {
	fprintf(stderr, "Could not allocate memory for file\n");
	sf_close(sndFile);
	return 1;
    }

    // Load the sample data into a array of floats...
    if (channels == 1)
	numFrames = sf_readf_float(sndFile, samplebuffer, sndInfo.frames);
    else			// channels>1
    {
	// we have at least two channels, but we only want one.
	// drop all channels except the first (left) one...
	long s = 0, t = 0;
	float *samplebuffertmp;
	samplebuffertmp = calloc(samplecount * 2, sizeof(float) * channels);	// too much
	if (samplebuffertmp == NULL)
	{
	    if (samplebuffer != NULL)
		free(samplebuffer);
	    fprintf(stderr, "Could not allocate memory for file\n");
	    sf_close(sndFile);
	    return 1;
	    // there should probably be some more allocation freeing in here, but
	    // who cleans up after themselves these days?
	}

	numFrames = sf_readf_float(sndFile, samplebuffertmp, sndInfo.frames);
	for (t = 0; t < numFrames; t++)
	{
	    samplebuffer[t] = samplebuffertmp[s];
	    s = s + channels;
	}
	free(samplebuffertmp);
    }

    // Check correct number of samples loaded
    if (numFrames != sndInfo.frames)
    {
	fprintf(stderr, "Did not read enough frames for source\n");
	sf_close(sndFile);
	free(samplebuffer);
	return 1;
    }

    sf_close(sndFile);
    return 0;
}

void usage(char *programname)
{
    fprintf(stderr, "%s %s %s\n", PROGNAME, __DATE__, __TIME__);
    fprintf(stderr, "Usage: %s -i INPUTFILE [-o OUTFORMAT] [-a #] [-s #] [-t] [-l #] [-h #]\n", programname);
    fprintf(stderr, "       where INPUTFILE is a mono or stereo WAV, OGG, or FLAC file.\n");
    fprintf(stderr, "             OUTPUTFORMAT is 0 for raw, 1 for bB, 2 for asm, 3 for 7800basic.\n");
    fprintf(stderr, "             -a # is freqency match algorithm. 0=peak (default) 1=fft1 2=fft2.\n");
    fprintf(stderr, "             -s # is frequency scale. e.g. -f 0.5 lowers frequencies in half.\n");
    fprintf(stderr, "             -t provides TIA data intended to play at 30Hz.\n");
    fprintf(stderr, "             -l # is a low pass filter, with value in hz. e.g. -l 4000\n");
    fprintf(stderr, "             -h # is a high pass filter, with value in hz. e.g. -l 120\n");
    fprintf(stderr, "\n");
}
