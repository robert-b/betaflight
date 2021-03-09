/**
 One Euro Filter, C version
 Jonathan Aceituno <join@oin.name>
 
 For details, see http://www.lifl.fr/~casiez/1euro
 */

#include <stdlib.h>
#include <math.h>

#include "maths.h"
#include "SF1eFilter.h"

void SFLowPassFilterInit(SFLowPassFilter *filter)
{
	filter->usedBefore = 0;
	filter->hatxprev = 0.0f;
	filter->xprev = 0.0f;
}

float SFLowPassFilterDo(SFLowPassFilter *filter, float x, float alpha)
{
	if (!filter->usedBefore)
	{
		filter->usedBefore = 1;
		filter->hatxprev = x;
	}

	const float hatx = alpha * x + (1.0f - alpha) * filter->hatxprev;
	filter->xprev = x;
	filter->hatxprev = hatx;
	return hatx;
}


void SF1eFilterInit(SF1eFilter *filter)
{
//	filter->rate = filter->config.rate;
//	filter->lastTime = 0;
	SFLowPassFilterInit(&(filter->xfilt));
	SFLowPassFilterInit(&(filter->dxfilt));

	filter->alpha = SF1eFilterAlpha(filter->config.rate, filter->config.derivativeCutoffFrequency);
}

float SF1eFilterDo(SF1eFilter *filter, float x)
{
	float dx = 0.0f;
	
	if (filter->xfilt.usedBefore)
	{
		dx = (x - filter->xfilt.xprev) * filter->config.rate;
	}
	
	const float edx = SFLowPassFilterDo(&(filter->dxfilt), dx, filter->alpha);
	filter->fc = filter->config.minCutoffFrequency + filter->config.cutoffSlope * fabsf(edx);
	//filter->fc = MIN(filter->fc, filter->fmax);
	return SFLowPassFilterDo(&(filter->xfilt), x, SF1eFilterAlpha(filter->config.rate, filter->fc));
}


float SF1eFilterAlpha(float rate, float cutoff)
{
	const float tau = 1.0f / (2.0f * M_PI * cutoff);
//	const float te  = 1.0f / filter->rate;
//
//	return 1.0f / (1.0f + tau / te);
	return 1.0f / (1.0f + tau * rate);
}

