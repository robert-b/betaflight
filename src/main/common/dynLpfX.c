
#include <stdbool.h>

#include "platform.h"
#include "math.h"
#include "maths.h"

#include "dynLpfx.h"
#include "SF1eFilter.h"

#include "fc/rc.h"

#include "build/debug.h"

#include "common/filter.h"

#include "sensors/gyro.h"

#include "fc/rc_controls.h"

//DEFINITIONS
//-----------
#define DYNLPF2_HYTEREIS  (2.0f)  //Value in °/s
#define TWO_PI 			  (2.0f * M_PIf)
#define BASE_LPF_HZ    	  (70.0f)


SF1eFilterConfiguration config;
SF1eFilter dynlpfX[3];

alphaBetaGammaFilter_t abgFilter[3];
alphaBetaGammaFilter_t abgFilterDTerm[3];

//TYPES
//-----


typedef struct
{
	pt1Filter_t pt1;          //PT1 filter

	float Fc;                 //Cutoff freq
	pt1Filter_t pt1Fc;    //PT1 on Fc

	bool Dyn_Fc;         //Dynamic E or Fixed E
} dynlpf2_t;


#if (defined(STM32F7) || defined(STM32F4))
#define MAX_WINDOW_SIZE 64
#else
#define MAX_WINDOW_SIZE 32
#endif

typedef struct kalman_s {
    uint32_t w;    // window size
    float q;       // process noise covariance
    float r;       // measurement noise covariance
    float p;       // estimation error covariance matrix
    float x;       // state
    float lastX;   // previous state

    float window[MAX_WINDOW_SIZE];
    float variance;
    float varianceSum;
    float mean;
    float meanSum;
    float windowSizeInverse;
    uint32_t windowIndex;
//    pt1Filter_t lp_filter;
//    float oldSetPoint;
    float updateRate;
} FastKalmanFilter_t;

//VARIABLES
//---------
    static dynlpf2_t dynLpf[3];
    static FastKalmanFilter_t fkFilter[3];

    static float Fmin_init,Fmax;
    static float throttleThreshold;
    static float throttleGain;
    static float dynFcThreshold;
    static float dynGainOnError;
    static flight_dynamics_index_t gyroDebugAxis;
    static float gyroDt;


// 0.3f caused a real low cutoff frequency of around 20 hz
//const float r_weight = 0.67f;

// Proper fast two-state Kalman
void fastKalmanInit(FastKalmanFilter_t *filter, float q, uint32_t w, float dT)
{
    if ( w > MAX_WINDOW_SIZE)
    {
    	w = MAX_WINDOW_SIZE;
    }

    memset(filter, 0, sizeof(FastKalmanFilter_t));
    // original: filter->q     = q * 0.000001f; // add multiplier to make tuning easier - not enough on the dynamics.
    filter->q     = q * 0.0001f; // add multiplier to make tuning easier
    filter->p     = q * 0.001f;    // add multiplier to make tuning easier
    filter->r     = q * 0.001f;    // add multiplier to make tuning easier
    filter->w     = w;
    filter->windowSizeInverse = 1.0f/(w - 1);
   	filter->updateRate = 1.0f / dT;
}

FAST_CODE float fastKalmanUpdate(FastKalmanFilter_t *filter, float input, float setPoint)
{
	static float e, Error;

	const float filteredValue = filter->x;

	float Average = fabsf(setPoint + filteredValue) * 0.5f;		// less input noise when using the state
	if (Average == 0.0f)
	{
		Average = 0.00001f;
	}

    // project the state ahead using acceleration
    filter->x += (filter->x - filter->lastX);

    // update last state
    filter->lastX = filter->x;

	// figure out how much to boost or reduce our error in the estimate based on setPoint target.
	// this should be close to 0 as we approach the setPoint and really high the further away we are from the setPoint.


	Error = fabsf(setPoint - input);
	e = 1.0f + Error / Average; 	// e is the error noise in percent of the signal

    // prediction update
    filter->p = filter->p + filter->q * e;

    // measurement update
    const float k = filter->p / (filter->p + filter->r);
    filter->x += k * (input - filter->x);
    filter->p = (1.0f - k) * filter->p;

    // update variance
    filter->window[filter->windowIndex] = input;

    filter->meanSum += filter->window[filter->windowIndex];
    filter->varianceSum = filter->varianceSum + (filter->window[filter->windowIndex] * filter->window[filter->windowIndex]);

    filter->windowIndex++;
    if (filter->windowIndex >= filter->w)
    {
        filter->windowIndex = 0;
    }

    filter->meanSum -= filter->window[filter->windowIndex];
    filter->varianceSum = filter->varianceSum - (filter->window[filter->windowIndex] * filter->window[filter->windowIndex]);

    filter->mean = filter->meanSum * filter->windowSizeInverse;
    filter->variance = fabsf(filter->varianceSum * filter->windowSizeInverse - (filter->mean * filter->mean));
    filter->r = sqrtf(filter->variance);

    return filter->x;
}


//////////////////////////////
//                          //
//       DYN PT1 INIT       //
//                          //
//////////////////////////////
void init_dynLpfxDTerm(float dT, uint16_t alpha, uint8_t filter_type)
{
    // ABG filter
	ABGInit(&abgFilterDTerm[0], alpha, dT, filter_type);
	ABGInit(&abgFilterDTerm[1], alpha, dT, filter_type);
	ABGInit(&abgFilterDTerm[2], alpha, dT, filter_type);
}

//////////////////////////////
//                          //
//       DYN PT1 INIT       //
//                          //
//////////////////////////////
void init_dynLpfx(void)
{
    gyroDt = gyro.targetLooptime * 1e-6f;

    // ABG filter
	ABGInit(&abgFilter[0], gyroConfig()->dynlpfx_alpha, gyroDt, gyroConfig()->dynlpfx_abg_filter_type);
	ABGInit(&abgFilter[1], gyroConfig()->dynlpfx_alpha, gyroDt, gyroConfig()->dynlpfx_abg_filter_type);
	ABGInit(&abgFilter[2], gyroConfig()->dynlpfx_alpha, gyroDt, gyroConfig()->dynlpfx_abg_filter_type);

    //Init PT1
    float gain = pt1FilterGain(gyroConfig()->dynlpfx_fmin, gyroDt);
    pt1FilterInit(&dynLpf[0].pt1, gain);
    pt1FilterInit(&dynLpf[1].pt1, gain);
    pt1FilterInit(&dynLpf[2].pt1, gain);

    //Fc filter
    gain = pt1FilterGain(gyroConfig()->dynlpfx_fc_fc, gyroDt);
    pt1FilterInit(&dynLpf[0].pt1Fc, gain);
    pt1FilterInit(&dynLpf[1].pt1Fc, gain);
    pt1FilterInit(&dynLpf[2].pt1Fc, gain);

    dynLpf[0].Fc = gyroConfig()->dynlpfx_fmin;
    dynLpf[1].Fc = gyroConfig()->dynlpfx_fmin;
    dynLpf[2].Fc = gyroConfig()->dynlpfx_fmin;

    dynLpf[0].Dyn_Fc = false;
    dynLpf[1].Dyn_Fc = false;
    dynLpf[2].Dyn_Fc = false;

    Fmax                 = (float)gyroConfig()->dynlpfx_fmax;         //PT1 maxFc in Hz
    Fmin_init            = (float)gyroConfig()->dynlpfx_fmin;         //PT1 min Fc in Hz
   
    throttleThreshold    = (float)gyroConfig()->dynlpfx_throttle_threshold;
    throttleGain         = (float)gyroConfig()->dynlpfx_throttle_gain;

    dynFcThreshold       = (float)(gyroConfig()->dynlpfx_center_threshold);   //Min Setpoint & Gyro value to rise PT1 Fc
    dynGainOnError       = (float)(gyroConfig()->dynlpfx_gain);
    gyroDebugAxis        = gyroConfig()->gyro_filter_debug_axis;    

	config.rate = 1.0f / gyroDt;
	config.cutoffSlope = 1.0f / (float)gyroConfig()->dynlpfx_cutoffSlope;
	config.derivativeCutoffFrequency = gyroConfig()->dynlpfx_fc_fc;
	config.minCutoffFrequency = Fmin_init;

	dynlpfX[0].fmax = Fmax;
	dynlpfX[1].fmax = Fmax;
	dynlpfX[2].fmax = Fmax;

	dynlpfX[0].config = config;
	dynlpfX[1].config = config;
	dynlpfX[2].config = config;

	SF1eFilterInit(&dynlpfX[0]);
	SF1eFilterInit(&dynlpfX[1]);
	SF1eFilterInit(&dynlpfX[2]);

    // iir filter
	fastKalmanInit(&fkFilter[0], (float)gyroConfig()->dynlpfx_Q, 32, gyroDt);
	fastKalmanInit(&fkFilter[1], (float)gyroConfig()->dynlpfx_Q, 32, gyroDt);
	fastKalmanInit(&fkFilter[2], (float)gyroConfig()->dynlpfx_Q, 32, gyroDt);
}

//////////////////////////////
//                          //
//      DYN LPF PROCESS     //
//         on ratio         //
//                          //
//////////////////////////////
FAST_CODE float dynlpf2_process_type1(dynlpf2_t* filter, float input, float target) {

float newFc, Fmin;
float throttle;
float Average;

Fmin = Fmin_init;
throttle  = (rcCommand[THROTTLE] - 1000.0f) * 0.1f; //Throttle scaled to [0-100]
const float gyroDt = gyro.targetLooptime * 1e-6f;

    //Compute average between setpoint and Gyro
    //-----------------------------------------
        Average = (target + input) * 0.5f;

    //Check if setpoint or gyro are high enought to compute "e" ratio
    //---------------------------------------------------------------
        if(filter->Dyn_Fc) {
            if ( (float)(fabs(Average)) < (dynFcThreshold - DYNLPF2_HYTEREIS) ) {
                filter->Dyn_Fc = false;
            }
        }else{
            //Enable Dyn_Fc when stick or Quad move
            if ( (float)(fabs(Average)) > (dynFcThreshold + DYNLPF2_HYTEREIS) ) {
                filter->Dyn_Fc = true;
            }
        }

    //Rise Fmin according to Throttle;
    //--------------------------------
        if(throttle > throttleThreshold){
            Fmin += (throttle - throttleThreshold) * throttleGain;
        }


    //Compute e & Fc
    //--------------
        if(filter->Dyn_Fc) {

            //Avoid division by 0.0f
                if(target == 0.0f)              { target = 0.00001f; }
                if(filter->pt1.state == 0.0f)   { filter->pt1.state = 0.0001f; }

            //Compute e factor
                float Error, e;
                Error =  (float)fabs( target - input );
                e =  Error / Average;                           //Compute ratio between Error and average. e is image of noise in % of signal

            //New freq
                newFc = Fmin + dynGainOnError * 100.0f  * powf(e, 3.0f);  //"e" power 3 and multiply by a gain

        } else {
                newFc  = Fmin;
        }

    //Limit & Filter newFc
    //---------------------
        //Low Limit
        if(newFc < Fmin)  { newFc  = Fmin; }

        //Interne high Limit
        if(newFc > 1000.0f)  { newFc  = 1000.0f; }

        //Filter the cut-off freq ;)
        newFc = pt1FilterApply(&filter->pt1Fc, newFc);

        //User high Limit
        if(newFc > Fmax)  { newFc  = Fmax; }

    //Update PT1 filter
    //------------------
        pt1FilterUpdateCutoff(&filter->pt1, pt1FilterGain(newFc, gyroDt));
        filter->Fc = newFc;

    //Apply filter
    //------------
        float output = pt1FilterApply(&filter->pt1, input);

 return output;
}


//////////////////////////////
//                          //
//      DYN LPF PROCESS     //
//         on error         //
//                          //
//////////////////////////////
FAST_CODE float dynlpf2_process_type2(dynlpf2_t* filter, float input, float target) {

float newFc, Fmin;
float throttle;

Fmin = Fmin_init;
throttle  = (rcCommand[THROTTLE] - 1000.0f) * 0.1f; //Throttle scaled to [0-100]
const float gyroDt = gyro.targetLooptime * 1e-6f;

    //Rise Fmin according to Throttle;
    //--------------------------------
        if(throttle > throttleThreshold){
            Fmin += (throttle - throttleThreshold) * throttleGain;
        }

    //Compute new Fc according to error.
    //----------------------------------
        float Error, Min, scaleFactor;

        //Compute scale factor : scaleFactor = 1 near to 0°/s and 0,25 next to 1000°/s
            Min = (float)MIN(input, target);
            scaleFactor = 1.0f / (1.0f + 3.0f * (Min / 1000.0f));

        //Compute Error
            Error =  (float)fabs( target - input ) * scaleFactor;

        //newFc
            float e = Error * dynGainOnError * 0.01f;
            newFc = Fmin + powf(e, 4.0f);               //"e" power 4 to amplify big error.


    //Limit & Filter newFc
    //---------------------
        //Low Limit
        if(newFc < Fmin)  { newFc  = Fmin; }

        //Interne high Limit
        if(newFc > 1000.0f)  { newFc  = 1000.0f; }

        //Filter the cut-off freq ;)
        newFc = pt1FilterApply(&filter->pt1Fc, newFc);

        //User high Limit
        if(newFc > Fmax)  { newFc  = Fmax; }

    //Update PT1 filter
    //------------------
        pt1FilterUpdateCutoff(&filter->pt1, pt1FilterGain(newFc, gyroDt));
        filter->Fc = newFc;

    //Apply filter
    //------------
        float output = pt1FilterApply(&filter->pt1, input);

 return output;
}

float dynLpfxApplyDTerm(int axis, float input)
{
	return  alphaBetaGammaApply(&abgFilterDTerm[axis], input);
}

float dynLpfxApply(int axis, float input)
{

	float output;
	const float target = getSetpointRate(axis);

  //Apply filter if filter is enable.
    if (gyroConfig()->dynlpfx_enable != 0)
    {
    	switch (gyroConfig()->dynlpfx_type)
    	{
			default:
    		case 0:
    			// abg prefiltering filter
    			output = dynlpf2_process_type1(&dynLpf[axis], input, target);
    			break;
    		case 1:
    			// abg prefiltering filter
    			output = dynlpf2_process_type2(&dynLpf[axis], input, target);
    			break;
    		case 2:
    			// abg filter only
    			output =  alphaBetaGammaApply(&abgFilter[axis], input);
    			break;
    		case 3:
    			// abg filter only
    			output = fastKalmanUpdate(&fkFilter[axis], input, target);
    			break;
    		case 4:
    			//input =  alphaBetaGammaApply(&abgFilter[axis], input);
    			output = SF1eFilterDo(&dynlpfX[axis], input);
    			break;
    	}
    }
    else
    {
      output = input;
      dynLpf[axis].Fc = 0;  //To show filter is disable in blackbox.
    }


  //Blackbox
    if (axis == gyroDebugAxis)
    {
    	switch (gyroConfig()->dynlpfx_type)
    	{
			default:
    		case 0:
    		case 1:
    	        DEBUG_SET(DEBUG_DYN_LPFX, 0, (int16_t)(lrintf(input)));
    	        DEBUG_SET(DEBUG_DYN_LPFX, 1, (int16_t)(lrintf(output)));
    	        DEBUG_SET(DEBUG_DYN_LPFX, 2, (int16_t)(lrintf(dynLpf[axis].Fc)));
    	        DEBUG_SET(DEBUG_DYN_LPFX, 3, (int16_t)(lrintf(target)));
    			break;
    		case 2:
    			DEBUG_SET(DEBUG_DYN_LPFX, 0, dynlpfX[axis].fc);
                DEBUG_SET(DEBUG_DYN_LPFX, 1, (int16_t)(lrintf(output)));
    			break;
    		case 3:
                DEBUG_SET(DEBUG_DYN_LPFX, 0, (int16_t)(lrintf(fkFilter[axis].r * 1000)));
    			break;
    		case 4:
    			DEBUG_SET(DEBUG_DYN_LPFX, 0, dynlpfX[axis].fc);
    			break;
    	}
    }

    return output;
}
