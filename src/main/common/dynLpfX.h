#pragma once



#define DEFAULT_DYNLPFX_ENABLE              1       //Enable DYN_LPF2 by default


#define DEFAULT_DYNLPFX_FMIN                75.0f   //Fmin in Hz
#define DEFAULT_DYNLPFX_FMAX               200.0f   //user Fmax in Hz
#define DEFAULT_DYNLPFX_DS_REV             (256)	// derivate slope, reversed

#define DEFAULT_DYNLPFX_GAIN                70      //Gain
#define DEFAULT_DYNLPFX_CENTER_THRESHOLD    10.0f   //Value in °/s

#define DEFAULT_DYNLPFX_THROTTLE_THRESHOLD  35      //Throttle in %
#define DEFAULT_DYNLPFX_THROTTLE_GAIN       12      // 12Hz / % throrrle over 35%
#define DEFAULT_DYNLPFX_Q                  (450)

#define DEFAULT_DYNLPFX_FC_FC              7.0f     //Cut of freq on FC value


#define DEFAULT_DYNLPFX_TYPE            3 //Default

#define DEFAULT_DYNLPFX_ALPHA       (750)

extern void init_dynLpfx(void);
extern float dynLpfxApply(int axis, float input);

extern void init_dynLpfxDTerm(float dT, uint16_t alpha, uint8_t filter_type);
extern float dynLpfxApplyDTerm(int axis, float input);
