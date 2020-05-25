#ifndef MODULES_AUDIO_PROCESSING_AEC_AEC_RESAMPLER_H_
#define MODULES_AUDIO_PROCESSING_AEC_AEC_RESAMPLER_H_

#include "aec_core.h"

//namespace webrtc {

enum { kResamplingDelay = 1 };
enum { kResamplerBufferSize = FRAME_LEN * 4 };

// Unless otherwise specified, functions return 0 on success and -1 on error.
void* WebRtcAec_CreateResampler();  // Returns NULL on error.
int WebRtcAec_InitResampler(void* resampInst, int deviceSampleRateHz);
void WebRtcAec_FreeResampler(void* resampInst);
#if 0
// Estimates skew from raw measurement.
int WebRtcAec_GetSkew(void* resampInst, int rawSkew, float* skewEst);

// Resamples input using linear interpolation.
void WebRtcAec_ResampleLinear(void* resampInst,
                              const float* inspeech,
                              size_t size,
                              float skew,
                              float* outspeech,
                              size_t* size_out);

//}  // namespace webrtc
#endif
#endif  // MODULES_AUDIO_PROCESSING_AEC_AEC_RESAMPLER_H_
