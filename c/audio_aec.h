#ifndef AUDIO_AEC_H_
#define AUDIO_AEC_H_
int Audio_Aec_init(uint32_t sampleRate, uint32_t channels, void *aecmInst[1]);
int Audio_Aec_Process(void *aecmInst[1], int16_t *far_frame, int16_t *near_inout[1], int16_t msInSndCardBuf, uint32_t channels);
void Audio_Aec_Free(void *aecmInst[1], uint32_t channels);
#endif

