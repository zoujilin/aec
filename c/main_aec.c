#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define DR_WAV_IMPLEMENTATION
#include "echo_cancellation.h"

/*****************************************/
// 输入：采样率，mic个数
// 输出：初始化指针数组
/*****************************************/
static int32_t *tmp_near[1];
static int32_t *tmp_far;
int Audio_Aec_init(uint32_t sampleRate, uint32_t channels, void *aecmInst[1])//,
{
    for(int j = 0; j < channels; j++) {
        aecmInst[j] = WebRtcAec_Create();
        if (aecmInst[j] == NULL)
            return -1;
        int status = WebRtcAec_Init(aecmInst[j], sampleRate, sampleRate);
        if (status != 0) {
            //printf("WebRtcAecm_Init fail\n");
            WebRtcAec_Free(aecmInst[j]);
            return -1;
        }
    }
	 for(int i=0;i<1;i++) {
		tmp_near[i] = (int32_t *)malloc(256*sizeof(int32_t));
	}
	tmp_far = (int32_t *)malloc(256*sizeof(int32_t));

    return 0;
}


/*****************************************/
// 输入：参考帧数据，各MIC帧数据(目前支持到4MIC)，参数指针数组，msInSndCardBuf目前固定为40，MIC个数
// 输出：处理后的MIC数据仍存放回输入MIC帧中
/*****************************************/

int Audio_Aec_Process(void *aecmInst[1], int16_t *far_frame, int16_t *near_inout[1],int16_t msInSndCardBuf, uint32_t channels)
{
    size_t samples = 256;//MIN(160, sampleRate / 100);
    
    for(int t=0; t<256; t++)
    {
        tmp_near[0][t] = near_inout[0][t];
        tmp_far[t]     = far_frame[t];
    }
    
    for (int j = 0; j < channels; j++) 
    {
        if (WebRtcAec_BufferFarend(aecmInst[j], tmp_far, samples) != 0) {//160->256 samples
            //printf("WebRtcAecm_BufferFarend() failed.");
            WebRtcAec_Free(aecmInst[j]);
            return -1;
        }

        int nRet = 0;
        nRet = WebRtcAec_Process(aecmInst[j], (const int32_t **)(&tmp_near[j]), 1,(int32_t **)(&tmp_near[j]), samples, msInSndCardBuf, 0);
        for(int t=0; t<256; t++)
        {
            near_inout[0][t] = (int16_t)tmp_near[0][t];//need to remove
        }
        if (nRet != 0) {
            //printf("failed in WebRtcAecm_Process\n");
            WebRtcAec_Free(aecmInst[j]);
            return -1;
        }
    }
    
}
void Audio_Aec_Free(void *aecmInst[1], uint32_t channels)
{
    for(int j = 0; j < channels; j++) {
        WebRtcAec_Free(aecmInst[j]);
    }
    for(int i=0;i<1;i++) {
    	free(tmp_near[i]);
	}
	free(tmp_far);
}

