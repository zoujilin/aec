#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
//#include "profiler.h"
#define DR_WAV_IMPLEMENTATION
#include "audio_aec.h"
#include "ring_buffer.h"
#include "aec_core.h"
#include "test_aec_data.h"

//#define CORE_M3

#ifdef CORE_M3
#include "iet_systick.h"
#include "BACH.h"
#include "iet_timer.h"
#include "iet_uart.h"
extern volatile unsigned long SysTickCnt_test;
#else
//#include "timing.h"//fix a compile bug at win10
#include "dr_wav.h"
#endif
#include "arm_math.h"
#include "arm_common_tables.h"
#include "arm_const_structs.h"
extern RingBuffer* buffer_;
uint32_t start;
//extern FFT_RorW_TypeDef fft_struct;
//extern uint32_t ecrd;
//extern int32_t g_fft_buf[128*2*2]  __attribute__((aligned(16)));//real,imag
extern arm_rfft_instance_q31  S_fft128;
extern arm_rfft_instance_q31  S_ifft128;
int fftin[260],fftout[260];
extern int fft_n;
extern int FilterAdaptation_n;
#ifdef CORE_M3
#include "cmsis_os2.h"
#include "iet_fft.h"
void uart_init(void)
{
    iet_systick_init();
    asr_uart_open(UART_NUM_0, UART_DATA_LEN_8, UART_STOP_BIT_1, UART_PARITY_NONE, 115200, 1);
    iet_printf_enable();
}

typedef uint32_t (*asr_wait_done)(uint32_t ui4_timeout);

static osEventFlagsId_t  os_fft_evt_id = NULL;


/*extern data*/
asr_wait_done  fft_finsh = NULL;

#define ASR_FFT_DONE_FLAGS                  (0x00000001U)
/*static function declare*/
static void asr_fft_callback(void);
static int32_t asr_fft_create_evt(void);
static uint32_t asr_fft_wait_done(uint32_t ui4_timeout);
static void asr_fft_set_callback(asr_wait_done wait_fft_done);

static uint8_t ui1_create_evt_flag = 0;

/*function implement*/
void asr_fft_int(void)
{
    uint32_t ui4_ret;

    if (ui1_create_evt_flag == 0)
    {
        ui4_ret = asr_fft_create_evt();
        ui1_create_evt_flag = 1;
    }

    ui4_ret = iet_fft_init(FFT_INT);

    ui4_ret = iet_fft_control(FFT_CMD_ISR, asr_fft_callback);

    asr_fft_set_callback(asr_fft_wait_done);
}

static int32_t asr_fft_create_evt(void)
{

    os_fft_evt_id = osEventFlagsNew(NULL);

    return 0;
}

static void asr_fft_callback(void)
{
    osEventFlagsSet(os_fft_evt_id, ASR_FFT_DONE_FLAGS);
#ifdef ASR_DEBUG_TRACK_DATA_FLOW
    ++ui4_fft_send_down_times;
#endif
}

static uint32_t asr_fft_wait_done(uint32_t ui4_timeout)
{
    uint32_t ui4_flags;
    static uint32_t ui4_fft_timeout_times = 0;

    ui4_flags = osEventFlagsWait(os_fft_evt_id, ASR_FFT_DONE_FLAGS, osFlagsWaitAny, ui4_timeout);

    if (ui4_flags == ASR_FFT_DONE_FLAGS)
    {
#ifdef ASR_DEBUG_TRACK_DATA_FLOW
        ++ui4_fft_receive_down_times;
#endif
        return 0;
    }
    else if (ui4_flags == osFlagsErrorTimeout)
    {
        ++ui4_fft_timeout_times;
        iet_printf_msg("\r\nfft error %d times !\r\n",ui4_fft_timeout_times);
        return 1;
    }
    return 0;

}

static void asr_fft_set_callback(asr_wait_done wait_fft_done)
{
    if (wait_fft_done != NULL)
    {
        fft_finsh = wait_fft_done;
    }
}
#endif


void Audio_FrontEnd_Sim(char *near_file, char *far_file, char *out_file) {

    uint32_t t0,t1;
    uint32_t sampleRate = 16000;
    uint64_t inSampleCount = 0;
    uint32_t inChannels = 1;
    uint32_t ref_sampleRate = 0;
    uint64_t ref_inSampleCount = 0;
    uint32_t ref_inChannels = 0;
    
#ifdef CORE_M3
	uart_init();

    iet_printf_msg("\r\n#### aec test start #### \r\n");

    osKernelInitialize();

    asr_fft_int();


#endif
    
    //ref_inSampleCount=32*256;
    //inSampleCount=ref_inSampleCount*inChannels;//for test
    //uint32_t i = 0, j = 0, inSamplePerChan = inSampleCount/inChannels;

#if 1
    int16_t *near = wavRead_int16(near_file, &sampleRate, &inSampleCount, &inChannels);
    int16_t *far = wavRead_int16(far_file, &ref_sampleRate, &ref_inSampleCount, &ref_inChannels);
    
    //ref_inSampleCount=32*256;//for test
    //inSampleCount=ref_inSampleCount*inChannels;//for test
    #if 0
    FILE *fp_o = fopen("neardata32.txt", "wt");
    FILE *fp_f = fopen("fardata32.txt", "wt");
    for (int i = 0; i < ref_inSampleCount; i++) {
      if ((i+1) % 32 == 0) {
        fprintf(fp_o, "\n");
        fprintf(fp_f, "\n");
      }
      fprintf(fp_o, "%d ,", near[i]);
      fprintf(fp_f, "%d ,", far[i]);
    }
    fclose(fp_o);
    fclose(fp_f);
    #endif
    uint32_t i = 0, j = 0, inSamplePerChan = inSampleCount/inChannels;
    int16_t *outBuffer = (int16_t *) malloc(inSamplePerChan*sizeof(int16_t));
    int16_t *near_frame = (int16_t *)malloc(inSampleCount*sizeof(int16_t));
    int16_t *far_frame = (int16_t *)malloc(ref_inSampleCount*sizeof(int16_t));
    for(i = 0; i < inSamplePerChan; i++) {
        //for(j = 0; j < inChannels; j++) {
            near_frame[i+0*inSamplePerChan] = near[0+i*1]; ////*4到int16最大值
        //}
    }
    for(int i = 0; i < ref_inSampleCount; i++) {
        far_frame[i] = far[i]; ///32767.0;
    }
#endif

    
    //inChannels=1;//for test
    void *aecmInst[1];
    Audio_Aec_init(sampleRate, inChannels, (void **)&aecmInst[0]);


#if 0
int fftin512[1028],fftout512[1028];
arm_rfft_instance_q31  S_fft512;
arm_rfft_instance_q31  S_ifft512;
arm_rfft_init_q31(&S_fft512, 512, 0, 1);
arm_rfft_init_q31(&S_ifft512, 512, 1, 1); 

    for(int i = 0; i < 1024; i++) {
        fftin512[i]=i+1000;
    }
    
    #ifdef CORE_M3
    //hard_rfft_int32(0, fftin512, fftout512);
    //memset(fftin512, 0, 1028*4);
    //hard_rfft_int32(1, fftout512, fftin512);
    hard_rfft_comp(0, fftin512, fftout512);
    memset(fftin512, 0, 1028*4);
    hard_rfft_comp(1, fftout512, fftin512);
    for(int i = 0; i < 514+514; i++) {
        iet_printf_msg("%d, ",fftin512[i]);
        if((i+1)%32==0)
        {
            iet_printf_msg("\r\n");
        }
    }
    iet_printf_msg("\r\n");
    #endif
    arm_rfft_q31(&S_fft512, &fftin512[512], fftout512);
    arm_rfft_q31(&S_ifft512, fftout512, fftin512);
    //arm_rfft_q31(&S_fft128, fftin512, fftout512);

#endif
    
#ifdef CORE_M3    
    for(int i = 0; i < 256; i++) {
        fftin[i]=near_frame[i];
    }
    //hard_rfft_int32(1, fftin, fftout);
    for(int i = 0; i < 256; i++) {
        //iet_printf_msg("%d, ",fftout[i]);
        if((i+1)%32==0)
        {
            //iet_printf_msg("\r\n");
        }
    }
//realFFT_forward(fftout, fftin, 128);//test fft
//arm_rfft_q31(&S_fft128, (q31_t *)fftin, (q31_t *)fftout);//test fft
//iet_printf_msg("fftout= %d, %d\r\n",fftout[0], fftout[2]);
//realFFT_Inverse(fftout, fftin, 128 );//test ifft
arm_rfft_q31(&S_ifft128, fftin, fftout);//test ifft
for(int i = 0; i < 128; i++) {
    //iet_printf_msg("%d, ",fftout[i]);
    if((i+1)%32==0)
    {
        //iet_printf_msg("\r\n");
    }
}

#endif
    int16_t echoMode = 0; // 0, 1, 2, 3 (default), 4
    int16_t msInSndCardBuf = 0;

    int16_t *near_input[1] = {near_frame};

    int16_t *far_input = far_frame;
    int16_t *sim_nearbuf[1];
    for(i=0; i<1; i++)
    {
        sim_nearbuf[i]=(int16_t *)malloc(256*sizeof(int16_t));
    }
    size_t nTotal = (inSamplePerChan / 256);
    #ifdef CORE_M3
    iet_printf_msg("nTotal=%d\r\n",nTotal);
    start = osKernelGetTickCount();
	//start = SysTickCnt_test;
	#else
    printf("nTotal = %d, Channels = %d\n", nTotal, inChannels);
    //nTotal = 584; //236 527 584
	//double startTime = now();
	#endif
	
    for (int n = 0; n < nTotal; n++) //
    {
        if(n==582)//227 137 250 235
        {
            sampleRate = 16000;
        }
        //fft_n = 0;
        //FilterAdaptation_n = 0;
        for(int t=0; t<256; t++)
        {
            sim_nearbuf[0][t] = near_input[0][t+n*256];
        }

        Audio_Aec_Process((void **)&aecmInst[0], far_input, (int16_t **)&sim_nearbuf[0], msInSndCardBuf, inChannels);

        for(int t=0; t<256; t++)//AEC result
        {
            near_input[0][t+n*256] = sim_nearbuf[0][t];
            fftout[t]              = sim_nearbuf[0][t];
            if(abs(fftout[t])>6000)
            {
                sampleRate = 16000;
            }
        }
        far_input += 256;
        //iet_printf_msg("FilterAdaptation_n = %d  \r\n", FilterAdaptation_n);
    }
    #ifdef CORE_M3
    iet_printf_msg("aec opt time =%d \r\n", osKernelGetTickCount()-start);
    iet_printf_msg("\r\n#### test end #### \r\n");
    for(int i = 0; i < ref_inSampleCount; i++) {
        //iet_printf_msg("%d, ",near_frame[i]);
        if((i+1)%32==0)
        {
            //iet_printf_msg("\r\n");
        }
    }
    #else
    //printf("aec time = %.4f ms  \r\n", (now()-startTime)*1000);
    wavWrite_int16(out_file, near_frame, sampleRate, 256*nTotal, inChannels);
    #endif
    for(i=0; i<1; i++)
    {
        free(sim_nearbuf[i]);
    }
    Audio_Aec_Free((void **)&aecmInst[0], inChannels);
    WebRtc_FreeBuffer(buffer_);
    #ifdef CORE_M3
    osKernelStart();
    #else
    #if 0
    FILE *fp3 = fopen("aec_out.txt", "wt");////////////////aec结果保存
    for(int j = 0; j < ref_inSampleCount; j++)
    {
    	fprintf(fp3,"%d, ",near_frame[j]);
    	if((j+1)%32==0)
        {
            fprintf(fp3,"\n");
        }
    }
    fclose(fp3);
    #endif
    #if 1
    free(near_frame);
    free(far_frame);
    free(outBuffer);
    free(near);
    free(far);
    #endif
    #endif

    
}

int main(int argc, char *argv[]) {
    //printf("test audio_frontEnd_process start \n");

    char *far_file  = "./wavfile/ref.wav";//"./wavfile/xmos_3_farend.wav";
    char *near_file = "./wavfile/mic.wav";//"./wavfile/xmos_3_nearend.wav";//for aec
    char *out_file = "./wavfile/out_test64_1.wav";
    Audio_FrontEnd_Sim(near_file, far_file, out_file);
    
    return 0;
}


