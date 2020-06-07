/*
 * The core AEC algorithm, which is presented with time-aligned signals.
 */

#include "aec_core.h"
#include "echo_cancellation.h"
#include <math.h>
#include <stddef.h>  // size_t
#include <stdlib.h>
#include <string.h>
//#include <algorithm>
#include <assert.h>

//#include "profiler.h"

#include "fft_c.h"
//extern "C" {
#include "ring_buffer.h"
//} 
#include "signal_processing_library.h"
#include "aec_common.h"
//#include "aec_core_optimized_methods.h"
#include "delay_estimator_wrapper.h"
#include "arch.h"
uint32_t t0,t1;
static const uint32_t kMaxSeedUsed = 0x80000000;
const static int32_t aecq0 = 0;//-7
static int32_t aecq = -15;
//static int32_t aecfarq = 0;
static int32_t aecnearq = -31;
//static int32_t aecsq = aecq0;
//static int32_t init_flag_q = 0;

#ifdef CORE_M3
#include "iet_systick.h"
#include "BACH.h"
#include "iet_timer.h"
#include "iet_uart.h"
extern volatile unsigned long SysTickCnt_test;
uint32_t start0;
#endif
extern  AecCore AecCorestruct;
extern  Aec Aecstruct;
#include "arm_math.h"
#include "arm_common_tables.h"
#include "arm_const_structs.h"

uint32_t start;

//#define HARD_FFT
//#define INT_MODE

#ifdef HARD_FFT
#include "iet_fft.h"
FFT_RorW_TypeDef fft_struct;
uint32_t ecrd;
int32_t g_fft_buf[512*2*2]  __attribute__((aligned(16)))={0};//real,imag
//#else

#endif
arm_rfft_instance_q31  S_fft128;
arm_rfft_instance_q31  S_ifft128;
static uint32_t IncreaseSeed(uint32_t* seed) {
  seed[0] = (seed[0] * ((int32_t)69069) + 1) & (kMaxSeedUsed - 1);
  return seed[0];
}

int16_t WebRtcSpl_RandU(uint32_t* seed) {
  return (int16_t)(IncreaseSeed(seed) >> 16);
}

// Creates an array of uniformly distributed variables.
int16_t WebRtcSpl_RandUArray(int16_t* vector,
                             int16_t vector_length,
                             uint32_t* seed) {
  int i;
  for (i = 0; i < vector_length; i++) {
    vector[i] = WebRtcSpl_RandU(seed);
  }
  return vector_length;
}

//namespace webrtc {
//namespace {
enum  DelaySource {
  kSystemDelay,    // The delay values come from the OS.
  kDelayAgnostic,  // The delay values come from the DA-AEC.
};

const int kMinDelayLogValue = -200;
const int kMaxDelayLogValue = 200;
const int kNumDelayLogBuckets = 100;

//}  // namespace

// Buffer size (samples)
//static const size_t kBufferSizeBlocks = 250;  // 1 second of audio in 16 kHz.
#define kBufferSizeBlocks 250
// Metrics
static const size_t kSubCountLen = 4;
static const size_t kCountLen = 50;
static const int kDelayMetricsAggregationWindow = 1250;  // 5 seconds at 16 kHz.

// Divergence metric is based on audio level, which gets updated every
// |kSubCountLen + 1| * PART_LEN samples. Divergence metric takes the statistics
// of |kDivergentFilterFractionAggregationWindowSize| audio levels. The
// following value corresponds to 1 second at 16 kHz.
static const int kDivergentFilterFractionAggregationWindowSize = 50;

// Quantities to control H band scaling for SWB input
//static const float cnScaleHband = 0.4f;  // scale for comfort noise in H band.
// Initial bin for averaging nlp gain in low band
static const int freqAvgIc = PART_LEN / 2;

// Matlab code to produce table:
// win = sqrt(hanning(127)); win = [0 ; win(1:64)];
// fprintf(1, '\t%.14f, %.14f, %.14f,\n', win);
const int16_t  WebRtcAec_sqrtHanning[128]  __attribute__((aligned(16))) = {//Q15
        0,   804,  1608,  2411,  3212,  4011,  4808,  5602,
	 6393,  7180,  7962,  8740,  9512, 10279, 11039, 11793,
	12540, 13279, 14010, 14733, 15447, 16151, 16846, 17531,
	18205, 18868, 19520, 20160, 20788, 21403, 22006, 22595,
	23170, 23732, 24279, 24812, 25330, 25833, 26320, 26791,
	27246, 27684, 28106, 28511, 28899, 29269, 29622, 29957,
	30274, 30572, 30853, 31114, 31357, 31581, 31786, 31972,
	32138, 32286, 32413, 32522, 32610, 32679, 32729, 32758,
	32767, 32758, 32729, 32679, 32610, 32522, 32413, 32286,//
	32138, 31972, 31786, 31581, 31357, 31114, 30853, 30572,
	30274, 29957, 29622, 29269, 28899, 28511, 28106, 27684,
	27246, 26791, 26320, 25833, 25330, 24812, 24279, 23732,
	23170, 22595, 22006, 21403, 20788, 20160, 19520, 18868,
	18205, 17531, 16846, 16151, 15447, 14733, 14010, 13279,
	12540, 11793, 11039, 10279,  9512,  8740,  7962,  7180,
	 6393,  5602,  4808,  4011,  3212,  2411,  1608,   804,
};

// Matlab code to produce table:
// weightCurve = [0 ; 0.3 * sqrt(linspace(0,1,64))' + 0.1];
// fprintf(1, '\t%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,\n', weightCurve);

ALIGN16_BEG const int ALIGN16_END WebRtcAec_weightCurve[65] = {//Q16
        0, 6554, 9031, 10057, 10844, 11508, 12092, 12621, 
	13107, 13560, 13985, 14387, 14769, 15134, 15485, 15822, 
	16147, 16462, 16767, 17063, 17351, 17631, 17905, 18172, 
	18433, 18689, 18939, 19184, 19425, 19661, 19893, 20121, 
	20345, 20566, 20783, 20997, 21208, 21416, 21621, 21823, 
	22023, 22220, 22414, 22607, 22797, 22984, 23170, 23354, 
	23535, 23715, 23893, 24069, 24243, 24416, 24587, 24756, 
	24924, 25090, 25255, 25418, 25580, 25741, 25900, 26058, 
	26214,
};
// Matlab code to produce table:
// overDriveCurve = [sqrt(linspace(0,1,65))' + 1];
// fprintf(1, '\t%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,\n', overDriveCurve);
ALIGN16_BEG const int ALIGN16_END WebRtcAec_overDriveCurve[65] = {//Q9
    512, 576, 603, 623, 640, 655, 669, 681, 
	693, 704, 714, 724, 734, 743, 751, 760, 
	768, 776, 784, 791, 798, 805, 812, 819, 
	826, 832, 838, 845, 851, 857, 863, 868, 
	874, 880, 885, 891, 896, 901, 907, 912, 
	917, 922, 927, 932, 937, 941, 946, 951, 
	955, 960, 965, 969, 974, 978, 982, 987, 
	991, 995, 999, 1004, 1008, 1012, 1016, 1020, 
	1024,
    };

// Delay Agnostic AEC parameters, still under development and may change.
//static const float kDelayQualityThresholdMax = 0.07f;
//static const float kDelayQualityThresholdMin = 0.01f;
static const int kInitialShiftOffset = 5;
#if !defined(WEBRTC_ANDROID)
static const int kDelayCorrectionStart = 1500;  // 10 ms chunks
#endif

//Twiddles for fft/ifft 2*N*3/4  N=128
//static float twd_128[192] = {1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 0.998795, -0.049068, 0.995185, -0.098017, 0.989177, -0.146730, 0.995185, -0.098017, 0.980785, -0.195090, 0.956940, -0.290285, 0.989177, -0.146730, 0.956940, -0.290285, 0.903989, -0.427555, 0.980785, -0.195090, 0.923880, -0.382683, 0.831470, -0.555570, 0.970031, -0.242980, 0.881921, -0.471397, 0.740951, -0.671559, 0.956940, -0.290285, 0.831470, -0.555570, 0.634393, -0.773010, 0.941544, -0.336890, 0.773010, -0.634393, 0.514103, -0.857729, 0.923880, -0.382683, 0.707107, -0.707107, 0.382683, -0.923880, 0.903989, -0.427555, 0.634393, -0.773010, 0.242980, -0.970031, 0.881921, -0.471397, 0.555570, -0.831470, 0.098017, -0.995185, 0.857729, -0.514103, 0.471397, -0.881921, -0.049068, -0.998795, 0.831470, -0.555570, 0.382683, -0.923880, -0.195090, -0.980785, 0.803208, -0.595699, 0.290285, -0.956940, -0.336890, -0.941544, 0.773010, -0.634393, 0.195090, -0.980785, -0.471397, -0.881921, 0.740951, -0.671559, 0.098017, -0.995185, -0.595699, -0.803208, 0.707107, -0.707107, 0.000000, -1.000000, -0.707107, -0.707107, 0.671559, -0.740951, -0.098017, -0.995185, -0.803208, -0.595699, 0.634393, -0.773010, -0.195090, -0.980785, -0.881921, -0.471397, 0.595699, -0.803208, -0.290285, -0.956940, -0.941544, -0.336890, 0.555570, -0.831470, -0.382683, -0.923880, -0.980785, -0.195090, 0.514103, -0.857729, -0.471397, -0.881921, -0.998795, -0.049068, 0.471397, -0.881921, -0.555570, -0.831470, -0.995185, 0.098017, 0.427555, -0.903989, -0.634393, -0.773010, -0.970031, 0.242980, 0.382683, -0.923880, -0.707107, -0.707107, -0.923880, 0.382683, 0.336890, -0.941544, -0.773010, -0.634393, -0.857729, 0.514103, 0.290285, -0.956940, -0.831470, -0.555570, -0.773010, 0.634393, 0.242980, -0.970031, -0.881921, -0.471397, -0.671559, 0.740951, 0.195090, -0.980785, -0.923880, -0.382683, -0.555570, 0.831470, 0.146730, -0.989177, -0.956940, -0.290285, -0.427555, 0.903989, 0.098017, -0.995185, -0.980785, -0.195090, -0.290285, 0.956940, 0.049068, -0.998795, -0.995185, -0.098017, -0.146730, 0.989177, };


// Target suppression levels for nlp modes.
// log{0.001, 0.00001, 0.00000001}
//static const float kTargetSupp1[3] = {-6.9f, -11.5f, -18.4f};
static const int kTargetSupp[3] = {-231525581, -385875968, -617401549};//Q25

// Two sets of parameters, one for the extended filter mode.
static const int kExtendedMinOverDrive[3] = {1536, 3072, 7680};//Q9
static const int kNormalMinOverDrive[3] = {512, 1024, 2560};
const long long WebRtcAec_kExtendedSmoothingCoefficients[2][2] = {{7373, 819},{7537, 655}};//{{0.9f, 0.1f},{0.92f, 0.08f}}; //Q13
const long long WebRtcAec_kNormalSmoothingCoefficients[2][2] =   {{7373, 819},{7619, 573}};//{{0.9f, 0.1f},{0.93f, 0.07f}};
const float fWebRtcAec_kExtendedSmoothingCoefficients[2][2] = {{0.9f, 0.1f},{0.92f, 0.08f}}; //Q13
const float fWebRtcAec_kNormalSmoothingCoefficients[2][2] =   {{0.9f, 0.1f},{0.93f, 0.07f}};
// Number of partitions forming the NLP's "preferred" bands.
enum { kPrefBandSize = 24 };

//WebRtcAecFilterFar WebRtcAec_FilterFar;
//WebRtcAecScaleErrorSignal WebRtcAec_ScaleErrorSignal;
//WebRtcAecFilterAdaptation WebRtcAec_FilterAdaptation;
//WebRtcAecOverdrive WebRtcAec_Overdrive;
//WebRtcAecSuppress WebRtcAec_Suppress;
//WebRtcAecComputeCoherence WebRtcAec_ComputeCoherence;
//WebRtcAecUpdateCoherenceSpectra WebRtcAec_UpdateCoherenceSpectra;
//WebRtcAecStoreAsComplex WebRtcAec_StoreAsComplex;
//WebRtcAecPartitionDelay WebRtcAec_PartitionDelay;
//WebRtcAecWindowData WebRtcAec_WindowData;

// TODO(minyue): Due to a legacy bug, |framelevel| and |averagelevel| use a
// window, of which the length is 1 unit longer than indicated. Remove "+1" when
// the code is refactored.
RingBuffer* buffer_;
ReInit() {
  WebRtc_InitBuffer(buffer_);
}
Insert(const int32_t block[PART_LEN]) {
  WebRtc_WriteBuffer(buffer_, block, 1);
}
ExtractExtendedBlock(int32_t extended_block[PART_LEN2]) {
  int32_t* block_ptr = NULL;
  int32_t* extend_ptr = NULL;
  //assert(0 < AvaliableSpace());
//int32_t intbuf[PART_LEN2];
  // Extract the previous block.
  WebRtc_MoveReadPtr(buffer_, -1);
  extend_ptr = (int32_t *)&extended_block[0];
  size_t read_elements = WebRtc_ReadBuffer(
      buffer_, (void**)(&block_ptr),extend_ptr, 1);  
  if (read_elements == 0u) {
    memset(&extended_block[0], 0, PART_LEN);
  } else if (block_ptr != extend_ptr) {
    #if 1
    for(int i=0; i<PART_LEN; i++)
    {
        extended_block[i] = block_ptr[i];
    }
    #else
    memcpy(&extended_block[0], block_ptr, PART_LEN * sizeof(int32_t));
    #endif
  }

  // Extract the current block.
  extend_ptr = (int32_t *)&extended_block[PART_LEN];
  read_elements =
      WebRtc_ReadBuffer(buffer_, (void**)(&block_ptr),
                        extend_ptr, 1);
  if (read_elements == 0u) {
    memset(&extended_block[PART_LEN], 0, PART_LEN);
  } else if (block_ptr != extend_ptr) {
    #if 1
    for(int i=0; i<PART_LEN; i++)
    {
        extended_block[PART_LEN+i] = block_ptr[i];
    }
    #else
    memcpy(&extended_block[PART_LEN], block_ptr, PART_LEN * sizeof(int32_t));
    #endif
  }
}

int AdjustSize(int buffer_size_decrease) {
  return WebRtc_MoveReadPtr(buffer_, buffer_size_decrease);
}

int Size() {
  return (int)(WebRtc_available_read(buffer_));
}

int AvaliableSpace() {
  return WebRtc_available_write(buffer_);
}
#if 0
DivergentFilterFraction::DivergentFilterFraction()
    : count_(0), occurrence_(0), fraction_(-1.0) {}

void DivergentFilterFraction::Reset() {
  Clear();
  fraction_ = -1.0;
}

void DivergentFilterFraction::AddObservation(const PowerLevel& nearlevel,
                                             const PowerLevel& linoutlevel,
                                             const PowerLevel& nlpoutlevel) {
  const float near_level = nearlevel.framelevel.GetLatestMean();
  const float level_increase =
      linoutlevel.framelevel.GetLatestMean() - near_level;
  const bool output_signal_active =
      nlpoutlevel.framelevel.GetLatestMean() > 40.0 * nlpoutlevel.minlevel;
  // Level increase should be, in principle, negative, when the filter
  // does not diverge. Here we allow some margin (0.01 * near end level) and
  // numerical error (1.0). We count divergence only when the AEC output
  // signal is active.
  if (output_signal_active && level_increase > std::max(0.01 * near_level, 1.0))
    occurrence_++;
  ++count_;
  if (count_ == kDivergentFilterFractionAggregationWindowSize) {
    fraction_ = static_cast<float>(occurrence_) /
                kDivergentFilterFractionAggregationWindowSize;
    Clear();
  }
}

float DivergentFilterFraction::GetLatestFraction() const {
  return fraction_;
}
#endif

// TODO(minyue): Moving some initialization from WebRtcAec_CreateAec() to ctor.
//AecCore::AecCore(int instance_index) {}

//AecCore::~AecCore() {}

static int CmpFloat(const void* a, const void* b) {
  const int* da = (const int*)a;
  const int* db = (const int*)b;

  //return (*da > *db) - (*da < *db);
  return (*da - *db);
}
#ifdef HARD_FFT
#define FFTN 128 //128 512
void hard_rfft_comp(int fft_flag, int *fftin, int *fftout)
{
    //in:先x1(n)的128点，再x2(n)的128点. result: 先x1(n)的复数65点，再x2(n)的复数65点
    if(fft_flag == 0)//fft
    {
        fft_struct.addr = &g_fft_buf[0];	
    	fft_struct.len = FFTN * 2;
    	#if 1
    	for(int j = 0; j < FFTN; j++)
        {
    		int k=j<<1;
    		g_fft_buf[k] = fftin[j];
    		g_fft_buf[k+1] = fftin[FFTN+j];
        }	
        #endif
        ecrd = iet_fft_init(FFT_POLL);
		if(ecrd != 0)
	    {
	        iet_printf_err("fft_poll iet_fft_init = 0x%x\r\n", ecrd);
	    }	
		ecrd = iet_fft_control(FFT_CMD_WRITE, &fft_struct);
        if(ecrd != E_OK)
	    {
	        iet_printf_err("fft_poll FFT_CMD_WRITE = 0x%x\r\n", ecrd);
	    }
	    ecrd = iet_fft_start();
		if(ecrd != E_OK)
	    {
	        iet_printf_err("fft_poll iet_fft_start = 0x%x\r\n", ecrd);
	    }
	    while(!iet_fft_if_done())
	    {
	    }

	    fftout[0] = g_fft_buf[0];           fftout[1] = 0;
	    fftout[FFTN+2+0] = g_fft_buf[1];       fftout[FFTN+2+1] = 0;
	    fftout[FFTN] = g_fft_buf[FFTN+1];     fftout[FFTN+1] = 0;
	    fftout[FFTN+2+FFTN] = g_fft_buf[FFTN+1]; fftout[FFTN+2+FFTN+1] = 0;
	    for(int j = 1; j < FFTN/2; j++)
        {
            int k=j<<1;
    		fftout[k] =   (g_fft_buf[k] + g_fft_buf[FFTN*2-k])>>1;
    		fftout[k+1] = (g_fft_buf[k+1] - g_fft_buf[FFTN*2-k+1])>>1;

    		fftout[FFTN+2+k] =   (g_fft_buf[k+1] + g_fft_buf[FFTN*2-k+1])>>1;
    		//fftout[FFTN+2+k+1] = (g_fft_buf[k] - g_fft_buf[FFTN*2-k])>>1;//公式是否漏了负号
    		fftout[FFTN+2+k+1] = (g_fft_buf[FFTN*2-k] - g_fft_buf[k])>>1;
        }
    }
    
    //in:先x1(n)的复数65点，再x2(n)的复数65点. result: 先x1(n)的128点，再x2(n)的128点
    else if(fft_flag == 1)//ifft 
    {
        fft_struct.addr = &g_fft_buf[0];	
    	fft_struct.len = FFTN * 2;

    	for(int j = 0; j < FFTN/2+1; j++)
    	{
    	    int k=j<<1;
    	    g_fft_buf[k] = fftin[k]-fftin[FFTN+2+k+1];
    	    g_fft_buf[k+1] = fftin[k+1]+fftin[FFTN+2+k];
    	}
    	for(int j = FFTN/2+1; j < FFTN; j++)
    	{
            int k=j<<1;
            g_fft_buf[k] = fftin[(FFTN-j)*2]+fftin[FFTN+2+(FFTN-j)*2+1];//注意，数据后半段不能直接取共轭
            g_fft_buf[k+1] = -fftin[(FFTN-j)*2+1]+fftin[FFTN+2+(FFTN-j)*2];
         }
        
        ecrd = iet_fft_init(FFT_POLL);
        if(ecrd != 0)
        {
            iet_printf_err("ifft_poll iet_fft_init = 0x%x\r\n", ecrd);
        }
    	ecrd = iet_fft_control(FFT_CMD_WRITE, &fft_struct);
        if(ecrd != E_OK)
        {
            iet_printf_err("ifft_poll FFT_CMD_WRITE = 0x%x\r\n", ecrd);
        }
        ecrd = iet_ifft_start();
        if(ecrd != E_OK)
        {
            iet_printf_err("ifft_poll iet_fft_start = 0x%x\r\n", ecrd);
        }
    	while(!iet_fft_if_done())
        {
        }
    	for(int j = 0; j < FFTN; j++)
        {
        	int k=j<<1;
    		fftout[j] = g_fft_buf[k];
    		fftout[FFTN+j] = g_fft_buf[k+1];
           
        }
    }
}
void hard_rfft_int32(int fft_flag, int *fftin, int *fftout)
{
    if(fft_flag == 0)//fft
    {
        fft_struct.addr = &g_fft_buf[0];	
    	fft_struct.len = FFTN * 2;
    	#if 1
    	for(int j = 0; j < FFTN; j++)
        {
    		int k=j<<1;
    		g_fft_buf[k] = fftin[j];
    		g_fft_buf[k+1] = 0;
        }	
        #endif
        ecrd = iet_fft_init(FFT_POLL);
		if(ecrd != 0)
	    {
	        iet_printf_err("fft_poll iet_fft_init = 0x%x\r\n", ecrd);
	    }	
		
		ecrd = iet_fft_control(FFT_CMD_WRITE, &fft_struct);
	    if(ecrd != E_OK)
	    {
	        iet_printf_err("fft_poll FFT_CMD_WRITE = 0x%x\r\n", ecrd);
	    }
	    ecrd = iet_fft_start();
	    if(ecrd != E_OK)
	    {
	        iet_printf_err("fft_poll iet_fft_start = 0x%x\r\n", ecrd);
	    }
		
	    while(!iet_fft_if_done())
	    {
	    }
	    #if 1
	    for(int j = 0; j < (FFTN+2); j++)
        {
    		fftout[j] = g_fft_buf[j];
        }
        #endif
    }
    else if(fft_flag == 1)//ifft
    {
        fft_struct.addr = &g_fft_buf[0];	
    	fft_struct.len = FFTN * 2;
    	#if 1
    	memcpy(g_fft_buf, fftin, sizeof(int)*(FFTN+2));
    	for(int j = FFTN/2+1; j < FFTN; j++)
    	{
            int k=j<<1;
            g_fft_buf[k] = g_fft_buf[FFTN*2-k];
            g_fft_buf[k+1] = -g_fft_buf[FFTN*2-k+1];
         }
#endif
        
        ecrd = iet_fft_init(FFT_POLL);
    	if(ecrd != 0)
        {
            iet_printf_err("ifft_poll iet_fft_init = 0x%x\r\n", ecrd);
        }	
    	ecrd = iet_fft_control(FFT_CMD_WRITE, &fft_struct);
        if(ecrd != E_OK)
        {
            iet_printf_err("ifft_poll FFT_CMD_WRITE = 0x%x\r\n", ecrd);
        }
        ecrd = iet_ifft_start();
        if(ecrd != E_OK)
        {
            iet_printf_err("ifft_poll iet_fft_start = 0x%x\r\n", ecrd);
        }
    	while(!iet_fft_if_done())
        {
        }
        #if 1
    	for(int j = 0; j < FFTN; j++)
        {
        	int k=j<<1;
    		fftout[j] = g_fft_buf[k];
           
        }
        #endif
    }
}
#endif
void FilterFar(int num_partitions,
                      int x_fft_buf_block_pos,
                      int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//Q0
                      int32_t h_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//aecq=-15
                      int32_t y_fft_hifi[PART_LEN3]//aecq0
                                     ) {
  int i,j;
  int shift=aecq0-aecq;
  for (i = 0; i < num_partitions; i++) {
      int xPos = (i + x_fft_buf_block_pos) * PART_LEN3;
      int pos = i * PART_LEN3;
      // Check for wrap
      if (i + x_fft_buf_block_pos >= num_partitions) {
        xPos -= num_partitions * (PART_LEN3);
      }

      for (j = 0; j < PART_LEN1; j++) 
      {//(a+bi)(c+di)
          y_fft_hifi[2*j] += ((long long)x_fft_buf_hifi[xPos + j*2] *h_fft_buf_hifi[pos + j*2] - (long long)x_fft_buf_hifi[xPos +j*2 + 1] * h_fft_buf_hifi[pos + j*2 + 1])>>(shift);
          y_fft_hifi[2*j + 1] += ((long long)x_fft_buf_hifi[xPos +j*2 + 1] * h_fft_buf_hifi[pos + j*2]+ (long long)x_fft_buf_hifi[xPos + j*2] * h_fft_buf_hifi[pos + j*2 + 1])>>(shift);
       }
  }
}
#if 0
#if 0
void ScaleErrorSignal(int mu,//Q31
                             int error_threshold,//Q31
                             long long x_pow[PART_LEN1],//Q0
                             int32_t ef_hifi[PART_LEN3],//in: Q0, out: Q31
                             long long ef_hifi_err[PART_LEN3]
) {
  int i;
  float ferror_threshold=1e-6f; 
  float abs_ef;
  float ef0, ef1;
  for (i = 0; i < (PART_LEN1); i++) {
    ef0 = (float)ef_hifi[2*i] /(x_pow[i] + 1e-10f);
    ef1 = (float)ef_hifi[2*i+1]/(x_pow[i] + 1e-10f);
    abs_ef = sqrtf(ef0 * ef0 + ef1 * ef1);

    if (abs_ef > ferror_threshold) {
      abs_ef = ferror_threshold / (abs_ef + 1e-10f);
      ef0 *= abs_ef;
      ef1 *= abs_ef;
    }

    // Stepsize factor
    ef0 = ef0* 0.4*2147483648;
    ef1 = ef1* 0.4*2147483648;
    ef_hifi_err[2*i] = ef0;
    ef_hifi_err[2*i+1]= ef1;
  }
}
#else
void ScaleErrorSignal(int mu,//Q31
                             int error_threshold,//Q31
                             long long x_pow[PART_LEN1],//Q0
                             int32_t ef_hifi[PART_LEN3],//in: Q0, out: Q31
                             long long ef_hifi_err[PART_LEN3]
) {
  int i;
  long long abs_ef;
  long long thres;
  //float fix_ef;
  long long thr_tmp = ((long long)error_threshold*mu);//Q62
  aecnearq = -31;
  for (i = 0; i < (PART_LEN1); i++) {
	  abs_ef = sqrtf((long long)ef_hifi[2*i] * ef_hifi[2*i] + (long long)ef_hifi[2*i+1] * ef_hifi[2*i+1]);
      thres = ((long long)x_pow[i]*error_threshold)>>31;
      if (abs_ef > thres) {
        ef_hifi_err[2*i] =    (((long long)ef_hifi[2*i]*thr_tmp) / (abs_ef))>>31;
        ef_hifi_err[2*i+1] = (((long long)ef_hifi[2*i+1]*thr_tmp) / (abs_ef))>>31;
        //fix_ef = (abs_ef*mu);// Q31
      }
      else
      {
          // Stepsize factor
          ef_hifi_err[2*i] = (((long long)ef_hifi[2*i]*mu)/(x_pow[i]+1));
          ef_hifi_err[2*i+1] = (((long long)ef_hifi[2*i+1]*mu)/(x_pow[i]+1));
      }

    }

}
#endif
#else
#define max(a, b)      ((a) > (b) ? (a) : (b))
#define min(a, b)      ((a) < (b) ? (a) : (b))
#define MAX_MOD(a, b) (max(abs(a), abs(b))+(min(abs(a),abs(b))>>2))
//#define MAX_MOD(a, b) (((max(abs(a), abs(b))<<2)+min(abs(a),abs(b)))>>2)
void ScaleErrorSignal(int mu,//Q31
                             int error_threshold,//Q31
                             long long x_pow[PART_LEN1],//Q0
                             int32_t ef_hifi[PART_LEN3],//in: Q0, out: Q31
                             long long ef_hifi_err[PART_LEN3]
) {
  int i;
  uint32_t abs_ef;
  long long thres;
  long long thr_tmp = ((long long)error_threshold*mu);//Q62
  aecnearq = -31;
  for (i = 0; i < (PART_LEN1); i++) {
	  abs_ef = MAX_MOD(ef_hifi[2*i], ef_hifi[2*i+1]);//sqrtf((long long)ef_hifi[2*i] * ef_hifi[2*i] + (long long)ef_hifi[2*i+1] * ef_hifi[2*i+1]);
      thres = ((long long)x_pow[i]*error_threshold)>>31;
      if (abs_ef > thres) {
        ef_hifi_err[2*i] =    (((long long)ef_hifi[2*i]*thr_tmp) / (abs_ef))>>31;
        ef_hifi_err[2*i+1] = (((long long)ef_hifi[2*i+1]*thr_tmp) / (abs_ef))>>31;
      }
      else
      {
          // Stepsize factor
          ef_hifi_err[2*i] = (((long long)ef_hifi[2*i]*mu)/(x_pow[i]+1));
          ef_hifi_err[2*i+1] = (((long long)ef_hifi[2*i+1]*mu)/(x_pow[i]+1));
      }

    }

}
#endif
 #ifdef HARD_FFT

typedef uint32_t (*asr_wait_done)(uint32_t ui4_timeout);
extern asr_wait_done  fft_finsh;

void AEC_fft(int* pData)
{
	fft_struct.addr = pData;
	fft_struct.len = 128 * 2;
	
	ecrd = iet_fft_control(FFT_CMD_WRITE, &fft_struct);
    if(ecrd != E_OK)
    {
        iet_printf_err("fft_poll FFT_CMD_WRITE = 0x%x\r\n", ecrd);
    }

    ecrd = iet_fft_start();
    if(ecrd != E_OK)
    {
        iet_printf_err("fft_poll iet_fft_start = 0x%x\r\n", ecrd);
    }

	
	int result = fft_finsh(1);
	if(result != 0){
		iet_printf_err("fft_finsh time out!!!!!!\r\n", ecrd);
	}
}


void hard_rfft_comp_INT(int fft_flag, int *fftin, int *fftout)
{
    //in:先x1(n)的128点，再x2(n)的128点. result: 先x1(n)的复数65点，再x2(n)的复数65点
    if(fft_flag == 0)//fft
    {
        //fft_struct.addr = &g_fft_buf[0];	
    	//fft_struct.len = 128 * 2;
    	#if 1
    	for(int j = 0; j < 128; j++)
        {
    		int k=j<<1;
    		g_fft_buf[k] = fftin[j];
    		g_fft_buf[k+1] = fftin[128+j];
        }	
        #endif
        AEC_fft(g_fft_buf);
        
	    fftout[0] = g_fft_buf[0];           fftout[1] = 0;
	    fftout[130+0] = g_fft_buf[1];       fftout[130+1] = 0;
	    fftout[128] = g_fft_buf[128+1];     fftout[128+1] = 0;
	    fftout[130+128] = g_fft_buf[128+1]; fftout[130+128+1] = 0;
	    for(int j = 1; j < 63; j++)
        {
            int k=j<<1;
    		fftout[k] =   (g_fft_buf[k] + g_fft_buf[256-k])>>1;
    		fftout[k+1] = (g_fft_buf[k+1] - g_fft_buf[256-k+1])>>1;

    		fftout[130+k] =   (g_fft_buf[k+1] + g_fft_buf[256-k+1])>>1;
    		fftout[130+k+1] = (g_fft_buf[k] - g_fft_buf[256-k])>>1;
        }
    }
    
    //in:先x1(n)的复数65点，再x2(n)的复数65点. result: 先x1(n)的128点，再x2(n)的128点
    else if(fft_flag == 1)//ifft 
    {
        //fft_struct.addr = &g_fft_buf[0];	
    	//fft_struct.len = 128 * 2;

    	for(int j = 0; j < 65; j++)
    	{
    	    int k=j<<1;
    	    g_fft_buf[k] = fftin[k]-fftin[130+k+1];
    	    g_fft_buf[k+1] = fftin[k+1]-fftin[130+k];
    	}
    	for(int j = 65; j < 128; j++)
    	{
            int k=j<<1;
            g_fft_buf[k] = g_fft_buf[256-k];
            g_fft_buf[k+1] = -g_fft_buf[256-k+1];
         }
        AEC_fft(g_fft_buf);
    	for(int j = 0; j < 128; j++)
        {
        	int k=j<<1;
    		fftout[j] = g_fft_buf[k];
    		fftout[128+j] = g_fft_buf[k+1];
           
        }
    }
}
void hard_rfft_int32_INT(int fft_flag, int *fftin, int *fftout)

{
    if(fft_flag == 0)//fft
    {
        //fft_struct.addr = &g_fft_buf[0];	
    	//fft_struct.len = 128 * 2;
    	#if 1
    	for(int j = 0; j < 128; j++)
        {
    		int k=j<<1;
    		g_fft_buf[k] = fftin[j];
    		g_fft_buf[k+1] = 0;
        }	
        #endif
        AEC_fft(g_fft_buf);

	    for(int j = 0; j < 130; j++)
        {
    		fftout[j] = g_fft_buf[j];
        }
    }
    else if(fft_flag == 1)//ifft
    {
        //fft_struct.addr = &g_fft_buf[0];	
    	//fft_struct.len = 128 * 2;
    	memcpy(g_fft_buf, fftin, sizeof(int)*130);
    	for(int j = 65; j < 128; j++)
    	{
            int k=j<<1;
            g_fft_buf[k] = g_fft_buf[256-k];
            g_fft_buf[k+1] = -g_fft_buf[256-k+1];
         }
        AEC_fft(g_fft_buf);
    	for(int j = 0; j < 128; j++)
        {
        	int k=j<<1;
    		fftout[j] = g_fft_buf[k];
           
        }
    }
}

 
void FilterAdaptation_comp(
    int num_partitions,
    int x_fft_buf_block_pos,
    int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//Q0
    int32_t e_fft_hifi[PART_LEN3],//Q31
    int32_t h_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3]//aecq
                       ) {
  int i, j;
  int32_t fft[PART_LEN2*2] __attribute__((aligned(16)));
  int32_t fft_hifi[PART_LEN3*2];
  int32_t fft_out[PART_LEN3*2];
  for (i = 0; i < num_partitions>>1; i++) {
     int xPos1 = (i*2 + x_fft_buf_block_pos) * (PART_LEN3);
     int pos1;
     int xPos2 = (i*2+1 + x_fft_buf_block_pos) * (PART_LEN3);
     int pos2;

     // Check for wrap
     if (i*2 + x_fft_buf_block_pos >= num_partitions) {
       xPos1 -= num_partitions * PART_LEN3;
     }
     if (i*2+1 + x_fft_buf_block_pos >= num_partitions) {
       xPos2 -= num_partitions * PART_LEN3;
     }

    pos1 = i*2 * PART_LEN3;
    pos2 = (i*2+1) * PART_LEN3;
    //1. in frequency domain, do inversefft for mu*e*x ; 2. then reserve part of result and scale fft in time domain; 3. transform the processed fft to frequency domain and update the filer coefficient
    #if 1
    for (j = 0; j < PART_LEN1; j++){ //(a-bi)*(c+di)
        int32_t err,eri;
        int k=j<<1;
        err = e_fft_hifi[k];
        eri = e_fft_hifi[k+1];
        fft_hifi[k] = (x_fft_buf_hifi[xPos1 + k]*err + x_fft_buf_hifi[xPos1 + k + 1]*eri)>>7;
        fft_hifi[k + 1] = (x_fft_buf_hifi[xPos1 + k]*eri-x_fft_buf_hifi[xPos1 + k + 1]*err)>>7;
        fft_hifi[PART_LEN3+k] = (x_fft_buf_hifi[xPos2 + k]*err + x_fft_buf_hifi[xPos2 + k + 1]*eri)>>7;
        fft_hifi[PART_LEN3+k + 1] = (x_fft_buf_hifi[xPos2 + k]*eri-x_fft_buf_hifi[xPos2 + k + 1]*err)>>7;
    }
    #ifndef INT_MODE 
    hard_rfft_comp(1, fft_hifi, fft);//hyli
    #else
    hard_rfft_comp_INT(1, fft_hifi, fft);
    #endif

    memset(fft + PART_LEN, 0, sizeof(int) * PART_LEN);
    memset(fft+PART_LEN2 + PART_LEN, 0, sizeof(int) * PART_LEN);

    #ifndef INT_MODE 
    hard_rfft_comp(0, fft, fft_out);//hyli
    #else
    hard_rfft_comp_INT(0, fft, fft_out);
    #endif

    for (j = 0; j < PART_LEN1; j++){
        int k=j<<1;
    	h_fft_buf_hifi[pos1 + k]     += (fft_out[k]>>(9));
    	h_fft_buf_hifi[pos1 + k + 1] += (fft_out[k + 1]>>(9));
    	h_fft_buf_hifi[pos2 + k]     += (fft_out[PART_LEN3+k]>>(9));
    	h_fft_buf_hifi[pos2 + k + 1] += (fft_out[PART_LEN3+k + 1]>>(9));
    }
    #else
    for (j = 0; j < PART_LEN1; j++){
        int32_t err,eri;
        int k=j<<1;
        err = e_fft_hifi[k];
        eri = e_fft_hifi[k+1];
    	h_fft_buf_hifi[pos1 + k]     += ((x_fft_buf_hifi[xPos1 + k]*err + x_fft_buf_hifi[xPos1 + k + 1]*eri)>>(16));
    	h_fft_buf_hifi[pos1 + k + 1] += ((x_fft_buf_hifi[xPos1 + k]*eri-x_fft_buf_hifi[xPos1 + k + 1]*err)>>(16));
    	h_fft_buf_hifi[pos2 + k]     += ((x_fft_buf_hifi[xPos2 + k]*err + x_fft_buf_hifi[xPos2 + k + 1]*eri)>>(16));
    	h_fft_buf_hifi[pos2 + k + 1] += ((x_fft_buf_hifi[xPos2 + k]*eri-x_fft_buf_hifi[xPos2 + k + 1]*err)>>(16));
    }
    #endif

  }

}
void FilterAdaptation_comp_parallel(
    int num_partitions,
    int x_fft_buf_block_pos,
    int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//Q0
    int32_t e_fft_hifi[PART_LEN3],//Q31
    int32_t h_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3]//aecq
                       )
{
  int i, j;
  int32_t fft[PART_LEN2*2] __attribute__((aligned(16)));
  int32_t fft_hifi[PART_LEN3*2];
  int32_t fft_out[PART_LEN3*2];
  int32_t fft_hifi2[PART_LEN3*2];
  for (i = 0; i < num_partitions>>2; i++) {
     int xPos1 = (i*4 + x_fft_buf_block_pos) * (PART_LEN3);
     int pos1;
     int xPos2 = (i*4+1 + x_fft_buf_block_pos) * (PART_LEN3);
     int pos2;

     // Check for wrap
     if (i*4 + x_fft_buf_block_pos >= num_partitions) {
       xPos1 -= num_partitions * PART_LEN3;
     }
     if (i*4+1 + x_fft_buf_block_pos >= num_partitions) {
       xPos2 -= num_partitions * PART_LEN3;
     }

    pos1 = i*4 * PART_LEN3;
    pos2 = (i*4+1) * PART_LEN3;
    //1. in frequency domain, do inversefft for mu*e*x ; 2. then reserve part of result and scale fft in time domain; 3. transform the processed fft to frequency domain and update the filer coefficient

    for (j = 0; j < PART_LEN1; j++){ //(a-bi)*(c+di)
        int32_t err,eri;
        int k=j<<1;
        err = e_fft_hifi[k];
        eri = e_fft_hifi[k+1];
        fft_hifi[k] = (x_fft_buf_hifi[xPos1 + k]*err + x_fft_buf_hifi[xPos1 + k + 1]*eri)>>7;
        fft_hifi[k + 1] = (x_fft_buf_hifi[xPos1 + k]*eri-x_fft_buf_hifi[xPos1 + k + 1]*err)>>7;
        fft_hifi[PART_LEN3+k] = (x_fft_buf_hifi[xPos2 + k]*err + x_fft_buf_hifi[xPos2 + k + 1]*eri)>>7;
        fft_hifi[PART_LEN3+k + 1] = (x_fft_buf_hifi[xPos2 + k]*eri-x_fft_buf_hifi[xPos2 + k + 1]*err)>>7;
    }

    hard_rfft_comp_INT(1, fft_hifi, fft);
    
//////
    int xPos3 = (i*4+2 + x_fft_buf_block_pos) * (PART_LEN3);
    int pos3;
    int xPos4 = (i*4+3 + x_fft_buf_block_pos) * (PART_LEN3);
    int pos4;

    // Check for wrap
    if (i*4+2 + x_fft_buf_block_pos >= num_partitions) {
        xPos3 -= num_partitions * PART_LEN3;
    }
    if (i*4+3 + x_fft_buf_block_pos >= num_partitions) {
        xPos4 -= num_partitions * PART_LEN3;
    }

    pos3 = (i*4+2) * PART_LEN3;
    pos4 = (i*4+3) * PART_LEN3;

    for (j = 0; j < PART_LEN1; j++){ //(a-bi)*(c+di)
        int32_t err,eri;
        int k=j<<1;
        err = e_fft_hifi[k];
        eri = e_fft_hifi[k+1];
        fft_hifi2[k] = (x_fft_buf_hifi[xPos3 + k]*err + x_fft_buf_hifi[xPos3 + k + 1]*eri)>>7;
        fft_hifi2[k + 1] = (x_fft_buf_hifi[xPos3 + k]*eri-x_fft_buf_hifi[xPos3 + k + 1]*err)>>7;
        fft_hifi2[PART_LEN3+k] = (x_fft_buf_hifi[xPos4 + k]*err + x_fft_buf_hifi[xPos4 + k + 1]*eri)>>7;
        fft_hifi2[PART_LEN3+k + 1] = (x_fft_buf_hifi[xPos4 + k]*eri-x_fft_buf_hifi[xPos4 + k + 1]*err)>>7;
    }
//////

    memset(fft + PART_LEN, 0, sizeof(int) * PART_LEN);
    memset(fft+PART_LEN2 + PART_LEN, 0, sizeof(int) * PART_LEN);

    hard_rfft_comp_INT(0, fft, fft_out);

    for (j = 0; j < PART_LEN1; j++){
        int k=j<<1;
    	h_fft_buf_hifi[pos1 + k]     += (fft_out[k]>>(9));
    	h_fft_buf_hifi[pos1 + k + 1] += (fft_out[k + 1]>>(9));
    	h_fft_buf_hifi[pos2 + k]     += (fft_out[PART_LEN3+k]>>(9));
    	h_fft_buf_hifi[pos2 + k + 1] += (fft_out[PART_LEN3+k + 1]>>(9));
    }

//////    
    hard_rfft_comp_INT(1, fft_hifi2, fft);
    
    memset(fft + PART_LEN, 0, sizeof(int) * PART_LEN);
    memset(fft+PART_LEN2 + PART_LEN, 0, sizeof(int) * PART_LEN);
    
    hard_rfft_comp_INT(0, fft, fft_out);

    for (j = 0; j < PART_LEN1; j++){
        int k=j<<1;
    	h_fft_buf_hifi[pos3 + k]     += (fft_out[k]>>(9));
    	h_fft_buf_hifi[pos3 + k + 1] += (fft_out[k + 1]>>(9));
    	h_fft_buf_hifi[pos4 + k]     += (fft_out[PART_LEN3+k]>>(9));
    	h_fft_buf_hifi[pos4 + k + 1] += (fft_out[PART_LEN3+k + 1]>>(9));
    }
  }

}
#endif
void FilterAdaptation(
    int num_partitions,
    int x_fft_buf_block_pos,
    int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//Q0
    long long e_fft_hifi[PART_LEN3],//Q31
    int32_t h_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3]//aecq
                       ) {
  int i, j;
  for (i = 0; i < num_partitions; i++) {
     int xPos = (i + x_fft_buf_block_pos) * (PART_LEN3);
     int pos;

     // Check for wrap
     if (i + x_fft_buf_block_pos >= num_partitions) {
       xPos -= num_partitions * PART_LEN3;
     }

     pos = i * PART_LEN3;
    //1. in frequency domain, do inversefft for mu*e*x ; 2. then reserve part of result and scale fft in time domain; 3. transform the processed fft to frequency domain and update the filer coefficient
    int32_t fft[PART_LEN2*2] __attribute__((aligned(16)));
    int32_t fft_hifi[PART_LEN3*2];
    int32_t fft_out[PART_LEN3*2];
    #if 0
    for (j = 0; j < PART_LEN1; j++){ //(a-bi)*(c+di)
        fft_hifi[2*j] = ((long long)x_fft_buf_hifi[xPos + j*2]*e_fft_hifi[2*j] +(long long) x_fft_buf_hifi[xPos + j*2 + 1]*e_fft_hifi[2*j+1])>>1;
        fft_hifi[2*j + 1] = ((long long)x_fft_buf_hifi[xPos + j*2]*e_fft_hifi[2*j+1]-(long long)x_fft_buf_hifi[xPos + j*2 + 1]*e_fft_hifi[2*j])>>1;
    }    
    #ifdef HARD_FFT
        #ifndef INT_MODE         
        hard_rfft_int32(1, fft_hifi, fft);//hyli
        #else
        hard_rfft_int32_INT(1, fft_hifi, fft);
        #endif
    #else
    //realFFT_Inverse( fft, fft_hifi, 128 );
    arm_rfft_q31(&S_ifft128, fft_hifi, fft);//hyli
    #endif
    memset(fft + PART_LEN, 0, sizeof(int) * PART_LEN);
    #ifdef HARD_FFT
        #ifndef INT_MODE         
        hard_rfft_int32(0, fft, fft_out);//hyli
        #else
        hard_rfft_int32_INT(0, fft, fft_out);
        #endif
    #else
    //realFFT_forward( fft_out, fft, 128 );
    arm_rfft_q31(&S_fft128, fft, fft_out);//hyli
    #endif
    for (j = 0; j < PART_LEN1; j++){
    	h_fft_buf_hifi[pos + j*2 ]     += (fft_out[2*j]);
    	h_fft_buf_hifi[pos + j*2 + 1] += (fft_out[2*j + 1]);
    }
    #else
    for (j = 0; j < PART_LEN1; j++){
        long long tmp1 = (((long long)x_fft_buf_hifi[xPos + j*2]*e_fft_hifi[2*j] + (long long)x_fft_buf_hifi[xPos + j*2 + 1]*e_fft_hifi[2*j+1])>>(16));
        long long tmp2 = (((long long)x_fft_buf_hifi[xPos + j*2]*e_fft_hifi[2*j+1]-(long long)x_fft_buf_hifi[xPos + j*2 + 1]*e_fft_hifi[2*j])>>(16));
    	h_fft_buf_hifi[pos + j*2 ]     += tmp1;
    	h_fft_buf_hifi[pos + j*2 + 1] += tmp2;
    }
    #endif


  }

}

int log_16q16(int x) {
  int t,y;
  if ( x <= 0 )
        return (-1145600);

  y=0xa65af;
  if(x<0x00008000) x<<=16,              y-=0xb1721;
  if(x<0x00800000) x<<= 8,              y-=0x58b91;
  if(x<0x08000000) x<<= 4,              y-=0x2c5c8;
  if(x<0x20000000) x<<= 2,              y-=0x162e4;
  if(x<0x40000000) x<<= 1,              y-=0x0b172;
  t=x+(x>>1); if((t&0x80000000)==0) x=t,y-=0x067cd;
  t=x+(x>>2); if((t&0x80000000)==0) x=t,y-=0x03920;
  t=x+(x>>3); if((t&0x80000000)==0) x=t,y-=0x01e27;
  t=x+(x>>4); if((t&0x80000000)==0) x=t,y-=0x00f85;
  t=x+(x>>5); if((t&0x80000000)==0) x=t,y-=0x007e1;
  t=x+(x>>6); if((t&0x80000000)==0) x=t,y-=0x003f8;
  t=x+(x>>7); if((t&0x80000000)==0) x=t,y-=0x001fe;
  x=0x80000000-x;
  y-=x>>15;
  return y;
  }


#if 1
static int PartitionDelay(
    int num_partitions,
    int32_t h_fft_buf[kExtendedNumPartitions * PART_LEN3]
                       ) {
  // Measures the energy in each filter partition and returns the partition with
  // highest energy.
  // TODO(bjornv): Spread computational cost by computing one partition per
  // block?
  long long wfEnMax = 0;
  int i;
  int delay = 0;
  for (i = 0; i < num_partitions; i++) {
     int j;
     int pos = i * PART_LEN3;
     long long wfEn = 0;
     for (j = 0; j < PART_LEN1; j++) {
       wfEn += (long long)h_fft_buf[pos + 2*j] *  h_fft_buf[pos + 2*j] +
               (long long)h_fft_buf[pos + 2*j + 1] * h_fft_buf[pos + 2*j + 1];
     }

     if (wfEn > wfEnMax) {
       wfEnMax = wfEn;
       delay = i;
     }
   }
  return delay;
}
#endif
#if 0
// Update metric with 10 * log10(numerator / denominator).
static void UpdateLogRatioMetric(Stats* metric,
                                 float numerator,
                                 float denominator) {
  assert(metric);
  assert(numerator >= 0);
  assert(denominator >= 0);

  const float log_numerator = log10(numerator + 1e-10f);
  const float log_denominator = log10(denominator + 1e-10f);
  metric->instant = 10.0f * (log_numerator - log_denominator);

  // Max.
  if (metric->instant > metric->max)
    metric->max = metric->instant;

  // Min.
  if (metric->instant < metric->min)
    metric->min = metric->instant;

  // Average.
  metric->counter++;
  // This is to protect overflow, which should almost never happen.
  assert(0 != metric->counter);
  metric->sum += metric->instant;
  metric->average = metric->sum / metric->counter;

  // Upper mean.
  if (metric->instant > metric->average) {
    metric->hicounter++;
    // This is to protect overflow, which should almost never happen.
    assert(0 != metric->hicounter);
    metric->hisum += metric->instant;
    metric->himean = metric->hisum / metric->hicounter;
  }
}
#endif

#if 0
// Puts fft output data into a complex valued array.
__inline static void StoreAsComplex(const float* data,
                                    float data_complex[2][PART_LEN1]) {
  int i;
  data_complex[0][0] = data[0];
  data_complex[1][0] = 0;
  for (i = 1; i < PART_LEN; i++) {
    data_complex[0][i] = data[2 * i];
    data_complex[1][i] = data[2 * i + 1];
  }
  data_complex[0][PART_LEN] = data[1];
  data_complex[1][PART_LEN] = 0;
}
#endif
void ComputeCoherence(const CoherenceState* coherence_state,
                             uint16_t* cohde,
                             uint16_t* cohxd) {
  // Subband coherence
  int de,xd;
  for (int i = 0; i < PART_LEN1; i++) {
    long long sxd0 = coherence_state->sxd[i][0]>>2;
    long long sxd1 = coherence_state->sxd[i][1]>>2;
    long long sde0 = coherence_state->sde[i][0];
    long long sde1 = coherence_state->sde[i][1];
    long long sx = coherence_state->sx[i];
    long long sd = coherence_state->sd[i];
    long long se = coherence_state->se[i];
    uint64_t  sdep, sdse, sxdp, sxsd;
    
long long lsd= sxd0 * sxd0+ sxd1 * sxd1;
lsd = (sx * sd);
    sdep = (sde0 * sde0 + sde1 * sde1);
    sdse = (sd * se);
    cohde[i] = ((sdep)/((sdse>>8)+1))<<8;

    sxdp = (sxd0 * sxd0 + sxd1 * sxd1);
    sxsd = (sx * sd);
    cohxd[i] = ((sxdp<<4)/((sxsd>>8)+1))<<8;
    
  }
}
#if 0
static void GetHighbandGain(const float* lambda, float* nlpGainHband) {
  int i;

  *nlpGainHband = 0.0f;
  for (i = freqAvgIc; i < PART_LEN1 - 1; i++) {
    *nlpGainHband += lambda[i];
  }
  *nlpGainHband /= static_cast<float>(PART_LEN1 - 1 - freqAvgIc);
}

static void GenerateComplexNoise(uint32_t* seed, float noise[2][PART_LEN1]) {
  const float kPi2 = 6.28318530717959f;
  int16_t randW16[PART_LEN];
  WebRtcSpl_RandUArray(randW16, PART_LEN, seed);

  noise[0][0] = 0;
  noise[1][0] = 0;
  for (size_t i = 1; i < PART_LEN1; i++) {
    float tmp = kPi2 * randW16[i - 1] / 32768.f;
    noise[0][i] = cosf(tmp);
    noise[1][i] = -sinf(tmp);
  }
  noise[1][PART_LEN] = 0;
}

static void ComfortNoise(bool generate_high_frequency_noise,
                         uint32_t* seed,
                         float e_fft[2][PART_LEN1],
                         float high_frequency_comfort_noise[2][PART_LEN1],
                         const float* noise_spectrum,
                         const float* suppressor_gain) {
  float complex_noise[2][PART_LEN1];

  GenerateComplexNoise(seed, complex_noise);

  // Shape, scale and add comfort noise.
  for (int i = 1; i < PART_LEN1; ++i) {
    float noise_scaling =
        sqrtf(WEBRTC_SPL_MAX(1 - suppressor_gain[i] * suppressor_gain[i], 0)) *
        sqrtf(noise_spectrum[i]);
    e_fft[0][i] += noise_scaling * complex_noise[0][i];
    e_fft[1][i] += noise_scaling * complex_noise[1][i];
  }

  // Form comfort noise for higher frequencies.
  if (generate_high_frequency_noise) {
    // Compute average noise power and nlp gain over the second half of freq
    // spectrum (i.e., 4->8khz).
    int start_avg_band = PART_LEN1 / 2;
    float upper_bands_noise_power = 0.f;
    float upper_bands_suppressor_gain = 0.f;
    for (int i = start_avg_band; i < PART_LEN1; ++i) {
      upper_bands_noise_power += sqrtf(noise_spectrum[i]);
      upper_bands_suppressor_gain +=
          sqrtf(WEBRTC_SPL_MAX(1 - suppressor_gain[i] * suppressor_gain[i], 0));
    }
    upper_bands_noise_power /= (PART_LEN1 - start_avg_band);
    upper_bands_suppressor_gain /= (PART_LEN1 - start_avg_band);

    // Shape, scale and add comfort noise.
    float noise_scaling = upper_bands_suppressor_gain * upper_bands_noise_power;
    high_frequency_comfort_noise[0][0] = 0;
    high_frequency_comfort_noise[1][0] = 0;
    for (int i = 1; i < PART_LEN1; ++i) {
      high_frequency_comfort_noise[0][i] = noise_scaling * complex_noise[0][i];
      high_frequency_comfort_noise[1][i] = noise_scaling * complex_noise[1][i];
    }
    high_frequency_comfort_noise[1][PART_LEN] = 0;
  } else {
    memset(high_frequency_comfort_noise, 0,
           2 * PART_LEN1 * sizeof(high_frequency_comfort_noise[0][0]));
  }
}
static void InitLevel(PowerLevel* level) {
  const float kBigFloat = 1E17f;
  level->averagelevel.Reset();
  level->framelevel.Reset();
  level->minlevel = kBigFloat;
}

static void InitStats(Stats* stats) {
  stats->instant = kOffsetLevel;
  stats->average = kOffsetLevel;
  stats->max = kOffsetLevel;
  stats->min = kOffsetLevel * (-1);
  stats->sum = 0;
  stats->hisum = 0;
  stats->himean = kOffsetLevel;
  stats->counter = 0;
  stats->hicounter = 0;
}

static void InitMetrics(AecCore* self) {
  self->stateCounter = 0;
  InitLevel(&self->farlevel);
  InitLevel(&self->nearlevel);
  InitLevel(&self->linoutlevel);
  InitLevel(&self->nlpoutlevel);

  InitStats(&self->erl);
  InitStats(&self->erle);
  InitStats(&self->aNlp);
  InitStats(&self->rerl);

  self->divergent_filter_fraction.Reset();
}

float CalculatePower(const int32_t* in, size_t num_samples) {
  size_t k;
  float energy = 0.0f;
//#ifdef HIFI
//  energy = vec_powerf(in,num_samples);
//#else
  for (k = 0; k < num_samples; ++k) {
    energy += in[k] * in[k];
  }
//#endif
  return energy / num_samples;
}

static void UpdateLevel(PowerLevel* level, float power) {
  level->framelevel.AddValue(power);
  if (level->framelevel.EndOfBlock()) {
    const float new_frame_level = level->framelevel.GetLatestMean();
    if (new_frame_level > 0) {
      if (new_frame_level < level->minlevel) {
        level->minlevel = new_frame_level;  // New minimum.
      } else {
        level->minlevel *= (1 + 0.001f);  // Small increase.
      }
    }
    level->averagelevel.AddValue(new_frame_level);
  }
}

static void UpdateMetrics(AecCore* aec) {
  const float actThresholdNoisy = 8.0f;
  const float actThresholdClean = 40.0f;

  const float noisyPower = 300000.0f;

  float actThreshold;

  if (aec->echoState) {  // Check if echo is likely present
    aec->stateCounter++;
  }

  if (aec->linoutlevel.framelevel.EndOfBlock()) {
    aec->divergent_filter_fraction.AddObservation(
        aec->nearlevel, aec->linoutlevel, aec->nlpoutlevel);
  }

  if (aec->farlevel.averagelevel.EndOfBlock()) {
    if (aec->farlevel.minlevel < noisyPower) {
      actThreshold = actThresholdClean;
    } else {
      actThreshold = actThresholdNoisy;
    }

    const float far_average_level = aec->farlevel.averagelevel.GetLatestMean();

    // The last condition is to let estimation be made in active far-end
    // segments only.
    if ((aec->stateCounter > (0.5f * kCountLen * kSubCountLen)) &&
        (aec->farlevel.framelevel.EndOfBlock()) &&
        (far_average_level > (actThreshold * aec->farlevel.minlevel))) {
      // ERL: error return loss.
      const float near_average_level =
          aec->nearlevel.averagelevel.GetLatestMean();
      UpdateLogRatioMetric(&aec->erl, far_average_level, near_average_level);

      // A_NLP: error return loss enhanced before the nonlinear suppression.
      const float linout_average_level =
          aec->linoutlevel.averagelevel.GetLatestMean();
      UpdateLogRatioMetric(&aec->aNlp, near_average_level,
                           linout_average_level);

      // ERLE: error return loss enhanced.
      const float nlpout_average_level =
          aec->nlpoutlevel.averagelevel.GetLatestMean();
      UpdateLogRatioMetric(&aec->erle, near_average_level,
                           nlpout_average_level);
    }

    aec->stateCounter = 0;
  }
}

static void UpdateDelayMetrics(AecCore* self) {
  int i = 0;
  int delay_values = 0;
  int median = 0;
  int lookahead = WebRtc_lookahead(self->delay_estimator);
  const int kMsPerBlock = PART_LEN / (self->mult * 8);
  int64_t l1_norm = 0;

  if (self->num_delay_values == 0) {
    // We have no new delay value data. Even though -1 is a valid |median| in
    // the sense that we allow negative values, it will practically never be
    // used since multiples of |kMsPerBlock| will always be returned.
    // We therefore use -1 to indicate in the logs that the delay estimator was
    // not able to estimate the delay.
    self->delay_median = -1;
    self->delay_std = -1;
    self->fraction_poor_delays = -1;
    return;
  }

  // Start value for median count down.
  delay_values = self->num_delay_values >> 1;
  // Get median of delay values since last update.
  for (i = 0; i < kHistorySizeBlocks; i++) {
    delay_values -= self->delay_histogram[i];
    if (delay_values < 0) {
      median = i;
      break;
    }
  }
  // Account for lookahead.
  self->delay_median = (median - lookahead) * kMsPerBlock;

  // Calculate the L1 norm, with median value as central moment.
  for (i = 0; i < kHistorySizeBlocks; i++) {
    l1_norm += abs(i - median) * self->delay_histogram[i];
  }
  self->delay_std = static_cast<int>((l1_norm + self->num_delay_values / 2) /
                                     self->num_delay_values) *
                    kMsPerBlock;

  // Determine fraction of delays that are out of bounds, that is, either
  // negative (anti-causal system) or larger than the AEC filter length.
  {
    int num_delays_out_of_bounds = self->num_delay_values;
    const int histogram_length =
        sizeof(self->delay_histogram) / sizeof(self->delay_histogram[0]);
    for (i = lookahead; i < lookahead + self->num_partitions; ++i) {
      if (i < histogram_length)
        num_delays_out_of_bounds -= self->delay_histogram[i];
    }
    self->fraction_poor_delays =
        static_cast<float>(num_delays_out_of_bounds) / self->num_delay_values;
  }

  // Reset histogram.
  memset(self->delay_histogram, 0, sizeof(self->delay_histogram));
  self->num_delay_values = 0;
}
#endif

#if 0
static int SignalBasedDelayCorrection(AecCore* self) {
  int delay_correction = 0;
  int last_delay = -2;
  assert(self);
#if !defined(WEBRTC_ANDROID)
  // On desktops, turn on correction after |kDelayCorrectionStart| frames.  This
  // is to let the delay estimation get a chance to converge.  Also, if the
  // playout audio volume is low (or even muted) the delay estimation can return
  // a very large delay, which will break the AEC if it is applied.
  if (self->frame_count < kDelayCorrectionStart) {
    //self->data_dumper->DumpRaw("aec_da_reported_delay", 1, &last_delay);
    return 0;
  }
#endif

  // 1. Check for non-negative delay estimate.  Note that the estimates we get
  //    from the delay estimation are not compensated for lookahead.  Hence, a
  //    negative |last_delay| is an invalid one.
  // 2. Verify that there is a delay change.  In addition, only allow a change
  //    if the delay is outside a certain region taking the AEC filter length
  //    into account.
  // TODO(bjornv): Investigate if we can remove the non-zero delay change check.
  // 3. Only allow delay correction if the delay estimation quality exceeds
  //    |delay_quality_threshold|.
  // 4. Finally, verify that the proposed |delay_correction| is feasible by
  //    comparing with the size of the far-end buffer.
  last_delay = WebRtc_last_delay(self->delay_estimator);
  //self->data_dumper->DumpRaw("aec_da_reported_delay", 1, &last_delay);
  if ((last_delay >= 0) && (last_delay != self->previous_delay) &&
      (WebRtc_last_delay_quality(self->delay_estimator) >
       self->delay_quality_threshold)) {
    int delay = last_delay - WebRtc_lookahead(self->delay_estimator);
    // Allow for a slack in the actual delay, defined by a |lower_bound| and an
    // |upper_bound|.  The adaptive echo cancellation filter is currently
    // |num_partitions| (of 64 samples) long.  If the delay estimate is negative
    // or at least 3/4 of the filter length we open up for correction.
    const int lower_bound = 0;
    const int upper_bound = self->num_partitions * 3 / 4;
    const int do_correction = delay <= lower_bound || delay > upper_bound;
    if (do_correction == 1) {
      int available_read = self->farend_block_buffer_.Size();
      // With |shift_offset| we gradually rely on the delay estimates.  For
      // positive delays we reduce the correction by |shift_offset| to lower the
      // risk of pushing the AEC into a non causal state.  For negative delays
      // we rely on the values up to a rounding error, hence compensate by 1
      // element to make sure to push the delay into the causal region.
      delay_correction = -delay;
      delay_correction += delay > self->shift_offset ? self->shift_offset : 1;
      self->shift_offset--;
      self->shift_offset = (self->shift_offset <= 1 ? 1 : self->shift_offset);
      if (delay_correction > available_read - self->mult - 1) {
        // There is not enough data in the buffer to perform this shift.  Hence,
        // we do not rely on the delay estimate and do nothing.
        delay_correction = 0;
      } else {
        self->previous_delay = last_delay;
        ++self->delay_correction_count;
      }
    }
  }
  // Update the |delay_quality_threshold| once we have our first delay
  // correction.
  if (self->delay_correction_count > 0) {
    float delay_quality = WebRtc_last_delay_quality(self->delay_estimator);
    delay_quality =
        (delay_quality > kDelayQualityThresholdMax ? kDelayQualityThresholdMax
                                                   : delay_quality);
    self->delay_quality_threshold =
        (delay_quality > self->delay_quality_threshold
             ? delay_quality
             : self->delay_quality_threshold);
  }
  //self->data_dumper->DumpRaw("aec_da_delay_correction", 1, &delay_correction);

  return delay_correction;
}

static void RegressorPower(
    int num_partitions,
    int latest_added_partition,
    int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],
    float x_pow[PART_LEN1]) {
  assert(latest_added_partition < num_partitions);
  memset(x_pow, 0, PART_LEN1 * sizeof(x_pow[0]));

  int partition = latest_added_partition;

  int x_fft_buf_position = partition * PART_LEN3;
  for (int i = 0; i < num_partitions; ++i) {
      for (int bin = 0; bin < PART_LEN1; ++bin) {
        //x_pow[bin] += normf(x_fft_buf_hifi[x_fft_buf_position]);
    	  float re = x_fft_buf_hifi[x_fft_buf_position];
    	  float im = x_fft_buf_hifi[x_fft_buf_position + 1];
    	  x_pow[bin] += re * re + im * im;
       // ++x_fft_buf_position;
    	  x_fft_buf_position = x_fft_buf_position + 2;
      }

      ++partition;
      if (partition == num_partitions) {
        partition = 0;
        assert((num_partitions * PART_LEN3) == x_fft_buf_position);
        x_fft_buf_position = 0;
      }
    }

}
#endif
int FilterAdaptation_n = 0;
void EchoSubtraction(
    int num_partitions,
    int extended_filter_enabled,
    int* extreme_filter_divergence,
    int filter_step_size,
    int error_threshold,
    int32_t* x_fft_hifi,//Q0
    int* x_fft_buf_block_pos,
    int32_t x_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],
    int32_t* const y,//aecq0
    long long x_pow[PART_LEN1],//Q0
    int32_t h_fft_buf_hifi[kExtendedNumPartitions * PART_LEN3],//aecq
    int32_t echo_subtractor_output[PART_LEN]) {
    
	  int32_t s_fft_hifi[PART_LEN3*2] __attribute__((aligned(16)));
	  int32_t e_extended[PART_LEN2*2];
	  int32_t s_extended[PART_LEN2*2];
	  int32_t e[PART_LEN];
	  //float s_hifi[PART_LEN];
	  int32_t e_fft_hifi[PART_LEN3*2];//aecnearq
	  long long e_fft_hifi_err[PART_LEN3*2];
	  int i;
	  int32_t* s;//, *e;

//int32_t h_tmp[128];
  // Update the x_fft_buf block position.
  (*x_fft_buf_block_pos)--;
  if ((*x_fft_buf_block_pos) == -1) {
    *x_fft_buf_block_pos = num_partitions - 1;
  }
  // Buffer x_fft.
    memcpy(x_fft_buf_hifi + (*x_fft_buf_block_pos) * PART_LEN3, x_fft_hifi,
           sizeof(int32_t) * PART_LEN3);
    memset(s_fft_hifi, 0, sizeof(s_fft_hifi));


  // Conditionally reset the echo subtraction filter if the filter has diverged
  // significantly.
  if (!extended_filter_enabled && *extreme_filter_divergence) {
     memset(h_fft_buf_hifi, 0,
             kExtendedNumPartitions * PART_LEN3 * sizeof(int32_t));
     *extreme_filter_divergence = 0;
   }
   //t0 = osKernelGetTickCount();
   //for (int n = 0; n < 4*32; n++)
   {
  // Produce echo estimate s_fft.
 //memcpy(h_tmp,h_fft_buf_hifi,128*4);//////////for test
 FilterFar(num_partitions, *x_fft_buf_block_pos, x_fft_buf_hifi, h_fft_buf_hifi, s_fft_hifi);//0.33ms@M3
 }
 //iet_printf_msg("FilterFar time =%d \r\n", osKernelGetTickCount()-t0);
  // Compute the time-domain echo estimate s.
    //#ifdef HARD_FFT
    //hard_rfft_int32(1, s_fft_hifi, s_extended);//hyli
    //#else
    //realFFT_Inverse(s_extended, s_fft_hifi, 128 );
    arm_rfft_q31(&S_ifft128, s_fft_hifi, s_extended);//hyli
    //#endif
    
    s = &s_extended[PART_LEN];//Q0
    //e = &e_extended[PART_LEN];//
    // Compute the time-domain echo prediction error.
    for (i = 0; i < PART_LEN; ++i){
        e[i] = y[i] - s[i]; //y is nearend_block
    }

  // Compute the frequency domain echo prediction error.
  memset(e_extended, 0, sizeof(int32_t) * PART_LEN);
  memcpy(e_extended + PART_LEN, e, sizeof(int32_t) * PART_LEN);
  //#ifdef HARD_FFT
  //  hard_rfft_int32(0, e_extended, e_fft_hifi);//hyli
  //  #else
  //realFFT_forward(e_fft_hifi, e_extended, 128);
  arm_rfft_q31(&S_fft128, e_extended, e_fft_hifi);//hyli
  //#endif
   //t0 = osKernelGetTickCount();
   //for (int n = 0; n < 4*32; n++)
   {
    ScaleErrorSignal(filter_step_size, error_threshold, x_pow, e_fft_hifi, e_fft_hifi_err);//0.48ms@M3
    }
//iet_printf_msg("ScaleErrorSignal time =%d \r\n", osKernelGetTickCount()-t0);
 //start0 = SysTickCnt_test;
 //t0 = osKernelGetTickCount();
 //FilterAdaptation_n++;
 //for (int n = 0; n < 4*32; n++)
 {

 #ifdef HARD_FFT
    FilterAdaptation_comp(num_partitions, *x_fft_buf_block_pos,
                               x_fft_buf_hifi, e_fft_hifi, h_fft_buf_hifi);//3.14ms@M3
     #if 0
     if(FilterAdaptation_n<8)//8
     {
        FilterAdaptation_n++;
        for(int i = 0; i < PART_LEN3; i++) {
            iet_printf_msg("%d, ",e_fft_hifi[i]);
            if((i+1)%32==0)
            {
                iet_printf_msg("\r\n");
            }
        }
        iet_printf_msg("\r\n");
        for(int i = 0; i < kExtendedNumPartitions * PART_LEN3; i++) {
            iet_printf_msg("%d, ",h_fft_buf_hifi[i]);
            if((i+1)%32==0)
            {
                iet_printf_msg("\r\n");
            }
        }
        iet_printf_msg("\r\n");
     }
     #endif
 #else
 FilterAdaptation(num_partitions, *x_fft_buf_block_pos,
                               x_fft_buf_hifi, e_fft_hifi_err, h_fft_buf_hifi);
 #endif 

 
 }
//iet_printf_msg("FilterAdaptation time =%d \r\n", osKernelGetTickCount()-t0);
    memcpy(echo_subtractor_output, e, sizeof(int32_t) * PART_LEN);

}
#if 1
int fastPow3(int a, int b) {//Q16
    int result = 65536;
    while(b > 0) {
        if ((b & 1) != 0) {
            result = ((long long)result*a)>>16;
            b--;
        }
        a = ((long long)a*a)>>16;
        b >>= 1;
    }

    return result;
}
static void Overdrive(int overdrive_scaling,//Q9
                      const int hNlFb,//Q16
                      uint16_t hNl[PART_LEN1]) {//Q16
  for (int i = 0; i < PART_LEN1; ++i) {
    // Weight subbands
    if (hNl[i] > hNlFb) {
      hNl[i] = ((long long)WebRtcAec_weightCurve[i] * hNlFb +(long long)(65536 - WebRtcAec_weightCurve[i]) * hNl[i])>>16;
    }
    int scale = (overdrive_scaling * WebRtcAec_overDriveCurve[i]+131072)>>18;
    hNl[i] = fastPow3(hNl[i], scale);
    //hNl[i] = hNl[i]*hNl[i]*hNl[i]*hNl[i];
  }
}
static int NLP_enable=1;
void Suppress(const uint16_t hNl[PART_LEN1],//Q16
		int32_t efw_hifi[PART_LEN3]

		             ) {
  if(NLP_enable){  	             
      for (int i = 0; i < PART_LEN1; ++i) {
    	  // Suppress error signal
    	      efw_hifi[2*i] = ((long long)efw_hifi[2*i]*hNl[i])>>16;
    	      efw_hifi[2*i + 1] = (-(long long)efw_hifi[2*i + 1]*hNl[i])>>16;//取共轭给后续的ifft

    	      // Ooura fft returns incorrect sign on imaginary component. It matters here
    	      // because we are making an additive change with comfort noise.
    	     // efw_hifi[2*i + 1] *= -1; for next scaled inverse fft take the -1 sign
      }
  }
  else{
        for (int i = 0; i < PART_LEN1; ++i) {
    	  // Suppress error signal
    	      efw_hifi[2*i] = efw_hifi[2*i];
    	      efw_hifi[2*i + 1] = -efw_hifi[2*i + 1];//取共轭给后续的ifft
      }
  }
}
// Threshold to protect against the ill-effects of a zero far-end.
const int WebRtcAec_kMinFarendPSD = 15;

// Updates the following smoothed Power Spectral Densities (PSD):
//  - sd  : near-end
//  - se  : residual echo
//  - sx  : far-end
//  - sde : cross-PSD of near-end and residual echo
//  - sxd : cross-PSD of near-end and far-end
//
// In addition to updating the PSDs, also the filter diverge state is
// determined.
#if 0
static void UpdateCoherenceSpectra(int mult,
                                   bool extended_filter_enabled,
                                   int32_t efw_hifi[PART_LEN3],
								   int32_t dfw_hifi[PART_LEN3],
								   int32_t xfw_hifi[PART_LEN3],
                                   CoherenceState* coherence_state,
                                   short* filter_divergence_state,
                                   int* extreme_filter_divergence) {
  // Power estimate smoothing coefficients.
  const long long* ptrGCoh =
      extended_filter_enabled
          ? WebRtcAec_kExtendedSmoothingCoefficients[mult - 1]
          : WebRtcAec_kNormalSmoothingCoefficients[mult - 1];
  int i;
  long long sdSum = 0, seSum = 0, ltmp;
  int COQ = 0;
  int GQ = 13;
  for (i = 0; i < PART_LEN1; i++) {
     coherence_state->sd[i] = (((long long)coherence_state->sd[i]*ptrGCoh[0])>>GQ) + ((((long long)dfw_hifi[2*i] * dfw_hifi[2*i] + (long long)dfw_hifi[2*i + 1] * dfw_hifi[2*i + 1])*ptrGCoh[1])>>(GQ+COQ));//Q(-COQ)
     coherence_state->se[i] = (((long long)coherence_state->se[i]*ptrGCoh[0])>>GQ) + ((((long long)efw_hifi[2*i] * efw_hifi[2*i] + (long long)efw_hifi[2*i + 1] * efw_hifi[2*i + 1])*ptrGCoh[1])>>(GQ+COQ));
     // We threshold here to protect against the ill-effects of a zero farend.
     // The threshold is not arbitrarily chosen, but balances protection and
     // adverse interaction with the algorithm's tuning.
     // TODO(bjornv): investigate further why this is so sensitive.
     coherence_state->sx[i] = (((long long)coherence_state->sx[i]*ptrGCoh[0])>>GQ) + 
     ((WEBRTC_SPL_MAX((long long)xfw_hifi[2*i] * xfw_hifi[2*i] + (long long)xfw_hifi[2*i + 1] * xfw_hifi[2*i + 1],WebRtcAec_kMinFarendPSD)*ptrGCoh[1])>>(GQ+COQ));

     coherence_state->sde[i][0] = (((long long)coherence_state->sde[i][0]*ptrGCoh[0])>>GQ) + ((((long long)dfw_hifi[2*i] * efw_hifi[2*i] + (long long)dfw_hifi[2*i + 1] * efw_hifi[2*i + 1])*ptrGCoh[1])>>(GQ+COQ));//Q(-COQ)
     coherence_state->sde[i][1] = (((long long)coherence_state->sde[i][1]*ptrGCoh[0])>>GQ) + ((((long long)dfw_hifi[2*i] * efw_hifi[2*i + 1] - (long long)dfw_hifi[2*i + 1] * efw_hifi[2*i])*ptrGCoh[1])>>(GQ+COQ));

     coherence_state->sxd[i][0] = (((long long)coherence_state->sxd[i][0]*ptrGCoh[0])>>GQ) + ((((long long)dfw_hifi[2*i] * xfw_hifi[2*i] + (long long)dfw_hifi[2*i + 1] * xfw_hifi[2*i + 1])*ptrGCoh[1])>>(GQ+COQ));
     coherence_state->sxd[i][1] = (((long long)coherence_state->sxd[i][1]*ptrGCoh[0])>>GQ) + ((((long long)dfw_hifi[2*i] * xfw_hifi[2*i + 1] - (long long)dfw_hifi[2*i + 1] * xfw_hifi[2*i])*ptrGCoh[1])>>(GQ+COQ));

     sdSum += coherence_state->sd[i];
     seSum += coherence_state->se[i];
   }

  // Divergent filter safeguard update.
  *filter_divergence_state =
      (*filter_divergence_state ? 269 : 256) * seSum > (sdSum<<8);//1.05f : 1.0f

  // Signal extreme filter divergence if the error is significantly larger
  // than the nearend (13 dB).
  *extreme_filter_divergence = ((seSum<<8) > (5107 * sdSum));//19.95f
}
#endif
#if 1
static void UpdateCoherenceSpectra(int mult,
                                   bool extended_filter_enabled,
                                   int32_t efw[PART_LEN3],
								   int32_t dfw[PART_LEN3],
								   int32_t xfw[PART_LEN3],
                                   CoherenceState* coherence_state,
                                   short* filter_divergence_state,
                                   int* extreme_filter_divergence) {
  // Power estimate smoothing coefficients.
  const long long* ptrGCoh =
      extended_filter_enabled
          ? WebRtcAec_kExtendedSmoothingCoefficients[mult - 1]
          : WebRtcAec_kNormalSmoothingCoefficients[mult - 1];
  int i;
  long long sdSum = 0, seSum = 0, ltmp;
  int COQ = 0;
  int GQ = 13;
  long long tmpe = ptrGCoh[1], tmpse, tmpd, tmpsd;
  for (i = 0; i < PART_LEN1; i++) {
    tmpd = (ptrGCoh[0] * coherence_state->sd[i]);
    tmpsd = (((long long)dfw[2*i] * dfw[2*i]) + ((long long)dfw[2*i+1] * dfw[2*i+1]));
    tmpsd = ptrGCoh[1] *tmpsd;
    coherence_state->sd[i] =
        ((ptrGCoh[0] * coherence_state->sd[i]) +
        (ptrGCoh[1] * (((long long)dfw[2*i] * dfw[2*i]) + ((long long)dfw[2*i+1] * dfw[2*i+1]))))>>GQ;
    coherence_state->se[i] =
        ((ptrGCoh[0] * coherence_state->se[i]) +
        (ptrGCoh[1] * (((long long)efw[2*i] * efw[2*i]) + ((long long)efw[2*i+1] * efw[2*i+1]))))>>GQ;


    coherence_state->sx[i] =
        ((ptrGCoh[0] * coherence_state->sx[i]) +
        (ptrGCoh[1] *
            WEBRTC_SPL_MAX(((long long)xfw[2*i] * xfw[2*i]) + ((long long)xfw[2*i+1] * xfw[2*i+1]),
                           WebRtcAec_kMinFarendPSD)))>>GQ;

    coherence_state->sde[i][0] =
        ((ptrGCoh[0] * coherence_state->sde[i][0]) +
        (ptrGCoh[1] * (((long long)dfw[2*i] * efw[2*i]) + ((long long)dfw[2*i+1] * efw[2*i+1]))))>>GQ;
    coherence_state->sde[i][1] =
        ((ptrGCoh[0] * coherence_state->sde[i][1]) +
        (ptrGCoh[1] * (((long long)dfw[2*i] * efw[2*i+1]) - ((long long)dfw[2*i+1] * efw[2*i]))))>>GQ;

    coherence_state->sxd[i][0] =
        ((ptrGCoh[0] * coherence_state->sxd[i][0]) +
        (ptrGCoh[1] * (((long long)dfw[2*i] * xfw[2*i]) + ((long long)dfw[2*i+1] * xfw[2*i+1]))))>>GQ;
    coherence_state->sxd[i][1] =
        ((ptrGCoh[0] * coherence_state->sxd[i][1]) +
        (ptrGCoh[1] * (((long long)dfw[2*i] * xfw[2*i+1]) - ((long long)dfw[2*i+1] * xfw[2*i]))))>>GQ;

    sdSum += coherence_state->sd[i];
    seSum += coherence_state->se[i];
  }

  // Divergent filter safeguard update.
  *filter_divergence_state =
      (*filter_divergence_state ? 269 : 256) * seSum > (sdSum<<8);//1.05f : 1.0f

  // Signal extreme filter divergence if the error is significantly larger
  // than the nearend (13 dB).
  *extreme_filter_divergence = ((seSum<<8) > (5107 * sdSum));//19.95f
}
#else
static void UpdateCoherenceSpectra(int mult,
                                   bool extended_filter_enabled,
                                   int32_t efw[PART_LEN3],
								   int32_t dfw[PART_LEN3],
								   int32_t xfw[PART_LEN3],
                                   CoherenceState* coherence_state,
                                   short* filter_divergence_state,
                                   int* extreme_filter_divergence) 
{
  // Power estimate smoothing coefficients.
  const float* ptrGCoh =
      extended_filter_enabled
          ? fWebRtcAec_kExtendedSmoothingCoefficients[mult - 1]
          : fWebRtcAec_kNormalSmoothingCoefficients[mult - 1];
  int i;
  float sdSum = 0, seSum = 0;
long long tmpe, tmpse;
  for (i = 0; i < PART_LEN1; i++) {
    coherence_state->sd[i] =
        ptrGCoh[0] * coherence_state->sd[i] +
        ptrGCoh[1] * (((long long)dfw[2*i] * dfw[2*i]) + ((long long)dfw[2*i+1] * dfw[2*i+1]));
    coherence_state->se[i] =
        ptrGCoh[0] * coherence_state->se[i] +
        ptrGCoh[1] * (((long long)efw[2*i] * efw[2*i]) + ((long long)efw[2*i+1] * efw[2*i+1]));
        tmpe = ((long long)efw[2*i] * efw[2*i]) + ((long long)efw[2*i+1] * efw[2*i+1]);
        tmpe = ptrGCoh[1] * tmpe;
        tmpse = ptrGCoh[0] * coherence_state->se[i];
    // We threshold here to protect against the ill-effects of a zero farend.
    // The threshold is not arbitrarily chosen, but balances protection and
    // adverse interaction with the algorithm's tuning.
    // TODO(bjornv): investigate further why this is so sensitive.
    coherence_state->sx[i] =
        ptrGCoh[0] * coherence_state->sx[i] +
        ptrGCoh[1] *
            WEBRTC_SPL_MAX(((long long)xfw[2*i] * xfw[2*i]) + ((long long)xfw[2*i+1] * xfw[2*i+1]),
                           WebRtcAec_kMinFarendPSD);

    coherence_state->sde[i][0] =
        ptrGCoh[0] * coherence_state->sde[i][0] +
        ptrGCoh[1] * (((long long)dfw[2*i] * efw[2*i]) + ((long long)dfw[2*i+1] * efw[2*i+1]));
    coherence_state->sde[i][1] =
        ptrGCoh[0] * coherence_state->sde[i][1] +
        ptrGCoh[1] * (((long long)dfw[2*i] * efw[2*i+1]) - ((long long)dfw[2*i+1] * efw[2*i]));

    coherence_state->sxd[i][0] =
        ptrGCoh[0] * coherence_state->sxd[i][0] +
        ptrGCoh[1] * (((long long)dfw[2*i] * xfw[2*i]) + ((long long)dfw[2*i+1] * xfw[2*i+1]));
    coherence_state->sxd[i][1] =
        ptrGCoh[0] * coherence_state->sxd[i][1] +
        ptrGCoh[1] * (((long long)dfw[2*i] * xfw[2*i+1]) - ((long long)dfw[2*i+1] * xfw[2*i]));

    sdSum += coherence_state->sd[i];
    seSum += coherence_state->se[i];
  }

  // Divergent filter safeguard update.
  *filter_divergence_state =
      (*filter_divergence_state ? 1.05f : 1.0f) * seSum > sdSum;

  // Signal extreme filter divergence if the error is significantly larger
  // than the nearend (13 dB).
  *extreme_filter_divergence = (seSum > (19.95f * sdSum));
}

#endif
// Window time domain data to be used by the fft.
 void WindowData(int32_t* x_windowed, int32_t* x) {
  int i;
  for (i = 0; i < PART_LEN; i++) {
    x_windowed[i]            = ((long long)x[i] * WebRtcAec_sqrtHanning[i])>>15;
    x_windowed[PART_LEN + i] = ((long long)x[PART_LEN + i] * WebRtcAec_sqrtHanning[PART_LEN - i])>>15;
  }
}

static void FormSuppressionGain(AecCore* aec,
                                uint16_t cohde[PART_LEN1],//Q16
                                uint16_t cohxd[PART_LEN1],//Q16
                                uint16_t hNl[PART_LEN1]) {//Q16
  int hNlDeAvg, hNlXdAvg;
  uint16_t hNlPref[kPrefBandSize];
  int hNlFb = 0, hNlFbLow = 0;//Q16
  const int prefBandSize = kPrefBandSize / aec->mult;
  //const float prefBandQuant = 0.75f, prefBandQuantLow = 0.5f;
  const int minPrefBand = 4 / aec->mult;
  // Power estimate smoothing coefficients.
  const int* min_overdrive = aec->extended_filter_enabled
                                   ? kExtendedMinOverDrive
                                   : kNormalMinOverDrive;

  hNlXdAvg = 0;
  for (int i = minPrefBand; i < prefBandSize + minPrefBand; ++i) {
    hNlXdAvg += cohxd[i];
  }
  hNlXdAvg /= prefBandSize;
  hNlXdAvg = 65536 - hNlXdAvg;

  hNlDeAvg = 0;
  for (int i = minPrefBand; i < prefBandSize + minPrefBand; ++i) {
    hNlDeAvg += cohde[i];
  }
  hNlDeAvg /= prefBandSize;

  if (hNlXdAvg < 49152 && hNlXdAvg < aec->hNlXdAvgMin) {//0.75
    aec->hNlXdAvgMin = hNlXdAvg;//Q16
  }

  if (hNlDeAvg > 64225 && hNlXdAvg > 58982) {//0.98 0.9
    aec->stNearState = 1;
  } else if (hNlDeAvg < 62259 || hNlXdAvg < 52429) {//0.95 0.8
    aec->stNearState = 0;
  }

  if (aec->hNlXdAvgMin == 65536) {
    aec->echoState = 0;
    aec->overDrive = min_overdrive[aec->nlp_mode];

    if (aec->stNearState == 1) {
      memcpy(hNl, cohde, sizeof(hNl[0]) * PART_LEN1);
      hNlFb = hNlDeAvg;
      hNlFbLow = hNlDeAvg;
    } else {
      for (int i = 0; i < PART_LEN1; ++i) {
        hNl[i] = 65536 - cohxd[i];
        hNl[i] = MAX(hNl[i], 0);
      }
      hNlFb = hNlXdAvg;
      hNlFbLow = hNlXdAvg;
    }
  } else {
    if (aec->stNearState == 1) {
      aec->echoState = 0;
      memcpy(hNl, cohde, sizeof(hNl[0]) * PART_LEN1);
      hNlFb = hNlDeAvg;
      hNlFbLow = hNlDeAvg;
    } else {
      aec->echoState = 1;
      for (int i = 0; i < PART_LEN1; ++i) {
        hNl[i] = WEBRTC_SPL_MIN(cohde[i], 65536 - cohxd[i]);
        hNl[i] = MAX(hNl[i], 0);
      }

      // Select an order statistic from the preferred bands.
      // TODO(peah): Using quicksort now, but a selection algorithm may be
      // preferred.
      memcpy(hNlPref, &hNl[minPrefBand], sizeof(hNl[0]) * prefBandSize);
      qsort(hNlPref, prefBandSize, sizeof(hNl[0]), CmpFloat);
      hNlFb =
          hNlPref[(((prefBandSize - 1)*3)>>2)];
      hNlFbLow = hNlPref[((prefBandSize - 1)>>1)];
    }
  }

  // Track the local filter minimum to determine suppression overdrive.
  if (hNlFbLow < 39322 && hNlFbLow < aec->hNlFbLocalMin) {//Q16 0.6
    aec->hNlFbLocalMin = hNlFbLow;
    aec->hNlFbMin = hNlFbLow;
    aec->hNlNewMin = 1;
    aec->hNlMinCtr = 0;
  }
  aec->hNlFbLocalMin =
      WEBRTC_SPL_MIN(aec->hNlFbLocalMin + 52 / aec->mult, 65536);//0.0008
  aec->hNlXdAvgMin = WEBRTC_SPL_MIN(aec->hNlXdAvgMin + 39 / aec->mult, 65536);//0.0006

  if (aec->hNlNewMin == 1) {
    aec->hNlMinCtr++;
  }
  if (aec->hNlMinCtr == 2) {
    aec->hNlNewMin = 0;
    aec->hNlMinCtr = 0;
    aec->overDrive = WEBRTC_SPL_MAX(
        kTargetSupp[aec->nlp_mode] /(log_16q16(aec->hNlFbMin + 1) + 1), min_overdrive[aec->nlp_mode]);
    //aec->overDrive = WEBRTC_SPL_MAX(
    //    (kTargetSupp1[aec->nlp_mode]*512 /static_cast<float>(log(aec->hNlFbMin/65536.0 + 1e-10f) + 1e-10f)), min_overdrive[aec->nlp_mode]);
  }

  // Smooth the overdrive.
  if ((aec->overDrive) < aec->overdrive_scaling) {
    aec->overdrive_scaling =
        ((507 * aec->overdrive_scaling + 5 * aec->overDrive)>>9);//0.99f 0.01
  } else {
    aec->overdrive_scaling =
        ((461 * aec->overdrive_scaling + 51 * aec->overDrive)>>9);//0.9 0.1
  }

  // Apply the overdrive.
  Overdrive(aec->overdrive_scaling, hNlFb, hNl);
}

void EchoSuppression(
                            AecCore* aec,
                            int32_t* nearend_extended_block_lowest_band,//aecq0
                            int32_t farend_extended_block[PART_LEN2],//aecq0
                            int32_t* echo_subtractor_output,//aecsq=aecq0
                            int32_t output[NUM_HIGH_BANDS_MAX + 1][PART_LEN]) {//Q0

	  int32_t efw_hifi[PART_LEN3*2];
	  int32_t xfw_hifi[PART_LEN3*2];
	  int32_t dfw_hifi[PART_LEN3*2];
	  //float comfortNoiseHband_hifi[PART_LEN3];


  int32_t fft[PART_LEN2*2];
  //float nlpGainHband;
  int i;
  size_t j;

  // Coherence and non-linear filter
  uint16_t cohde[PART_LEN1], cohxd[PART_LEN1];//Q16
  uint16_t hNl[PART_LEN1];

  // Filter energy
  const int delayEstInterval = 10 * aec->mult;

  int32_t* xfw_ptr = NULL;
  // Update eBuf with echo subtractor output.
  memcpy(aec->eBuf + PART_LEN, echo_subtractor_output,
         sizeof(int32_t) * PART_LEN);
 if(NLP_enable){  
  // Analysis filter banks for the echo suppressor.
  
  // Windowed near-end ffts.
  WindowData(fft, nearend_extended_block_lowest_band);//aecq0
#ifdef HARD_FFT
        #ifndef INT_MODE         
        hard_rfft_int32(0, fft, dfw_hifi);//hyli
        #else
        hard_rfft_int32_INT(0, fft, dfw_hifi);
        #endif
    #else
  //realFFT_forward(dfw_hifi, fft, 128);
  arm_rfft_q31(&S_fft128, fft, dfw_hifi);//hyli
  #endif
  
  // Windowed echo suppressor output ffts.
  WindowData(fft, aec->eBuf);//aecq0
#ifdef HARD_FFT
        #ifndef INT_MODE         
        hard_rfft_int32(0, fft, efw_hifi);//hyli
        #else
        hard_rfft_int32_INT(0, fft, efw_hifi);
        #endif
    #else
  //realFFT_forward(efw_hifi, fft, 128);
  arm_rfft_q31(&S_fft128, fft, efw_hifi);//hyli
  #endif
  // NLP
#if 0
    int sum_zeron = 0;
    for( i = 0; i < PART_LEN2; ++i)
    {
        if(abs(farend_extended_block[i])> 50)
        {
            break;
        }
        sum_zeron++;
    }
    if(sum_zeron == 128)
    {
        memset(xfw_hifi, 0, PART_LEN3*sizeof(int32_t));
    }
    else
	#endif
    {
        // Convert far-end partition to the frequency domain with windowing.
        WindowData(fft, farend_extended_block);//aecq0  
        #ifdef HARD_FFT
            #ifndef INT_MODE         
            hard_rfft_int32(0, fft, xfw_hifi);//hyli
            #else
            hard_rfft_int32_INT(0, fft, xfw_hifi);
            #endif
        #else
        //realFFT_forward(xfw_hifi, fft, 128);
        arm_rfft_q31(&S_fft128, fft, xfw_hifi);//hyli
        #endif
    }
  
  // image part take conjunction for the result from core ooura_fft/fft_realf_ie are conjuncted
  for( i = 0; i < PART_LEN1; ++i)
    {
	  efw_hifi[2*i + 1] = -efw_hifi[2*i + 1];
	  xfw_hifi[2*i + 1] = -xfw_hifi[2*i + 1];
	  dfw_hifi[2*i + 1] = -dfw_hifi[2*i + 1];
    }   
  xfw_ptr = &xfw_hifi[0];
  // Buffer far.
  memcpy(aec->xfwBuf_hifi, xfw_ptr, sizeof(int32_t) * PART_LEN3);
  aec->delayEstCtr++;
    if (aec->delayEstCtr == delayEstInterval) {
      aec->delayEstCtr = 0;
      aec->delayIdx = PartitionDelay(aec->num_partitions, aec->wfBuf_hifi);//aecq
    }
    // Use delayed far.
     memcpy(xfw_hifi, aec->xfwBuf_hifi + aec->delayIdx * PART_LEN3,
            sizeof(xfw_hifi[0])  * PART_LEN3);
}
else
{
  int shifted = 0;
  // Windowed echo suppressor output ffts.
  WindowData(fft, aec->eBuf);//aecq0
    #ifdef HARD_FFT
        #ifndef INT_MODE         
        hard_rfft_int32(0, fft, efw_hifi);//hyli
        #else
        hard_rfft_int32_INT(0, fft, efw_hifi);
        #endif
    #else
  //realFFT_forward(efw_hifi, fft, 128);
  arm_rfft_q31(&S_fft128, fft, efw_hifi);//hyli
  #endif
  // image part take conjunction for the result from core ooura_fft/fft_realf_ie are conjuncted
  for( i = 0; i < PART_LEN1; ++i)
    {
	  efw_hifi[2*i + 1] = -efw_hifi[2*i + 1];
    }   
}            
if(NLP_enable){  
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{
  UpdateCoherenceSpectra(aec->mult, aec->extended_filter_enabled == 1,
                                      efw_hifi, dfw_hifi, xfw_hifi, &aec->coherence_state,
                                      &aec->divergeState,
                                      &aec->extreme_filter_divergence);
}                                      
//iet_printf_msg("UpdateCoherenceSpectra time =%d \r\n", osKernelGetTickCount()-t0);        
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{
  ComputeCoherence(&aec->coherence_state, cohde, cohxd);
}                                      
//iet_printf_msg("ComputeCoherence time =%d \r\n", osKernelGetTickCount()-t0);  
//printf("cycle=%d,\n",GETCLOCK()-t0);
    //t0=GETCLOCK();
  // Select the microphone signal as output if the filter is deemed to have
  // diverged.
  if (aec->divergeState) {
 memcpy(efw_hifi, dfw_hifi, sizeof(efw_hifi[0]) *  PART_LEN3);
   }
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{   
  FormSuppressionGain(aec, cohde, cohxd, hNl);
}                                      
//iet_printf_msg("FormSuppressionGain time =%d \r\n", osKernelGetTickCount()-t0);  
 }
  //aec->data_dumper->DumpRaw("aec_nlp_gain", PART_LEN1, hNl);
  Suppress(hNl, efw_hifi);
  // Add comfort noise.
  // ComfortNoise(aec->num_bands > 1, &aec->seed, efw, comfortNoiseHband,
  //              aec->noisePow, hNl);
  // Inverse error fft.
  #ifdef HARD_FFT
    #ifndef INT_MODE         
        hard_rfft_int32(1, efw_hifi, fft);//hyli
        #else
        hard_rfft_int32_INT(1, efw_hifi, fft);
        #endif
    #else
  //realFFT_Inverse(fft, efw_hifi, 128);
  arm_rfft_q31(&S_ifft128, efw_hifi, fft);//hyli
    #endif
    
  for (i = 0; i < (PART_LEN); i++) {
    output[0][i] = (fft[i] * WebRtcAec_sqrtHanning[i] + aec->outBuf[i] * WebRtcAec_sqrtHanning[64+i])>>(15);//Q0
  }

  memcpy(aec->outBuf, &fft[PART_LEN], PART_LEN * sizeof(aec->outBuf[0]));
// optimization unnecessary

   #if 0
  // For H band
  if (aec->num_bands > 1) {
    // H band gain
    // average nlp over low band: average over second half of freq spectrum
    // (4->8khz)
    GetHighbandGain(hNl, &nlpGainHband);

    ScaledInverseFft(comfortNoiseHband_hifi, fft, 0);
    // compute gain factor
    for (j = 1; j < aec->num_bands; ++j) {
      for (i = 0; i < PART_LEN; i++) {
        output[j][i] = aec->previous_nearend_block[j][i] * nlpGainHband;
      }
    }

    // Add some comfort noise where Hband is attenuated.
    for (i = 0; i < PART_LEN; i++) {
      output[1][i] += cnScaleHband * fft[i];
    }

    // Saturate output to keep it in the allowed range.
    for (j = 1; j < aec->num_bands; ++j) {
      for (i = 0; i < PART_LEN; i++) {
        output[j][i] = WEBRTC_SPL_SAT(WEBRTC_SPL_WORD16_MAX, output[j][i],
                                      WEBRTC_SPL_WORD16_MIN);
      }
    }
  }
  #endif
  // Copy the current block to the old position.
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{ 
  memcpy(aec->eBuf, aec->eBuf + PART_LEN, sizeof(int32_t) * PART_LEN);
  if(NLP_enable){  
  int wlen = kExtendedNumPartitions * PART_LEN3 - PART_LEN3;
  //memmove(aec->xfwBuf_hifi + PART_LEN3, aec->xfwBuf_hifi, wlen*sizeof(int32_t));
    while(wlen--) 
    {
      *(aec->xfwBuf_hifi + PART_LEN3 + wlen) = *(aec->xfwBuf_hifi + wlen);
    }
  //memmove(aec->xfwBuf_hifi + PART_LEN3, aec->xfwBuf_hifi,
   //         sizeof(aec->xfwBuf_hifi) - sizeof(float) * PART_LEN3);
  }
 }                                      
//iet_printf_msg("memmove time =%d \r\n", osKernelGetTickCount()-t0);   
}

#endif
void ProcessNearendBlock(
    AecCore* aec,
    int32_t farend_extended_block_lowest_band[PART_LEN2],//aecq0
    int32_t nearend_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN],//aecq0
    int32_t output_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN]) {//Q0
  //size_t i;

  int32_t fft[PART_LEN2*2]  __attribute__((aligned(16)));
  int32_t nearend_extended_block_lowest_band[PART_LEN2];

  int32_t farend_fft_hifi[PART_LEN3*2];
  int32_t nearend_fft_hifi[PART_LEN3*2];
    int32_t echo_subtractor_output[PART_LEN];

 // float Pow_Part1[PART_LEN1];
 // float Pow_Part2[PART_LEN1];
 // float far_spectrum[PART_LEN1];
 // float near_spectrum[PART_LEN1];


  long long far_spectrum = 0;
  long long near_spectrum = 0;
  //float abs_far_spectrum[PART_LEN1];
  //float abs_near_spectrum[PART_LEN1];
  const int gPow[2] = {1843, 205};//const float gPow[2] = {0.9f, 0.1f}; Q11




  // Noise estimate constants.
  const int noiseInitBlocks = 500 * aec->mult;
  const int step = 819;//0.1f; //Q13
  const int ramp = 16387;//1.0002f; //Q14
  const int sramp = 6555;//Q16 step*ramp
  const int gInitNoise[2] = {32735, 33};//{0.999f, 0.001f}; //Q15

#if 0
  if (aec->metricsMode == 1) {
    // Update power levels
    UpdateLevel(
        &aec->farlevel,
        CalculatePower(&farend_extended_block_lowest_band[PART_LEN], PART_LEN));
    UpdateLevel(&aec->nearlevel,
                CalculatePower(&nearend_block[0][0], PART_LEN));
  }
#endif
  // Convert far-end signal to the frequency domain.
  memcpy(fft, farend_extended_block_lowest_band, sizeof(int32_t) * PART_LEN2);
//#ifdef HARD_FFT
    //hard_rfft_int32(0, fft, farend_fft_hifi);//hyli
//#else
  //realFFT_forward(farend_fft_hifi, fft, 128);
  arm_rfft_q31(&S_fft128, fft, farend_fft_hifi);//hyli
//#endif
    // Power smoothing.
    int i;
    #if 0
    if (aec->refined_adaptive_filter_enabled){
	  for (i = 0; i < PART_LEN1; ++i) {
	       //far_spectrum = (float)farend_fft_hifi[2*i] * farend_fft_hifi[2*i] +
	    	//	   (float)farend_fft_hifi[2*i+1] * farend_fft_hifi[2*i+1];
	       // Calculate the magnitude spectrum.
	       //abs_far_spectrum[i] = sqrtf(far_spectrum);
	     }

       RegressorPower(aec->num_partitions, aec->xfBufBlockPos, aec->xfBuf_hifi,
                       aec->xPow);
     } else
     #endif
     {

    	 for (i = 0; i < PART_LEN1; ++i) {
    	       far_spectrum = ((long long)farend_fft_hifi[2*i] * farend_fft_hifi[2*i] +
    	    		   (long long)farend_fft_hifi[2*i+1] * farend_fft_hifi[2*i+1]);
    	       aec->xPow[i] = (aec->xPow[i]*gPow[0] + far_spectrum* aec->num_partitions*gPow[1])>>11;
    	       // Calculate the magnitude spectrum.
    	       //abs_far_spectrum[i] = sqrtf(far_spectrum);
    	     }
    }

  // Form extended nearend frame.
  memcpy(&nearend_extended_block_lowest_band[0],
         &aec->previous_nearend_block[0][0], sizeof(int32_t) * PART_LEN);
  memcpy(&nearend_extended_block_lowest_band[PART_LEN], &nearend_block[0][0],
         sizeof(int32_t) * PART_LEN);
  #if 0

  // Convert near-end signal to the frequency domain.
  memcpy(fft, nearend_extended_block_lowest_band, sizeof(int32_t) * PART_LEN2);
  realFFT_forward(nearend_fft_hifi, fft, 128);

  for (i = 0; i < PART_LEN1; ++i) {
	    near_spectrum = (long long)nearend_fft_hifi[2*i] * nearend_fft_hifi[2*i] +
	    		(long long)nearend_fft_hifi[2*i + 1] * nearend_fft_hifi[2*i + 1];
	    aec->dPow[i] = (aec->dPow[i]*gPow[0] + near_spectrum*gPow[1])>>11;
	    // Calculate the magnitude spectrum.
	    //abs_near_spectrum[i] = sqrtf(near_spectrum);
	  }

  // Estimate noise power. Wait until dPow is more stable.
  if (aec->noiseEstCtr > 50) {
    for (i = 0; i < PART_LEN1; i++) {
      if (aec->dPow[i] < aec->dMinPow[i]) {
        aec->dMinPow[i] = aec->dPow[i] + (((aec->dMinPow[i] - aec->dPow[i])*sramp)>>16);
      } else {
        aec->dMinPow[i] = (aec->dMinPow[i]*ramp)>>14;
      }
    }
  }
  // Smooth increasing noise power from zero at the start,
  // to avoid a sudden burst of comfort noise.
  if (aec->noiseEstCtr < noiseInitBlocks) {
    aec->noiseEstCtr++;
    for (i = 0; i < PART_LEN1; i++) {
      if (aec->dMinPow[i] > aec->dInitMinPow[i]) {
        aec->dInitMinPow[i] = (aec->dInitMinPow[i]*gInitNoise[0] + aec->dMinPow[i]*gInitNoise[1])>>15;
      } else {
        aec->dInitMinPow[i] = aec->dMinPow[i];
      }
    }
    //aec->noisePow = aec->dInitMinPow;
  } else {
    //aec->noisePow = aec->dMinPow;
  }
  #endif
    #if 0
  // Block wise delay estimation used for logging
  if (aec->delay_logging_enabled) {
    if (WebRtc_AddFarSpectrumFloat(aec->delay_estimator_farend,
                                   abs_far_spectrum, PART_LEN1) == 0) {
      int delay_estimate = WebRtc_DelayEstimatorProcessFloat(
          aec->delay_estimator, abs_near_spectrum, PART_LEN1);
      if (delay_estimate >= 0) {
        // Update delay estimate buffer.
        aec->delay_histogram[delay_estimate]++;
        aec->num_delay_values++;
      }
      if (aec->delay_metrics_delivered == 1 &&
          aec->num_delay_values >= kDelayMetricsAggregationWindow) {
        UpdateDelayMetrics(aec);
      }
    }
  }
  #endif
//#endif
  // Perform echo subtraction.
#if 1
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{
  EchoSubtraction(
          aec->num_partitions, aec->extended_filter_enabled,
          &aec->extreme_filter_divergence, aec->filter_step_size,
          aec->error_threshold, farend_fft_hifi, &aec->xfBufBlockPos, aec->xfBuf_hifi,
          &nearend_block[0][0], aec->xPow, aec->wfBuf_hifi, echo_subtractor_output);
}
//iet_printf_msg("EchoSubtraction time =%d \r\n", osKernelGetTickCount()-t0);          
#endif
#if 0
  if (aec->metricsMode == 1) {
    UpdateLevel(&aec->linoutlevel,
                CalculatePower(echo_subtractor_output, PART_LEN));
  }
#endif


#if 0
if(aec->frame_count == 583)
{
    for(int j = 0; j < 64; j++)
    {
        printf("%d, ",echo_subtractor_output[j]);
        if((j+1)%64==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sd[j]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.se[j]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sx[j]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sde[j][0]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sde[j][1]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sxd[j][0]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    for(int j = 0; j < 65; j++)
    {
        printf("%lld, ",aec->coherence_state.sxd[j][1]);
        if((j+1)%65==0)
        {
            printf("\n");
        }
    }
    
    for(int j = 0; j < 128; j++)
    {
        printf("%d, ",aec->eBuf[j]);
        if((j+1)%64==0)
        {
            printf("\n");
        }
    }
}
#endif
//t0 = osKernelGetTickCount();
//for (int n = 0; n < 4*32; n++)
{
  // Perform echo suppression.
  EchoSuppression(aec, nearend_extended_block_lowest_band,farend_extended_block_lowest_band, echo_subtractor_output,output_block);
}
//iet_printf_msg("EchoSuppression time =%d \r\n", osKernelGetTickCount()-t0);
#if 0
  if (aec->metricsMode == 1) {
    UpdateLevel(&aec->nlpoutlevel,
                CalculatePower(&output_block[0][0], PART_LEN));
    UpdateMetrics(aec);
  }
#endif

  // for HIFI disable NLP
  #if 0
  for (i = 0; i < PART_LEN; i++)
  {
    output_block[0][i] = echo_subtractor_output[i]>>(-aecq0);//aecsq
  }
  #endif
  
  // Store the nearend signal until the next frame.
  for (i = 0; i < aec->num_bands; ++i) {
    memcpy(&aec->previous_nearend_block[i][0], &nearend_block[i][0],
           sizeof(int32_t) * PART_LEN);
  }

  // aec->data_dumper->DumpWav("aec_out", PART_LEN, &output_block[0][0],
  //                           std::min(aec->sampFreq, 16000), 1);
}
AecCore* WebRtcAec_CreateAec() {
  AecCore* aec = &AecCorestruct;

  aec->nearend_buffer_size = 0;
  memset(&aec->nearend_buffer[0], 0, sizeof(aec->nearend_buffer));
  // Start the output buffer with zeros to be able to produce
  // a full output frame in the first frame.
  //aec->output_buffer_size = PART_LEN - (FRAME_LEN - PART_LEN);
  memset(&aec->output_buffer[0], 0, sizeof(aec->output_buffer));

  aec->delay_estimator_farend =
      WebRtc_CreateDelayEstimatorFarend(PART_LEN1, kHistorySizeBlocks);
  if (aec->delay_estimator_farend == NULL) {
    WebRtcAec_FreeAec(aec);
    return NULL;
  }
  // We create the delay_estimator with the same amount of maximum lookahead as
  // the delay history size (kHistorySizeBlocks) for symmetry reasons.
  aec->delay_estimator = WebRtc_CreateDelayEstimator(
      aec->delay_estimator_farend, kHistorySizeBlocks);
  if (aec->delay_estimator == NULL) {
    WebRtcAec_FreeAec(aec);
    return NULL;
  }
#ifdef WEBRTC_ANDROID
  aec->delay_agnostic_enabled = 1;  // DA-AEC enabled by default.
  // DA-AEC assumes the system is causal from the beginning and will self adjust
  // the lookahead when shifting is required.
  WebRtc_set_lookahead(aec->delay_estimator, 0);
#else
  aec->delay_agnostic_enabled = 0;
  WebRtc_set_lookahead(aec->delay_estimator, kLookaheadBlocks);
#endif
  aec->extended_filter_enabled = 0;
  aec->refined_adaptive_filter_enabled = 0;

  // Assembly optimization
  //WebRtcAec_FilterFar = FilterFar;
  //WebRtcAec_ScaleErrorSignal = ScaleErrorSignal;
  //WebRtcAec_FilterAdaptation = FilterAdaptation;
  //WebRtcAec_Overdrive = Overdrive;
  //WebRtcAec_Suppress = Suppress;
  //WebRtcAec_ComputeCoherence = ComputeCoherence;
  //WebRtcAec_UpdateCoherenceSpectra = UpdateCoherenceSpectra;
  //WebRtcAec_StoreAsComplex = StoreAsComplex;
  //WebRtcAec_PartitionDelay = PartitionDelay;
  //WebRtcAec_WindowData = WindowData;

// #if defined(WEBRTC_ARCH_X86_FAMILY)
//   if (WebRtc_GetCPUInfo(kSSE2)) {
//     WebRtcAec_InitAec_SSE2();
//   }
// #endif

#if defined(MIPS_FPU_LE)
  WebRtcAec_InitAec_mips();
#endif

#if defined(WEBRTC_HAS_NEON)
  WebRtcAec_InitAec_neon();
#endif

  return aec;
}

void WebRtcAec_FreeAec(AecCore* aec) {
  WebRtc_FreeDelayEstimator(aec->delay_estimator);
  WebRtc_FreeDelayEstimatorFarend(aec->delay_estimator_farend);

  //delete aec;
}

static void SetAdaptiveFilterStepSize(AecCore* aec) {
  // Extended filter adaptation parameter.
  // TODO(ajm): No narrowband tuning yet.
  const int kExtendedMu = 858993459;//Q31

  if (aec->refined_adaptive_filter_enabled) {
    aec->filter_step_size = 107374182;
  } else {
    if (aec->extended_filter_enabled) {
      aec->filter_step_size = kExtendedMu;
    } else {
      if (aec->sampFreq == 8000) {
        aec->filter_step_size = 1288490189;
      } else {
        aec->filter_step_size = 1073741824;
      }
    }
  }
}

static void SetErrorThreshold(AecCore* aec) {
  // Extended filter adaptation parameter.
  // TODO(ajm): No narrowband tuning yet.
  static const int kExtendedErrorThreshold = 2147;//1.0e-6f; Q31

  if (aec->extended_filter_enabled) {
    aec->error_threshold = kExtendedErrorThreshold;
  } else {
    if (aec->sampFreq == 8000) {
      aec->error_threshold = 4295;//2e-6f;
    } else {
      aec->error_threshold = 3221;//1.5e-6f;
    }
  }
}

int WebRtcAec_InitAec(AecCore* aec, int sampFreq) {
  int i;
  //#ifdef HARD_FFT
  //#else
    arm_rfft_init_q31(&S_fft128, 128, 0, 1);
    arm_rfft_init_q31(&S_ifft128, 128, 1, 1); 
  //#endif  
    //aec->data_dumper->InitiateNewSetOfRecordings();
    buffer_ = WebRtc_CreateBuffer(kBufferSizeBlocks, sizeof(int32_t) * PART_LEN);
  aec->sampFreq = sampFreq;

  SetAdaptiveFilterStepSize(aec);
  SetErrorThreshold(aec);

  if (sampFreq == 8000) {
    aec->num_bands = 1;
  } else {
    aec->num_bands = (size_t)(sampFreq / 16000);
  }

  // Start the output buffer with zeros to be able to produce
  // a full output frame in the first frame.
  aec->output_buffer_size = PART_LEN - (FRAME_LEN - PART_LEN);//48;//
  memset(&aec->output_buffer[0], 0, sizeof(aec->output_buffer));
  aec->nearend_buffer_size = 0;
  memset(&aec->nearend_buffer[0], 0, sizeof(aec->nearend_buffer));

  // Initialize far-end buffer.
  ReInit();

  aec->system_delay = 0;

  if (WebRtc_InitDelayEstimatorFarend(aec->delay_estimator_farend) != 0) {
    return -1;
  }
  if (WebRtc_InitDelayEstimator(aec->delay_estimator) != 0) {
    return -1;
  }
  aec->delay_logging_enabled = 0;
  aec->delay_metrics_delivered = 0;
  memset(aec->delay_histogram, 0, sizeof(aec->delay_histogram));
  aec->num_delay_values = 0;
  aec->delay_median = -1;
  aec->delay_std = -1;
  //aec->fraction_poor_delays = -1.0f;

  aec->previous_delay = -2;  // (-2): Uninitialized.
  aec->delay_correction_count = 0;
  aec->shift_offset = kInitialShiftOffset;
  //aec->delay_quality_threshold = kDelayQualityThresholdMin;

  aec->num_partitions = kNormalNumPartitions;

  // Update the delay estimator with filter length.  We use half the
  // |num_partitions| to take the echo path into account.  In practice we say
  // that the echo has a duration of maximum half |num_partitions|, which is not
  // true, but serves as a crude measure.
  WebRtc_set_allowed_offset(aec->delay_estimator, aec->num_partitions / 2);
  // TODO(bjornv): I currently hard coded the enable.  Once we've established
  // that AECM has no performance regression, robust_validation will be enabled
  // all the time and the APIs to turn it on/off will be removed.  Hence, remove
  // this line then.
  WebRtc_enable_robust_validation(aec->delay_estimator, 1);
  aec->frame_count = 0;

  // Default target suppression mode.
  aec->nlp_mode = 1;

  // Sampling frequency multiplier w.r.t. 8 kHz.
  // In case of multiple bands we process the lower band in 16 kHz, hence the
  // multiplier is always 2.
  if (aec->num_bands > 1) {
    aec->mult = 2;
  } else {
    aec->mult = (int16_t)(aec->sampFreq) / 8000;
  }

  aec->farBufWritePos = 0;
  aec->farBufReadPos = 0;

  aec->inSamples = 0;
  aec->outSamples = 0;
  aec->knownDelay = 0;

  // Initialize buffers
  memset(aec->previous_nearend_block, 0, sizeof(aec->previous_nearend_block));
  memset(aec->eBuf, 0, sizeof(aec->eBuf));

  memset(aec->xPow, 0, sizeof(aec->xPow));
  memset(aec->dPow, 0, sizeof(aec->dPow));
  memset(aec->dInitMinPow, 0, sizeof(aec->dInitMinPow));
  //aec->noisePow = aec->dInitMinPow;
  aec->noiseEstCtr = 0;

  // Initial comfort noise power
  for (i = 0; i < PART_LEN1; i++) {
    aec->dMinPow[i] = 1e6;
  }

  // Holds the last block written to
  aec->xfBufBlockPos = 0;
  // TODO(peah): Investigate need for these initializations. Deleting them
  // doesn't change the output at all and yields 0.4% overall speedup.

   memset(aec->xfBuf_hifi, 0, sizeof(int32_t) * kExtendedNumPartitions * PART_LEN3);
   memset(aec->wfBuf_hifi, 0, sizeof(int32_t) * kExtendedNumPartitions * PART_LEN3);
   memset(aec->coherence_state.sde, 0, sizeof(complex_t) * PART_LEN1);
   memset(aec->coherence_state.sxd, 0, sizeof(complex_t) * PART_LEN1);
   memset(aec->xfwBuf_hifi, 0,
          sizeof(int32_t) * kExtendedNumPartitions * PART_LEN3);
   memset(aec->coherence_state.se, 0, sizeof(long long) * PART_LEN1);


  // To prevent numerical instability in the first block.
  for (i = 0; i < PART_LEN1; i++) {
    aec->coherence_state.sd[i] = 1;
  }
  for (i = 0; i < PART_LEN1; i++) {
    aec->coherence_state.sx[i] = 1;
  }

  //memset(aec->hNs, 0, sizeof(aec->hNs));
  memset(aec->outBuf, 0, sizeof(int32_t) * PART_LEN);

  aec->hNlFbMin = 65536;//Q16
  aec->hNlFbLocalMin = 65536;
  aec->hNlXdAvgMin = 65536;
  aec->hNlNewMin = 0;
  aec->hNlMinCtr = 0;
  aec->overDrive = 2;
  aec->overdrive_scaling = 2<<9;
  aec->delayIdx = 0;
  aec->stNearState = 0;
  aec->echoState = 0;
  aec->divergeState = 0;

  aec->seed = 777;
  aec->delayEstCtr = 0;

  aec->extreme_filter_divergence = 0;

  // Metrics disabled by default
  aec->metricsMode = 0;
  //InitMetrics(aec);
  // Extended Filter// for 128/160ms
  WebRtcAec_enable_extended_filter(aec, 1);//kNormalNumPartitions=12或4时，采用normal filter
  return 0;
}

void WebRtcAec_BufferFarendBlock(AecCore* aec, const int32_t* farend) {
  // Check if the buffer is full, and in that case flush the oldest data.
  if (AvaliableSpace() < 1) {
    AdjustSize(1);
  }
  Insert(farend);
}

int WebRtcAec_AdjustFarendBufferSizeAndSystemDelay(AecCore* aec,
                                                   int buffer_size_decrease) {
  int achieved_buffer_size_decrease =
      AdjustSize(buffer_size_decrease);
  aec->system_delay -= achieved_buffer_size_decrease * PART_LEN;
  return achieved_buffer_size_decrease;
}

void FormNearendBlock(
    size_t nearend_start_index,
    size_t num_bands,
    const int32_t* const* nearend_frame,
    size_t num_samples_from_nearend_frame,
    const int32_t nearend_buffer[NUM_HIGH_BANDS_MAX + 1]
                              [PART_LEN - (FRAME_LEN - PART_LEN)],
    int32_t nearend_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN]) {
  assert(num_samples_from_nearend_frame <= PART_LEN);
  const int num_samples_from_buffer = PART_LEN - num_samples_from_nearend_frame;

  if (num_samples_from_buffer > 0) {
    for (size_t i = 0; i < num_bands; ++i) {
        #if 1
        for(int j = 0; j<num_samples_from_buffer; j++)
        {
            nearend_block[i][j] = nearend_buffer[i][j];
        }
        #else
      memcpy(&nearend_block[i][0], &nearend_buffer[i][0],
             num_samples_from_buffer * sizeof(float));
        #endif
    }
  }

  for (size_t i = 0; i < num_bands; ++i) {
    #if 1
        for(int j = 0; j<num_samples_from_nearend_frame; j++)
        {
            nearend_block[i][num_samples_from_buffer+j] = nearend_frame[i][nearend_start_index+j];
        }
    #else
    memcpy(&nearend_block[i][num_samples_from_buffer],
           &nearend_frame[i][nearend_start_index],
           num_samples_from_nearend_frame * sizeof(float));
    #endif
  }
}

void BufferNearendFrame(
    size_t nearend_start_index,
    size_t num_bands,
    const int32_t* const* nearend_frame,
    size_t num_samples_to_buffer,
    int32_t nearend_buffer[NUM_HIGH_BANDS_MAX + 1]
                        [PART_LEN - (FRAME_LEN - PART_LEN)]) {
  for (size_t i = 0; i < num_bands; ++i) {
    memcpy(&nearend_buffer[i][0],
           &nearend_frame[i][nearend_start_index + FRAME_LEN -
                             num_samples_to_buffer],
           num_samples_to_buffer * sizeof(int32_t));
  }
}

void BufferOutputBlock(
    size_t num_bands,
    const int32_t output_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN],
    size_t* output_buffer_size,
    int32_t output_buffer[NUM_HIGH_BANDS_MAX + 1][2 * PART_LEN]) {
    for (size_t i = 0; i < num_bands; ++i) {
        #if 1
        for(int j = 0; j<PART_LEN; j++)
        {
            output_buffer[i][*output_buffer_size+j] = output_block[i][j];
        }
        #else
    memcpy(&output_buffer[i][*output_buffer_size], &output_block[i][0],
           PART_LEN * sizeof(float));
        #endif
    }
  (*output_buffer_size) += PART_LEN;
}

void FormOutputFrame(size_t output_start_index,
                     size_t num_bands,
                     size_t* output_buffer_size,
                     int32_t output_buffer[NUM_HIGH_BANDS_MAX + 1][2 * PART_LEN],
                     int32_t* const* output_frame) {
  assert(FRAME_LEN <= (*output_buffer_size));
  for (size_t i = 0; i < num_bands; ++i) {
    memcpy(&output_frame[i][output_start_index], &output_buffer[i][0],
           FRAME_LEN * sizeof(int32_t));
  }
  (*output_buffer_size) -= FRAME_LEN;
  if (*output_buffer_size > 0) {
    assert((2 * PART_LEN - FRAME_LEN) >= (*output_buffer_size));
    for (size_t i = 0; i < num_bands; ++i) {
      memcpy(&output_buffer[i][0], &output_buffer[i][FRAME_LEN],
             (*output_buffer_size) * sizeof(int32_t));
    }
  }
}

void WebRtcAec_ProcessFrames(AecCore* aec,
                             const int32_t* const* nearend,
                             size_t num_bands,
                             size_t num_samples,
                             int knownDelay,
                             int32_t* const* out) {
  //assert(num_samples == 80 || num_samples == 160);

  aec->frame_count++;
  // For each frame the process is as follows:
  // 1) If the system_delay indicates on being too small for processing a
  //    frame we stuff the buffer with enough data for 10 ms.
  // 2 a) Adjust the buffer to the system delay, by moving the read pointer.
  //   b) Apply signal based delay correction, if we have detected poor AEC
  //    performance.
  // 3) TODO(bjornv): Investigate if we need to add this:
  //    If we can't move read pointer due to buffer size limitations we
  //    flush/stuff the buffer.
  // 4) Process as many partitions as possible.
  // 5) Update the |system_delay| with respect to a full frame of FRAME_LEN
  //    samples. Even though we will have data left to process (we work with
  //    partitions) we consider updating a whole frame, since that's the
  //    amount of data we input and output in audio_processing.
  // 6) Update the outputs.

  // The AEC has two different delay estimation algorithms built in.  The
  // first relies on delay input values from the user and the amount of
  // shifted buffer elements is controlled by |knownDelay|.  This delay will
  // give a guess on how much we need to shift far-end buffers to align with
  // the near-end signal.  The other delay estimation algorithm uses the
  // far- and near-end signals to find the offset between them.  This one
  // (called "signal delay") is then used to fine tune the alignment, or
  // simply compensate for errors in the system based one.
  // Note that the two algorithms operate independently.  Currently, we only
  // allow one algorithm to be turned on.

  assert(aec->num_bands == num_bands);

  for (size_t j = 0; j < num_samples; j += FRAME_LEN) {
    // 1) At most we process |aec->mult|+1 partitions in 10 ms. Make sure we
    // have enough far-end data for that by stuffing the buffer if the
    // |system_delay| indicates others.
    if (aec->system_delay < FRAME_LEN) {
      // We don't have enough data so we rewind 10 ms.
      WebRtcAec_AdjustFarendBufferSizeAndSystemDelay(aec, -(aec->mult + 1));
    }

    if (!aec->delay_agnostic_enabled) {
      // 2 a) Compensate for a possible change in the system delay.

      // TODO(bjornv): Investigate how we should round the delay difference;
      // right now we know that incoming |knownDelay| is underestimated when
      // it's less than |aec->knownDelay|. We therefore, round (-32) in that
      // direction. In the other direction, we don't have this situation, but
      // might flush one partition too little. This can cause non-causality,
      // which should be investigated. Maybe, allow for a non-symmetric
      // rounding, like -16.
      int move_elements = (aec->knownDelay - knownDelay - 32) / PART_LEN;
      int moved_elements = AdjustSize(move_elements);
      // MaybeLogDelayAdjustment(moved_elements * (aec->sampFreq == 8000 ? 8 : 4),
      //                         DelaySource::kSystemDelay);
      aec->knownDelay -= moved_elements * PART_LEN;
    } 
    #if 0
    else {
      // 2 b) Apply signal based delay correction.
      int move_elements = SignalBasedDelayCorrection(aec);
      int moved_elements = aec->farend_block_buffer_.AdjustSize(move_elements);
      // MaybeLogDelayAdjustment(moved_elements * (aec->sampFreq == 8000 ? 8 : 4),
      //                         DelaySource::kDelayAgnostic);
      int far_near_buffer_diff =
          aec->farend_block_buffer_.Size() -
          (aec->nearend_buffer_size + FRAME_LEN) / PART_LEN;
      WebRtc_SoftResetDelayEstimator(aec->delay_estimator, moved_elements);
      WebRtc_SoftResetDelayEstimatorFarend(aec->delay_estimator_farend,
                                           moved_elements);
      // If we rely on reported system delay values only, a buffer underrun here
      // can never occur since we've taken care of that in 1) above.  Here, we
      // apply signal based delay correction and can therefore end up with
      // buffer underruns since the delay estimation can be wrong.  We therefore
      // stuff the buffer with enough elements if needed.
      if (far_near_buffer_diff < 0) {
        WebRtcAec_AdjustFarendBufferSizeAndSystemDelay(aec,
                                                       far_near_buffer_diff);
      }
    }
    #endif
//    static_assert(
//       16 == (FRAME_LEN - PART_LEN),
//       "These constants need to be properly related for this code to work");
    int32_t output_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN]={0};
    int32_t nearend_block[NUM_HIGH_BANDS_MAX + 1][PART_LEN];
    int32_t farend_extended_block_lowest_band[PART_LEN2];
	//if(j==192 && (aec->frame_count ==1))
	{
	    //init_flag_q = 1;
	}
    // Form and process a block of nearend samples, buffer the output block of
    // samples.
    ExtractExtendedBlock(farend_extended_block_lowest_band);
    FormNearendBlock(j, num_bands, nearend, PART_LEN - aec->nearend_buffer_size,aec->nearend_buffer, nearend_block);//copy PART_LEN from nearend_buffer/nearend to nearend_block
    ProcessNearendBlock(aec, farend_extended_block_lowest_band, nearend_block,output_block);
    BufferOutputBlock(num_bands, output_block, &aec->output_buffer_size,aec->output_buffer);//copy PART_LEN from output_block to output_buffer


    if ((FRAME_LEN - PART_LEN + aec->nearend_buffer_size) == PART_LEN) {
      // When possible (every fourth frame) form and process a second block of
      // nearend samples, buffer the output block of samples.
      ExtractExtendedBlock(
          farend_extended_block_lowest_band);
      FormNearendBlock(j + FRAME_LEN - PART_LEN, num_bands, nearend, PART_LEN,aec->nearend_buffer, nearend_block);
      ProcessNearendBlock(aec, farend_extended_block_lowest_band, nearend_block,output_block);
      BufferOutputBlock(num_bands, output_block, &aec->output_buffer_size,aec->output_buffer);

      // Reset the buffer size as there are no samples left in the nearend input
      // to buffer.
      aec->nearend_buffer_size = 0;
    } else {
      // Buffer the remaining samples in the nearend input.
      aec->nearend_buffer_size += FRAME_LEN - PART_LEN;
      BufferNearendFrame(j, num_bands, nearend, aec->nearend_buffer_size,
                         aec->nearend_buffer);
    }

    // 5) Update system delay with respect to the entire frame.
    aec->system_delay -= FRAME_LEN;

    // 6) Form the output frame.
    FormOutputFrame(j, num_bands, &aec->output_buffer_size, aec->output_buffer,
                    out);
  }
}
#if 0
int WebRtcAec_GetDelayMetricsCore(AecCore* self,
                                  int* median,
                                  int* std,
                                  float* fraction_poor_delays) {
  assert(self);
  assert(median);
  assert(std);

  if (self->delay_logging_enabled == 0) {
    // Logging disabled.
    return -1;
  }

  if (self->delay_metrics_delivered == 0) {
    UpdateDelayMetrics(self);
    self->delay_metrics_delivered = 1;
  }
  *median = self->delay_median;
  *std = self->delay_std;
  *fraction_poor_delays = self->fraction_poor_delays;

  return 0;
}

int WebRtcAec_echo_state(AecCore* self) {
  return self->echoState;
}

void WebRtcAec_GetEchoStats(AecCore* self,
                            Stats* erl,
                            Stats* erle,
                            Stats* a_nlp,
                            float* divergent_filter_fraction) {
  assert(erl);
  assert(erle);
  assert(a_nlp);
  *erl = self->erl;
  *erle = self->erle;
  *a_nlp = self->aNlp;
  *divergent_filter_fraction =
      self->divergent_filter_fraction.GetLatestFraction();
}
#endif
void WebRtcAec_SetConfigCore(AecCore* self,
                             int nlp_mode,
                             int metrics_mode,
                             int delay_logging) {
  assert(nlp_mode >= 0);
  assert(nlp_mode < 3);
  self->nlp_mode = nlp_mode;
  self->metricsMode = metrics_mode;
  if (self->metricsMode) {
    //InitMetrics(self);
  }
  // Turn on delay logging if it is either set explicitly or if delay agnostic
  // AEC is enabled (which requires delay estimates).
  self->delay_logging_enabled = delay_logging || self->delay_agnostic_enabled;
  if (self->delay_logging_enabled) {
    memset(self->delay_histogram, 0, sizeof(self->delay_histogram));
  }
}

void WebRtcAec_enable_delay_agnostic(AecCore* self, int enable) {
  self->delay_agnostic_enabled = enable;
}

int WebRtcAec_delay_agnostic_enabled(AecCore* self) {
  return self->delay_agnostic_enabled;
}

void WebRtcAec_enable_refined_adaptive_filter(AecCore* self, bool enable) {
  self->refined_adaptive_filter_enabled = enable;
  SetAdaptiveFilterStepSize(self);
  SetErrorThreshold(self);
}

bool WebRtcAec_refined_adaptive_filter_enabled(const AecCore* self) {
  return self->refined_adaptive_filter_enabled;
}

void WebRtcAec_enable_extended_filter(AecCore* self, int enable) {
  self->extended_filter_enabled = enable;
  SetAdaptiveFilterStepSize(self);
  SetErrorThreshold(self);
  self->num_partitions = enable ? kExtendedNumPartitions : kNormalNumPartitions;
  // Update the delay estimator with filter length.  See InitAEC() for details.
  WebRtc_set_allowed_offset(self->delay_estimator, self->num_partitions / 2);
}

int WebRtcAec_extended_filter_enabled(AecCore* self) {
  return self->extended_filter_enabled;
}

int WebRtcAec_system_delay(AecCore* self) {
  return self->system_delay;
}

void WebRtcAec_SetSystemDelay(AecCore* self, int delay) {
  assert(delay >= 0);
  self->system_delay = delay;
}
//}  // namespace webrtc
