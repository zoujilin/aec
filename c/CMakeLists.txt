cmake_minimum_required(VERSION 3.5)
project(aec)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
add_executable(aec main.c fft_c.cpp fft.cpp dtype.cpp aec_core.c aec_resampler.c delay_estimator.c
 delay_estimator_wrapper.c dot_product_with_scale.c echo_cancellation.c
 main_aec.c ring_buffer.c 
 arm_bitreversal.c arm_cfft_q31.c arm_cfft_radix2_init_q31.c arm_cfft_radix2_q31.c 
 arm_cfft_radix4_init_q31.c arm_cfft_radix4_q31.c arm_common_tables.c arm_const_structs.c 
 arm_rfft_init_q31.c arm_rfft_q31.c
 dtype.h fft.h fft_c.h fft_define.h aec_common.h aec_resampler.h arch.h audio_aec.h test_aec_data.h
 delay_estimator.h delay_estimator_internal.h delay_estimator_wrapper.h dot_product_with_scale.h
 dr_wav.h echo_cancellation.h ring_buffer.h signal_processing_library.h timing.h)

target_link_libraries(aec m)
