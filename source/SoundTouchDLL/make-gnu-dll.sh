#!/bin/sh
#
# This script compiles SoundTouch dynamic-link library for GNU environment 
# with wrapper functions that are easier to import to Java / Mono / etc
#

g++ -fPIC -shared -DDLL_EXPORTS -fvisibility=hidden -I../../include -o SoundTouchDll.so \
     SoundTouchDLL.cpp ../SoundTouch/*.cpp -I../SoundTouch -O3 -msse
