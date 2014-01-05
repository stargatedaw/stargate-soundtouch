////////////////////////////////////////////////////////////////////////////////
/// 
/// Sample rate transposer. Changes sample rate by using linear interpolation 
/// together with anti-alias filtering (first order interpolation with anti-
/// alias filtering should be quite adequate for this application)
///
/// Author        : Copyright (c) Olli Parviainen
/// Author e-mail : oparviai 'at' iki.fi
/// SoundTouch WWW: http://www.surina.net/soundtouch
///
////////////////////////////////////////////////////////////////////////////////
//
// Last changed  : $Date$
// File revision : $Revision: 4 $
//
// $Id$
//
////////////////////////////////////////////////////////////////////////////////
//
// License :
//
//  SoundTouch audio processing library
//  Copyright (c) Olli Parviainen
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
////////////////////////////////////////////////////////////////////////////////

#include <memory.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "RateTransposer.h"
#include "AAFilter.h"

using namespace soundtouch;



/// A linear samplerate transposer class that uses integer arithmetics.
/// for the transposing.
class LinearTransposerBase: public TransposerBase
{
protected:
    virtual int transposeStereo(SAMPLETYPE *dest, 
                        const SAMPLETYPE *src, 
                        uint numSamples) = 0;
    virtual int transposeMono(SAMPLETYPE *dest, 
                        const SAMPLETYPE *src, 
                        uint numSamples)  = 0;
    virtual int transposeMulti(SAMPLETYPE *dest, 
                        const SAMPLETYPE *src, 
                        uint numSamples) = 0;
public:
    virtual int transpose(FIFOSampleBuffer &dest, FIFOSampleBuffer &src);

    static LinearTransposerBase *newInstance();
};


/// A linear samplerate transposer class that uses integer arithmetics.
/// for the transposing.
class LinearTransposerInteger : public LinearTransposerBase
{
protected:
    int iSlopeCount;
    int iRate;
    SAMPLETYPE *sPrevSample;

    virtual void resetRegisters();

    virtual int transposeStereo(SAMPLETYPE *dest, 
                         const SAMPLETYPE *src, 
                         uint numSamples);
    virtual int transposeMono(SAMPLETYPE *dest, 
                       const SAMPLETYPE *src, 
                       uint numSamples);
    virtual int transposeMulti(SAMPLETYPE *dest, const SAMPLETYPE *src, uint numSamples);
public:
    LinearTransposerInteger();
    virtual ~LinearTransposerInteger();

    /// Sets new target rate. Normal rate = 1.0, smaller values represent slower 
    /// rate, larger faster rates.
    virtual void setRate(float newRate);
};


/// A linear samplerate transposer class that uses floating point arithmetics
/// for the transposing.
class LinearTransposerFloat : public LinearTransposerBase
{
protected:
    float fSlopeCount;
    SAMPLETYPE *sPrevSample;

    virtual void resetRegisters();

    virtual int transposeStereo(SAMPLETYPE *dest, 
                         const SAMPLETYPE *src, 
                         uint numSamples);
    virtual int transposeMono(SAMPLETYPE *dest, 
                       const SAMPLETYPE *src, 
                       uint numSamples);
    virtual int transposeMulti(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples);

public:
    LinearTransposerFloat();
    virtual ~LinearTransposerFloat();
};


TransposerBase::TransposerBase()
{
    numChannels = 0;
    rate = 1.0f;
}


TransposerBase::~TransposerBase()
{
}


void TransposerBase::setChannels(int channels)
{
    numChannels = channels;
    resetRegisters();
}


void TransposerBase::setRate(float newRate)
{
    rate = newRate;
}


// Operator 'new' is overloaded so that it automatically creates a suitable instance 
// depending on if we've a MMX/SSE/etc-capable CPU available or not.
/*
void * RateTransposer::operator new(size_t s)
{
    ST_THROW_RT_ERROR("Error in RateTransoser::new: don't use \"new TDStretch\" directly, use \"newInstance\" to create a new instance instead!");
    return newInstance();
}
*/


TransposerBase *TransposerBase::newInstance()
{
#ifdef SOUNDTOUCH_INTEGER_SAMPLES
    return ::new LinearTransposerInteger;
#else
    return ::new LinearTransposerFloat;
#endif
}


// Constructor
RateTransposer::RateTransposer() : FIFOProcessor(&outputBuffer)
{
    bUseAAFilter = TRUE;

    // Instantiates the anti-alias filter with default tap length
    // of 32
    pAAFilter = new AAFilter(32);
    pTransposer = TransposerBase::newInstance();
}



RateTransposer::~RateTransposer()
{
    delete pAAFilter;
    delete pTransposer;
}



/// Enables/disables the anti-alias filter. Zero to disable, nonzero to enable
void RateTransposer::enableAAFilter(BOOL newMode)
{
    bUseAAFilter = newMode;
}


/// Returns nonzero if anti-alias filter is enabled.
BOOL RateTransposer::isAAFilterEnabled() const
{
    return bUseAAFilter;
}


AAFilter *RateTransposer::getAAFilter()
{
    return pAAFilter;
}



// Sets new target iRate. Normal iRate = 1.0, smaller values represent slower 
// iRate, larger faster iRates.
void RateTransposer::setRate(float newRate)
{
    double fCutoff;

    pTransposer->setRate(newRate);

    // design a new anti-alias filter
    if (newRate > 1.0f) 
    {
        fCutoff = 0.5f / newRate;
    } 
    else 
    {
        fCutoff = 0.5f * newRate;
    }
    pAAFilter->setCutoffFreq(fCutoff);
}


// Outputs as many samples of the 'outputBuffer' as possible, and if there's
// any room left, outputs also as many of the incoming samples as possible.
// The goal is to drive the outputBuffer empty.
//
// It's allowed for 'output' and 'input' parameters to point to the same
// memory position.
/*
void RateTransposer::flushinputBuffer()
{
    if (inputBuffer.isEmpty()) return;

    outputBuffer.moveSamples(inputBuffer);
}
*/


// Adds 'nSamples' pcs of samples from the 'samples' memory position into
// the input of the object.
void RateTransposer::putSamples(const SAMPLETYPE *samples, uint nSamples)
{
    processSamples(samples, nSamples);
}


// Transposes sample rate by applying anti-alias filter to prevent folding. 
// Returns amount of samples returned in the "dest" buffer.
// The maximum amount of samples that can be returned at a time is set by
// the 'set_returnBuffer_size' function.
void RateTransposer::processSamples(const SAMPLETYPE *src, uint nSamples)
{
    uint count;

    if (nSamples == 0) return;

    // Store samples to input buffer
    inputBuffer.putSamples(src, nSamples);

    // If anti-alias filter is turned off, simply transpose without applying
    // the filter
    if (bUseAAFilter == FALSE) 
    {
        count = pTransposer->transpose(outputBuffer, inputBuffer);
        return;
    }

    assert(pAAFilter);

    // Transpose with anti-alias filter
    if (pTransposer->rate < 1.0f) 
    {
        // If the parameter 'Rate' value is smaller than 1, first transpose
        // the samples and then apply the anti-alias filter to remove aliasing.

        // Transpose the samples, store the result to end of "midBuffer"
        pTransposer->transpose(midBuffer, inputBuffer);

        // Apply the anti-alias filter for transposed samples in midBuffer
        pAAFilter->evaluate(outputBuffer, midBuffer);
    } 
    else  
    {
        // If the parameter 'Rate' value is larger than 1, first apply the
        // anti-alias filter to remove high frequencies (prevent them from folding
        // over the lover frequencies), then transpose.

        // Apply the anti-alias filter for samples in inputBuffer
        pAAFilter->evaluate(midBuffer, inputBuffer);

        // Transpose the AA-filtered samples in "midBuffer"
        pTransposer->transpose(outputBuffer, midBuffer);
    }
}


// Transposes the sample rate of the given samples using linear interpolation. 
// Returns the number of samples returned in the "dest" buffer
int LinearTransposerBase::transpose(FIFOSampleBuffer &dest, FIFOSampleBuffer &src)
{
    int numSrcSamples = src.numSamples();
    int sizeDemand = (int)((float)numSrcSamples / rate) + 8;
    int numOutput;
    SAMPLETYPE *psrc = src.ptrBegin();
    SAMPLETYPE *pdest = dest.ptrEnd(sizeDemand);

#ifndef USE_MULTICH_ALWAYS
    if (numChannels == 1)
    {
        numOutput = transposeMono(pdest, psrc, numSrcSamples);
    }
    else if (numChannels == 2) 
    {
        numOutput = transposeStereo(pdest, psrc, numSrcSamples);
    } 
    else 
#endif // USE_MULTICH_ALWAYS
    {
        assert(numChannels > 0);
        numOutput = transposeMulti(pdest, psrc, numSrcSamples);
    }
    dest.putSamples(numOutput);
    src.receiveSamples(numOutput);
    return numOutput;
}


// Sets the number of channels, 1 = mono, 2 = stereo
void RateTransposer::setChannels(int nChannels)
{
    assert(nChannels > 0);

    if (pTransposer->numChannels == nChannels) return;
    pTransposer->setChannels(nChannels);

    inputBuffer.setChannels(nChannels);
    midBuffer.setChannels(nChannels);
    outputBuffer.setChannels(nChannels);
}


// Clears all the samples in the object
void RateTransposer::clear()
{
    outputBuffer.clear();
    midBuffer.clear();
    inputBuffer.clear();
}


// Returns nonzero if there aren't any samples available for outputting.
int RateTransposer::isEmpty() const
{
    int res;

    res = FIFOProcessor::isEmpty();
    if (res == 0) return 0;
    return inputBuffer.isEmpty();
}


//////////////////////////////////////////////////////////////////////////////
//
// LinearTransposerInteger - integer arithmetic implementation
// 

/// fixed-point interpolation routine precision
#define SCALE    65536

// Constructor
LinearTransposerInteger::LinearTransposerInteger() : LinearTransposerBase()
{
    // Notice: use local function calling syntax for sake of clarity, 
    // to indicate the fact that C++ constructor can't call virtual functions.
    sPrevSample=0;
    resetRegisters();
    setRate(1.0f);
}


LinearTransposerInteger::~LinearTransposerInteger()
{
    if (sPrevSample) delete[] sPrevSample;
}


void LinearTransposerInteger::resetRegisters()
{
    iSlopeCount = 0;
    delete[] sPrevSample;
    sPrevSample = new SAMPLETYPE[numChannels];
    memset(sPrevSample, 0, numChannels * sizeof(SAMPLETYPE));
}



// Transposes the sample rate of the given samples using linear interpolation. 
// 'Mono' version of the routine. Returns the number of samples returned in 
// the "dest" buffer
int LinearTransposerInteger::transposeMono(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remain;
    LONG_SAMPLETYPE temp, vol1;

    if (nSamples == 0) return 0;  // no samples, no work

    remain = nSamples - 1;
    i = 0;

    // Process the last sample saved from the previous call first...
    while (iSlopeCount <= SCALE) 
    {
        vol1 = (LONG_SAMPLETYPE)(SCALE - iSlopeCount);
        temp = vol1 * sPrevSample[0] + iSlopeCount * src[0];
        dest[i] = (SAMPLETYPE)(temp / SCALE);
        i++;
        iSlopeCount += iRate;
    }
    // now always (iSlopeCount > SCALE)
    iSlopeCount -= SCALE;

    while (1)
    {
        while (iSlopeCount > SCALE) 
        {
            iSlopeCount -= SCALE;
            src ++;
            remain --;
            if (remain == 0) goto end;
        }
        vol1 = (LONG_SAMPLETYPE)(SCALE - iSlopeCount);
        temp = src[0] * vol1 + iSlopeCount * src[1];
        dest[i] = (SAMPLETYPE)(temp / SCALE);

        i++;
        iSlopeCount += iRate;
    }
end:
    // Store the last sample for the next round
    sPrevSample[0] = src[0];

    return i;
}


// Transposes the sample rate of the given samples using linear interpolation. 
// 'Stereo' version of the routine. Returns the number of samples returned in 
// the "dest" buffer
int LinearTransposerInteger::transposeStereo(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remain;
    LONG_SAMPLETYPE temp, vol1;

    if (nSamples == 0) return 0;  // no samples, no work

    remain = nSamples - 1;
    i = 0;

    // Process the last sample saved from the sPrevSampleLious call first...
    while (iSlopeCount <= SCALE) 
    {
        vol1 = (LONG_SAMPLETYPE)(SCALE - iSlopeCount);
        temp = vol1 * sPrevSample[0] + iSlopeCount * src[0];
        dest[2 * i] = (SAMPLETYPE)(temp / SCALE);
        temp = vol1 * sPrevSample[1] + iSlopeCount * src[1];
        dest[2 * i + 1] = (SAMPLETYPE)(temp / SCALE);
        i++;
        iSlopeCount += iRate;
    }
    // now always (iSlopeCount > SCALE)
    iSlopeCount -= SCALE;

    while (1)
    {
        while (iSlopeCount > SCALE) 
        {
            iSlopeCount -= SCALE;
            remain --;
            src += 2;
            if (remain == 0) goto end;
        }
        vol1 = (LONG_SAMPLETYPE)(SCALE - iSlopeCount);
        temp = src[0] * vol1 + iSlopeCount * src[2];
        dest[2 * i] = (SAMPLETYPE)(temp / SCALE);
        temp = src[1] * vol1 + iSlopeCount * src[3];
        dest[2 * i + 1] = (SAMPLETYPE)(temp / SCALE);

        i++;
        iSlopeCount += iRate;
    }
end:
    // Store the last sample for the next round
    sPrevSample[0] = src[0];
    sPrevSample[1] = src[1];

    return i;
}


int LinearTransposerInteger::transposeMulti(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remaining;
    LONG_SAMPLETYPE temp, vol1;

    if (nSamples == 0) return 0;  // no samples, no work

    remaining = nSamples - 1;
    i = 0;

    // Process the last sample saved from the sPrevSampleLious call first...
    while (iSlopeCount <= SCALE)
    {
        for (int c = 0; c < numChannels; c ++)
        {
            vol1 = (SCALE - iSlopeCount);
            temp = vol1 * sPrevSample[c] + iSlopeCount * src[c];
            *dest = (SAMPLETYPE)(temp / SCALE);
            dest ++;
        }
        i++;

        iSlopeCount += iRate;
    }
    // now always (iSlopeCount > SCALE)
    iSlopeCount -= SCALE;

    while (1)
    {
        while (iSlopeCount > SCALE) 
        {
            iSlopeCount -= SCALE;
            src += numChannels;
            remaining --;
            if (remaining == 0) goto end;
        }

        for (int c = 0; c < numChannels; c ++)
        {
            vol1 = (SCALE - iSlopeCount);
            temp = src[c] * vol1 + iSlopeCount * src[c + numChannels];
            *dest = (SAMPLETYPE)(temp / SCALE);
            dest++;
        }

        i++;
        iSlopeCount += iRate;
    }
end:
    // Store the last sample for the next round
    memcpy(sPrevSample, src, numChannels * sizeof(SAMPLETYPE));

    return i;
}

// Sets new target iRate. Normal iRate = 1.0, smaller values represent slower 
// iRate, larger faster iRates.
void LinearTransposerInteger::setRate(float newRate)
{
    iRate = (int)(newRate * SCALE + 0.5f);
    TransposerBase::setRate(newRate);
}


//////////////////////////////////////////////////////////////////////////////
//
// LinearTransposerFloat - floating point arithmetic implementation
// 
//////////////////////////////////////////////////////////////////////////////

// Constructor
LinearTransposerFloat::LinearTransposerFloat() : LinearTransposerBase()
{
    // Notice: use local function calling syntax for sake of clarity, 
    // to indicate the fact that C++ constructor can't call virtual functions.
    sPrevSample = NULL;
    resetRegisters();
    setRate(1.0f);
}


LinearTransposerFloat::~LinearTransposerFloat()
{
    delete[] sPrevSample;
}


void LinearTransposerFloat::resetRegisters()
{
    fSlopeCount = 0;
    delete[] sPrevSample;
    sPrevSample = new SAMPLETYPE[numChannels];
    memset(sPrevSample, 0, numChannels * sizeof(SAMPLETYPE));
}



// Transposes the sample rate of the given samples using linear interpolation. 
// 'Mono' version of the routine. Returns the number of samples returned in 
// the "dest" buffer
int LinearTransposerFloat::transposeMono(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remain;

    remain = nSamples - 1;
    i = 0;

    // Process the last sample saved from the previous call first...
    while (fSlopeCount <= 1.0f) 
    {
        dest[i] = (SAMPLETYPE)((1.0f - fSlopeCount) * sPrevSample[0] + fSlopeCount * src[0]);
        i++;
        fSlopeCount += rate;
    }
    fSlopeCount -= 1.0f;

    if (nSamples > 1)
    {
        while (1)
        {
            while (fSlopeCount > 1.0f) 
            {
                fSlopeCount -= 1.0f;
                src ++;
                remain --;
                if (remain == 0) goto end;
            }
            dest[i] = (SAMPLETYPE)((1.0f - fSlopeCount) * src[0] + fSlopeCount * src[1]);
            i++;
            fSlopeCount += rate;
        }
    }
end:
    // Store the last sample for the next round
    sPrevSample[0] = src[0];

    return i;
}


// Transposes the sample rate of the given samples using linear interpolation. 
// 'Mono' version of the routine. Returns the number of samples returned in 
// the "dest" buffer
int LinearTransposerFloat::transposeStereo(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remain;

    if (nSamples == 0) return 0;  // no samples, no work

    remain = nSamples - 1;
    i = 0;

    // Process the last sample saved from the sPrevSampleLious call first...
    while (fSlopeCount <= 1.0f) 
    {
        dest[2 * i] = (SAMPLETYPE)((1.0f - fSlopeCount) * sPrevSample[0] + fSlopeCount * src[0]);
        dest[2 * i + 1] = (SAMPLETYPE)((1.0f - fSlopeCount) * sPrevSample[1] + fSlopeCount * src[1]);
        i++;
        fSlopeCount += rate;
    }
    // now always (iSlopeCount > 1.0f)
    fSlopeCount -= 1.0f;

    if (nSamples > 1)
    {
        while (1)
        {
            while (fSlopeCount > 1.0f) 
            {
                fSlopeCount -= 1.0f;
                remain --;
                src += 2;
                if (remain == 0) goto end;
            }

            dest[2 * i] = (SAMPLETYPE)((1.0f - fSlopeCount) * src[0] 
                + fSlopeCount * src[2]);
            dest[2 * i + 1] = (SAMPLETYPE)((1.0f - fSlopeCount) * src[1] 
                + fSlopeCount * src[3]);

            i++;
            fSlopeCount += rate;
        }
    }
end:
    // Store the last sample for the next round
    sPrevSample[0] = src[0];
    sPrevSample[1] = src[1];

    return i;
}

int LinearTransposerFloat::transposeMulti(SAMPLETYPE *dest, const SAMPLETYPE *src, uint nSamples)
{
    int i, remaining;

    if (nSamples == 0) return 0;  // no samples, no work

    remaining = nSamples - 1;
    i = 0;

    // Process the last sample saved from the sPrevSampleLious call first...
    while (fSlopeCount <= 1.0f) 
    {
        for (int c = 0; c < numChannels; c ++)
        {
            *dest = (SAMPLETYPE)((1.0f - fSlopeCount) * sPrevSample[c] + fSlopeCount * src[c]);
            dest ++;
        }
        i++;
        fSlopeCount += rate;
    }
    // now always (iSlopeCount > 1.0f)
    fSlopeCount -= 1.0f;

    while (remaining > 0)
    {
        while (fSlopeCount > 1.0f) 
        {
            fSlopeCount -= 1.0f;
            src += numChannels;
            remaining --;
            if (remaining == 0) goto end;
        }

        for (int c = 0; c < numChannels; c ++)
        {
            *dest = (SAMPLETYPE)((1.0f - fSlopeCount) * src[c] 
                + fSlopeCount * src[c + numChannels]);
            dest++;
        }

        i++;
        fSlopeCount += rate;
    }

end:
    // Store the last sample for the next round
    memcpy(sPrevSample, src, numChannels * sizeof(SAMPLETYPE));

    return i;
}
