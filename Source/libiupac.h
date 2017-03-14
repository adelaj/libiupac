#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef const uint8_t *iupacPattern;
typedef uint8_t iupacEncodeLookup, iupacDecodeLookup;
typedef uint32_t iupacInclusionMask;

typedef void iupacStorage;

// Required by 'iupacEncodeLookupCreate'
inline iupacPattern *iupacPatternStd()
{
    extern const uint8_t *patternStd[];
    return (iupacPattern *) patternStd;
}

// Required by 'iupacDecode'
inline iupacDecodeLookup *iupacDecodeLookupStd()
{
    extern const uint8_t *patternStd[];
    return (iupacDecodeLookup *) patternStd[0];
}

// Nucleotide-wise test for inclusion ('x' in 'y'). Result is stored into a bit vector 'mask'.
void iupacTestInclusionMask(iupacStorage *x, iupacStorage *y, iupacInclusionMask *mask, size_t cnt);

// Test for inclusion ('x' in 'y'). Returns 'true' if inclusion holds. 
bool iupacTestInclusion(iupacStorage *x, iupacStorage *y, size_t cnt);

// Reverse complement of 'src'. Result is stored into 'dest'. The residual part of 'dest' is preserved.
void iupacReverseComplement(iupacStorage *dest, iupacStorage *src, size_t cnt);

// Encoding pattern which may be generated only once . Required by 'iupacEncode'.
iupacEncodeLookup *iupacEncodeLookupCreate(iupacPattern *pat);

//  Storage allocation for 'cnt' nucleotides
iupacStorage *iupacStorageCreate(size_t cnt);

// Storage disposal
void iupacStorageDispose(iupacStorage *stor);

// Convert string of length 'cnt' to bit representation. Residual of 'dest' remains unchanged.
void iupacEncode(uint8_t *src, iupacStorage *dest, size_t cnt, iupacEncodeLookup *lookup);

// Convert bit representation of 'cnt' nucleotides to string. 
void iupacDecode(iupacStorage *src, uint8_t *dest, size_t cnt, iupacDecodeLookup *lookup);
