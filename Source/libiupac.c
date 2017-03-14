#include "libiupac.h"
#include <stdlib.h>

#ifdef _MSC_VER
#   pragma warning(disable : 4200)
#   define alignof __alignof
#   define aligned_alloc(OFF, SZ) _aligned_malloc(SZ, OFF)
#   define aligned_free(PTR) _aligned_free(PTR)
#elif defined __GNUC__ || defined __clang__
#   include <stdalign.h>
#   define aligned_free(PTR) free(PTR) 
#endif

#define UINT8(X) ((uint8_t *) (X))

#ifndef __AVX__
#   define __AVX__
#endif

#ifdef __AVX__

#   include <immintrin.h>
#   define M128I(X) ((__m128i *) (X))

iupacStorage *iupacStorageCreate(size_t cnt)
{
    return aligned_alloc(alignof(__m128), ((cnt + 31) >> 5) * sizeof(__m128));
}

void iupacStorageDispose(iupacStorage *stor)
{
    aligned_free(stor);
}

void iupacTestInclusionMask(iupacStorage *x, iupacStorage *y, iupacInclusionMask *mask, size_t cnt)
{
    const __m128i
        x01 = _mm_set1_epi8(0x01),
        x0f = _mm_set1_epi8(0x0f),
        sha = _mm_set_epi32(0x0f0b0703, 0x0d090501, 0x0e0a0602, 0x0c080400);

    const size_t ceil = (cnt + 31) >> 5, rem = cnt & 31;
        
    for (size_t i = 0; i < ceil; i++)
    {
        __m128i 
            a = _mm_andnot_si128(_mm_load_si128(M128I(y) + i), _mm_load_si128(M128I(x) + i)),
            c = _mm_shuffle_epi8(
                    _mm_or_si128(
                        _mm_min_epu8(_mm_and_si128(a, x0f), x01),
                        _mm_slli_epi64(_mm_min_epu8(_mm_and_si128(_mm_srli_epi64(a, 4), x0f), x01), 1)
                    ),
                    sha
                );

        c = _mm_or_si128(_mm_slli_epi64(_mm_shuffle_epi32(c, 0x4e), 2), c);

        mask[i] = ~(uint32_t) _mm_cvtsi128_si32(_mm_or_si128(_mm_srli_epi64(c, 28), c));        
    }

    if (rem) mask[ceil - 1] &= (uint32_t) 0xffffffff >> (32 - rem);
}

static inline const uint8_t *partialMask(size_t off, size_t odd, size_t lo)
{
    static const uint8_t mask[] = {
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    return mask + (odd << 5) + (lo << 4) + odd + off;
}

bool iupacTestInclusion(iupacStorage *x, iupacStorage *y, size_t cnt)
{
    const size_t floor = cnt >> 5, rem = cnt & 31;

    for (size_t i = 0; i < floor; i++) if (!_mm_testc_si128(_mm_load_si128(M128I(y) + i), _mm_load_si128(M128I(x) + i))) return 0;

    if (rem)
    {
        return _mm_testc_si128(
            _mm_load_si128(M128I(y) + floor),
            _mm_and_si128(
                _mm_load_si128(M128I(x) + floor),
                _mm_loadu_si128((__m128i *) partialMask((32 - rem) >> 1, rem & 1, 1))
            )
        );       
    }

    return 1;
}

#   ifdef _MSC_VER
#       pragma warning(push)
#       pragma warning(disable : 4701)
#   endif

void iupacReverseComplement(iupacStorage *dest, iupacStorage *src, size_t cnt)
{
    const __m128i
        x0f = _mm_set1_epi8(0x0f),
        inv = _mm_set_epi8(0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f), // reversion shuffle
        pat = _mm_set_epi8(0x0f, 0x07, 0x0b, 0x03, 0x0d, 0x05, 0x09, 0x01, 0x0e, 0x06, 0x0a, 0x02, 0x0c, 0x04, 0x08, 0x00); // complement shuffle

    __m128i hi; // residual part of unit storage  

    const size_t ceil = (cnt + 31) >> 5, hceil = ceil >> 1, rem = cnt & 31;
    const size_t restb = (32 - rem) >> 1, floor = ceil - 1, remp = rem & 1;

    if (rem) hi = _mm_and_si128(_mm_load_si128(M128I(dest) + floor), _mm_loadu_si128(M128I(partialMask(restb, remp, 0))));
    
    for (size_t i = 0; i < hceil; i++)
    {
        __m128i a = M128I(src)[i], b = M128I(src)[floor - i];

        _mm_store_si128(
            M128I(dest) + floor - i,
            _mm_shuffle_epi8(
                _mm_or_si128(
                    _mm_shuffle_epi8(pat, _mm_and_si128(_mm_srli_epi64(a, 4), x0f)),
                    _mm_slli_epi64(_mm_shuffle_epi8(pat, _mm_and_si128(a, x0f)), 4)
                ),
                inv
            )
        );

        _mm_store_si128(
            M128I(dest) + i,
            _mm_shuffle_epi8(
                _mm_or_si128(
                    _mm_shuffle_epi8(pat, _mm_and_si128(_mm_srli_epi64(b, 4), x0f)),
                    _mm_slli_epi64(_mm_shuffle_epi8(pat, _mm_and_si128(b, x0f)), 4)
                ),
                inv
            )
        );
    }

    if (ceil & 1)
    {
        __m128i a = M128I(src)[hceil];

        _mm_store_si128(
            M128I(dest) + hceil,
            _mm_shuffle_epi8(
                _mm_or_si128(
                    _mm_shuffle_epi8(pat, _mm_and_si128(_mm_srli_epi64(a, 4), x0f)),
                    _mm_slli_epi64(_mm_shuffle_epi8(pat, _mm_and_si128(a, x0f)), 4)
                ),
                inv
            )
        ); 
    }
    
    if (rem)
    {
        static const uint8_t mask[] = {
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
            0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
        };
        
        if (remp) // Hard part
        {
            uint8_t lo = UINT8(dest)[restb] >> 4;
            
            for (size_t i = 0; i < floor; i++)
            {
                __m128i a = _mm_loadu_si128(M128I(UINT8(M128I(dest) + i) + restb + 1));
                                
                _mm_store_si128(
                    M128I(dest) + i,
                    _mm_or_si128(
                        _mm_insert_epi8(_mm_slli_si128(_mm_and_si128(_mm_srli_epi64(a, 4), x0f), 1), lo, 0),
                        _mm_slli_epi64(_mm_and_si128(a, x0f), 4)
                    )
                );

                lo = UINT8(M128I(dest) + i + 1)[restb] >> 4;
            }

            __m128i a = _mm_shuffle_epi8(_mm_load_si128(M128I(dest) + floor), _mm_loadu_si128(M128I(mask + restb + 1)));

            _mm_store_si128(
                M128I(dest) + floor,
                _mm_or_si128(
                    _mm_or_si128(
                        _mm_insert_epi8(_mm_slli_si128(_mm_and_si128(_mm_srli_epi64(a, 4), x0f), 1), lo, 0),
                        _mm_slli_epi64(_mm_and_si128(a, x0f), 4)
                    ),
                    hi
                )
            );
        }
        else // Easy part
        {
            for (size_t i = 0; i < floor; i++) _mm_store_si128(M128I(dest) + i, _mm_loadu_si128(M128I(UINT8(M128I(dest) + i) + restb)));
            
            _mm_store_si128(
                M128I(dest) + floor,
                _mm_or_si128(
                    _mm_shuffle_epi8(
                        _mm_load_si128(M128I(dest) + floor),
                        _mm_loadu_si128(M128I(mask + restb))
                    ),
                    hi
                )
            );
        }        
    }
}

#   ifdef _MSC_VER
#       pragma warning(pop)
#   endif

#else

struct iupacStorage { uint64_t ptr[]; };

#endif

const uint8_t *patternStd[] = {
    (uint8_t *) ".ACMGRSVTWYHKDBN",
    (uint8_t *) "-acmgrsvtwyhkdbn",
    (uint8_t *) "\0\x08Uu",
    NULL
};

extern inline iupacPattern *iupacPatternStd();
extern inline iupacDecodeLookup *iupacDecodeLookupStd();

iupacEncodeLookup *iupacEncodeLookupCreate(iupacPattern *pat)
{
    uint8_t *res = (uint8_t *) calloc(256, sizeof *res);
    if (!res) return NULL;

    for (const uint8_t **s = pat; *s; s++)
    {
        const uint8_t *t = *s;
        uint8_t i = 0;

        if (*t) do res[t[i]] = i; while (t[++i]);
        else while (t[i + 2]) res[t[i + 2]] = t[1], i++;
    }

    return res;
}

void iupacEncode(uint8_t *src, iupacStorage *dest, size_t cnt, iupacEncodeLookup *lookup)
{
    for (size_t pos = 0; pos + 1 < cnt; UINT8(dest)[pos >> 1] = lookup[src[pos]] | (uint8_t) (lookup[src[pos + 1]] << 4), pos += 2);
    if (cnt & 1)
    {
        size_t hcnt = cnt >> 1;
        UINT8(dest)[hcnt] = lookup[src[cnt - 1]] | (UINT8(dest)[hcnt] & 0xf0);
    }
}

void iupacDecode(iupacStorage *src, uint8_t *dest, size_t cnt, iupacDecodeLookup *lookup)
{
    for (size_t pos = 0; pos + 1 < cnt; pos += 2)
    {
        uint8_t tmp = UINT8(src)[pos >> 1];
        dest[pos] = lookup[tmp & 15];
        dest[pos + 1] = lookup[tmp >> 4];
    }
    if (cnt & 1) dest[cnt - 1] = lookup[UINT8(src)[cnt >> 1] & 15];
}
