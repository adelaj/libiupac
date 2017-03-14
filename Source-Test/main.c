#ifdef _MSC_VER
#   define _CRT_SECURE_NO_WARNINGS
#   pragma warning(disable : 4204; disable : 4221)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <libiupac.h>

#define countof(ARR) (sizeof (ARR) / sizeof (ARR)[0])

#ifdef max
#   undef max
#endif

#ifdef min
#   undef min
#endif

static inline size_t max(size_t a, size_t b) { return a > b ? a : b; }
static inline size_t min(size_t a, size_t b) { return a > b ? b : a; }

bool test01(size_t cnt, size_t rep)
{
    bool succ = 0;
    iupacEncodeLookup *tbl = iupacEncodeLookupCreate(iupacPatternStd());   
    iupacStorage *stor = iupacStorageCreate(cnt);
    uint8_t *text0 = malloc(cnt * sizeof *text0), *text1 = malloc(cnt * sizeof *text1);
    
    if (!(tbl && stor && text0 && text1)) goto error;

    for (size_t i = 0; i < rep; i++)
    {
        for (size_t j = 0; j < cnt; j++) text0[j] = iupacDecodeLookupStd()[(unsigned) rand() & 15];

        iupacEncode(text0, stor, cnt, tbl);
        iupacDecode(stor, text1, cnt, iupacDecodeLookupStd());

        if (memcmp(text0, text1, cnt)) goto error;
    }    

    succ = 1;

error:
    free(tbl);
    iupacStorageDispose(stor);
    free(text0);
    free(text1);

    return succ;
}

bool test02(size_t cnt, size_t rep)
{
    bool succ = 0;
    size_t ceil = (cnt + 31) >> 5;
    
    iupacEncodeLookup *tbl = iupacEncodeLookupCreate(iupacPatternStd());
    iupacStorage *stor0 = iupacStorageCreate(cnt), *stor1 = iupacStorageCreate(cnt);
    uint8_t *text0 = malloc(cnt * sizeof *text0), *text1 = malloc(cnt * sizeof *text1);
    uint32_t *mask0 = malloc(ceil * sizeof *mask0), *mask1 = malloc(ceil * sizeof *mask1);

    if (!(tbl && stor0 && stor1 && text0 && text1 && mask0 && mask1)) goto error;

    for (size_t i = 0; i < rep; i++)
    {
        for (size_t j = 0; j < cnt; j++) text0[j] = iupacDecodeLookupStd()[(unsigned) rand() & 15];
        iupacEncode(text0, stor0, cnt, tbl);

        // Test 1: full inclusion
        memset(mask0, 0, ceil * sizeof(uint32_t));
        
        for (size_t j = 0; j < cnt; j++) text1[j] = iupacDecodeLookupStd()[tbl[text0[j]] | ((unsigned) rand() & 15)], mask0[j >> 5] |= (uint32_t) 1 << (j & 31);
        iupacEncode(text1, stor1, cnt, tbl);

        if (!iupacTestInclusion(stor0, stor1, cnt)) goto error;
        
        iupacTestInclusionMask(stor0, stor1, mask1, cnt);
        if (memcmp(mask0, mask1, ceil * sizeof (uint32_t))) goto error;

        // Test 2: partial inclusion
        memset(mask0, 0, ceil * sizeof(uint32_t));
        
        for (size_t j = 0; j < cnt; j++)
        {
            uint8_t src = tbl[text0[j]], tmp = src & ((unsigned) rand() & 15);
            text1[j] = iupacDecodeLookupStd()[tmp];
            mask0[j >> 5] |= (~tmp & src) ? 0 : (uint32_t) 1 << (j & 31);
        }

        iupacEncode(text1, stor1, cnt, tbl);
                
        iupacTestInclusionMask(stor0, stor1, mask1, cnt);
        if (memcmp(mask0, mask1, ceil * sizeof(uint32_t))) goto error;
    }
    
    succ = 1;

error:
    free(tbl);
    iupacStorageDispose(stor0);
    iupacStorageDispose(stor1);
    free(text0);
    free(text1);
    free(mask0);
    free(mask1);

    return succ;
}

bool test03(size_t cnt, size_t rep)
{
    static uint8_t comp[] = {
        0b0000, 0b1000, 0b0100, 0b1100, 0b0010, 0b1010, 0b0110, 0b1110,
        0b0001, 0b1001, 0b0101, 0b1101, 0b0011, 0b1011, 0b0111, 0b1111       
    };
    
    bool succ = 0;
    
    iupacEncodeLookup *tbl = iupacEncodeLookupCreate(iupacPatternStd());
    iupacStorage *stor0 = iupacStorageCreate(cnt), *stor1 = iupacStorageCreate(cnt);
    uint8_t *text0 = malloc(cnt * sizeof *text0), *text1 = malloc(cnt * sizeof *text1);

    if (!(tbl && stor0 && stor1 && text0 && text1)) goto error;

    for (size_t i = 0; i < rep; i++)
    {
        // Test 1: full reverse complement
        for (size_t j = 0; j < cnt; j++)
        {
            uint8_t tmp = (unsigned) rand() & 15;                

            text0[j] = iupacDecodeLookupStd()[tmp];
            text1[cnt - j - 1] = iupacDecodeLookupStd()[comp[tmp]];
        }

        iupacEncode(text0, stor0, cnt, tbl);

        iupacReverseComplement(stor1, stor0, cnt);
        iupacDecode(stor1, text0, cnt, iupacDecodeLookupStd());
        
        if (memcmp(text0, text1, cnt)) goto error;

        // Test 2: partial reverse complement
        size_t part = ((size_t) rand() + (size_t) rand() * (size_t) RAND_MAX) % (cnt + 1);

        for (size_t j = 0; j < part; j++)
        {
            uint8_t tmp = (unsigned) rand() & 15;

            text0[j] = iupacDecodeLookupStd()[tmp];
            text1[part - j - 1] = iupacDecodeLookupStd()[comp[tmp]];
        }

        for (size_t j = part; j < cnt; j++) text0[j] = text1[j] = iupacDecodeLookupStd()[(unsigned) rand() & 15];

        iupacEncode(text0, stor0, cnt, tbl);
        iupacEncode(text1, stor1, cnt, tbl);

        iupacReverseComplement(stor1, stor0, part);
        iupacDecode(stor1, text0, cnt, iupacDecodeLookupStd());

        if (memcmp(text0, text1, cnt)) goto error;
    }
    
    succ = 1;

error:
    free(tbl);
    iupacStorageDispose(stor0);
    iupacStorageDispose(stor1);
    free(text0);
    free(text1);

    return succ;
}

int main()
{
    struct { bool(*callback)(size_t, size_t); size_t cntmin, cntmax, cntstp, repmul, repmax; const char *descr; } tests[] = {
        { .callback = test01, .cntmin = 0, .cntmax = 1000, .cntstp = 1, .repmul = 10, .repmax = 1000, .descr = "encoding/decoding" },
        { .callback = test02, .cntmin = 0, .cntmax = 1000, .cntstp = 1, .repmul = 10, .repmax = 1000, .descr = "fine/coarse inclusion" },
        { .callback = test03, .cntmin = 0, .cntmax = 1000, .cntstp = 1, .repmul = 10, .repmax = 1000, .descr = "reverse complement" }
    };
    
    srand((unsigned) time(NULL));

    for (size_t i = 0; i < countof(tests); i++)
    {
        bool res = 1;
        for (size_t j = tests[i].cntmin; res && j <= tests[i].cntmax; j += tests[i].cntstp) res &= tests[i].callback(j, min(j * tests[i].repmul, tests[i].repmax));
        printf("Test for %s %s.\n", tests[i].descr, res ? "succeeded" : "failed");        
    }
    return 0;
}