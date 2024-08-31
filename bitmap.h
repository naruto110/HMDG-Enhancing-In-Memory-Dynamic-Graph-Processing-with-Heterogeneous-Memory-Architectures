#ifndef __BITMAP_H_
#define __BITMAP_H_

#include "global.h"

#include <stdint.h>
#include <stdio.h>

#define BITS_PER_WORD 32
#define BITMAP_SIZE 1  // 256 bits / 32 bits per word = 8 words

typedef struct {
    uint32_t word;
} Bitmap32;

void initializeBitmap(Bitmap32* bitmap) {
    bitmap->word = 0;
}

void setBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word |= (1U << bitOffset);
}

void clearBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word &= ~(1U << bitOffset);
}

bool testBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    return (bitmap->word & (1U << bitOffset)) != 0;
}

void andBit(Bitmap32* des_bitmap, Bitmap32* src_bitmap) {
    des_bitmap->word = des_bitmap->word | src_bitmap->word;
}

void print_bitmap(Bitmap32* bitmap) {
    for (int i = 0; i < 32; i++) {
        if (testBit(bitmap, i)) {
            cout << 1 << " ";
        } else {
            cout << 0 << " ";
        }
    }
    cout << endl;
}



#endif

