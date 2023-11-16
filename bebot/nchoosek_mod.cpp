#include "../include/nchoosek_mod.h"

unsigned long long nchoosek_mod(int N, int k) {
    unsigned long long binom = 1;
    for (int j = 1; j <= k; ++j) {
        binom *= (N - (k - j));
        binom /= j;
    }
    return binom;
}
