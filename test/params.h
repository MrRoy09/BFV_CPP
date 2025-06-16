#pragma once

#include <cstdint>

typedef int64_t ll;

ll n = 1024;        // Polynomial degree
ll q = 132120577; // Ciphertext modulus
ll root = 73993;      // Primitive root for NTT
ll t = 16;        // Plaintext modulus
ll mu = 0;        // Mean for Gaussian noise
ll sigma = 0;     // Stddev for Gaussian noise (zero for comparison)
ll T = 16;