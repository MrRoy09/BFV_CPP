#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/polynomials.hpp"
#include "../include/BFV.hpp"
#include <iostream>
#include <vector>
#include <cmath>

TEST_CASE("BFV Encryption and Decryption", "[BFV]") {
    ll n = 64;                // Polynomial degree
    ll q = 998244353;         // Ciphertext modulus
    ll root = 3;              // Primitive root for NTT
    ll t = 17;                // Plaintext modulus
    ll mu = 0;                // Mean for Gaussian noise
    ll sigma = 2;             // Stddev for Gaussian noise

    std::vector<ll> np = {0}; // Placeholder for NTT precomputation
    BFV bfv(n, q, t, mu, sigma, root);

    bfv.sk_gen();
    bfv.pk_gen();

    Poly message(n, q, root, np);
    for (ll i = 0; i < 32; i++) {
        message.F[i] = i*i*i*i*i;
    }

    std::pair<Poly, Poly> ciphertext = bfv.Encryption(message);
    Poly decrypted = bfv.Decryption(ciphertext);

    for (ll i = 0; i < n; i++) {
        INFO("Mismatch at index " << i << ": expected " << message.F[i] % t
             << ", got " << decrypted.F[i] % t);
        REQUIRE(decrypted.F[i] % t == message.F[i] % t);
    }
}
