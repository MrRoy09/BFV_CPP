#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/polynomials.hpp"
#include "../include/BFV.hpp"
#include "params.h"
#include <iostream>
#include <vector>
#include <cmath>

TEST_CASE("BFV Encryption and Decryption", "[BFV]") {


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
        REQUIRE(decrypted.F[i] % t == message.F[i] % t);
    }
}
