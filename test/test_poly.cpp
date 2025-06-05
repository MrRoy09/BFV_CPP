#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/polynomials.hpp"
#include <vector>
#include <cstdint>

TEST_CASE("Poly randomize function", "[poly][random]")
{
    ll n = 64;
    ll q = 998244353;
    ll root = 3;
    std::vector<ll> ntt_params = {3, 8193, 5, 2457}; // Example values
    Poly poly(n, q, root, ntt_params);

    std::vector<ll> coff;
    bool randomized = false;

    SECTION("Uniform randomization")
    {
        poly.randomize(100, true, 0); // B = 100, domain=true, type=0
        REQUIRE(poly.F.size() == n);
        for (ll coeff : poly.F)
        {
            REQUIRE(coeff >= 0);
            REQUIRE(coeff < q);
        }
        REQUIRE(poly.inNTT == true);
        coff = poly.F;

        poly.randomize(100, true, 0); // Run again to check randomness
        REQUIRE(poly.F.size() == n);
        for (ll i = 0; i < n; i++)
        {
            if (poly.F[i] != coff[i])
                randomized = true;
        }
        REQUIRE(randomized == true);
    }

    SECTION("Gaussian randomization")
    {
        poly.randomize(0, false, 1, 0.0, 5.0); // type=1, mu=0, sigma=5
        REQUIRE(poly.F.size() == n);
        for (ll coeff : poly.F)
        {
            REQUIRE(coeff >= 0);
            REQUIRE(coeff < q);
        }
        REQUIRE(poly.inNTT == false);
        coff = poly.F;
        randomized = false;

        poly.randomize(100, true, 0); // Now test another uniform call
        REQUIRE(poly.F.size() == n);
        for (ll i = 0; i < n; i++)
        {
            if (poly.F[i] != coff[i])
                randomized = true;
        }
        REQUIRE(randomized == true);
    }
}
