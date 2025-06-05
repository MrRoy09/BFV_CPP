#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/BFV.hpp"

TEST_CASE("BFV Homomorphic Operations", "[bfv][homomorphic]")
{
    ll n = 64;
    ll q = 998244353;
    ll t = 17;
    ll mu = 0;
    ll sigma = 2;
    ll root = 3;
    ll T = 256;

    BFV bfv(n, q, t, mu, sigma, root);

    bfv.sk_gen();
    bfv.pk_gen();

    ll m1 = 10;
    ll m2 = 3;

    Poly pt1 = bfv.IntEncode(m1);
    Poly pt2 = bfv.IntEncode(m2);

    auto ct1 = bfv.Encryption(pt1);
    auto ct2 = bfv.Encryption(pt2);

    SECTION("Homomorphic Addition")
    {
        auto ct_add = bfv.HomomorphicAddition(ct1, ct2);
        Poly pt_dec = bfv.Decryption(ct_add);
        ll result = (((bfv.IntDecode(pt_dec)) % t) + t) % t;
        ll expected = (m1 + m2) % t;
        if (expected < 0)
            expected += t;

        INFO("Addition Test Failed");
        CAPTURE(m1, m2, result, expected);
        REQUIRE(result == expected);
    }

    SECTION("Homomorphic Subtraction")
    {
        auto ct_sub = bfv.HomomorphicSubtraction(ct1, ct2);
        Poly pt_dec = bfv.Decryption(ct_sub);
        ll result = (((bfv.IntDecode(pt_dec)) % t) + t) % t;
        ll expected = (m1 - m2) % t;
        if (expected < 0)
            expected += t;

        INFO("Subtraction Test Failed");
        CAPTURE(m1, m2, result, expected);
        REQUIRE(result == expected);
    }
}
