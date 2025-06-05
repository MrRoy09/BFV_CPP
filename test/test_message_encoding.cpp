#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/BFV.hpp"
#include "../include/polynomials.hpp"

TEST_CASE("BFV IntEncode and IntDecode", "[BFV][EncodeDecode]")
{
    int n = 64;
    int q = 998244353;
    int t = 17;
    int root = 3;
    int mu = 0;
    int sigma = 2;

    std::vector<ll> ntt_params = {0};
    BFV bfv(n, q, t, mu, sigma, root);

    bfv.sk_gen();
    bfv.pk_gen();

    SECTION("Positive integers")
    {
        for (ll i = 0; i <= 100; i++)
        {
            Poly encoded = bfv.IntEncode(i);
            ll decoded = bfv.IntDecode(encoded);
            REQUIRE(decoded == i);
        }
    }

    SECTION("Negative integers")
    {
        for (ll i = -100; i < 0; i++)
        {
            Poly encoded = bfv.IntEncode(i);
            ll decoded = bfv.IntDecode(encoded);
            REQUIRE(decoded == i);
        }
    }

    SECTION("Zero encoding")
    {
        Poly encoded = bfv.IntEncode(0);
        ll decoded = bfv.IntDecode(encoded);
        REQUIRE(decoded == 0);
    }

    SECTION("BIG NUMBER")
    {
        Poly encoded = bfv.IntEncode((72803712312323234));
        ll decoded = bfv.IntDecode(encoded);
        INFO(decoded);
        REQUIRE(72803712312323234==decoded);
    }

    SECTION("COMPLETE PIPELINE")
    {
        ll message = 0xfffffffffffffff1;
        Poly encoded = bfv.IntEncode(message);
        auto ct = bfv.Encryption(encoded);
        auto pt = bfv.Decryption(ct);
        ll decoded = bfv.IntDecode(pt);
        INFO(decoded);
        REQUIRE(decoded==message);

    }
}
