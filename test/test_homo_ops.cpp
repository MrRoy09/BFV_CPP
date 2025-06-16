#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/BFV.hpp"
#include "params.h"

TEST_CASE("BFV Homomorphic Operations", "[bfv][homomorphic]")
{

    BFV bfv(n, q, t, mu, sigma, root);

    bfv.sk_gen();
    bfv.pk_gen();

    bfv.EvalKeyGen(T);

    ll m1 = 50;
    ll m2 = 3;

    Poly pt1 = bfv.IntEncode(m1);
    Poly pt2 = bfv.IntEncode(m2);

    auto ct1 = bfv.Encryption(pt1);
    auto ct2 = bfv.Encryption(pt2);

    SECTION("Homomorphic Addition")
    {
        auto ct_add = bfv.HomomorphicAddition(ct1, ct2);
        Poly pt_dec = bfv.Decryption(ct_add);
        ll result = bfv.IntDecode(pt_dec);
        ll expected = (m1 + m2);
        INFO("Addition Test Failed");
        REQUIRE(result == expected);
    }

    SECTION("Homomorphic Subtraction")
    {
        auto ct_sub = bfv.HomomorphicSubtraction(ct1, ct2);
        Poly pt_dec = bfv.Decryption(ct_sub);

        std::cout << "Subtraction coeffs: ";

        ll result = (bfv.IntDecode(pt_dec));
        ll expected = (m1 - m2);

        INFO("Subtraction Test Failed");
        REQUIRE(result == expected);
    }

    SECTION("Homomorphic Multiplication")
    {

        auto ct_mul = bfv.HomomorphicMultiplication(ct1, ct2);

        Poly sk2 = bfv.sk * bfv.sk;
        Poly m_before_scaling = ct_mul[0] + (ct_mul[1] * bfv.sk) + (ct_mul[2] * sk2);

        auto ct_rel = bfv.Relinearization(ct_mul);
        Poly pt_dec2 = bfv.DecryptionV2(ct_mul);
        Poly pt_dec = bfv.Decryption(ct_rel);

        ll result = (((bfv.IntDecode(pt_dec))));
        ll result2 = bfv.IntDecode(pt_dec2);
        ll expected = (m1 * m2);
        INFO("Multiplication Test Failed");
        REQUIRE(result2 == expected);
        REQUIRE(result == expected);
    }
}
