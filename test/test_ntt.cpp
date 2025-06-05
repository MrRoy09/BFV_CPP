#define CATCH_CONFIG_MAIN
#include "../include/catch_amalgamated.hpp"
#include "../include/ntt.h"

const int q = 998244353; // Modulus
const int root = 3;      // Primitive root for q

TEST_CASE("Basic polynomial multiplication with NTT", "[ntt]")
{
    vector<long long> A = {1, 2, 3}; // A(x) = 1 + 2x + 3x^2
    vector<long long> B = {4, 5, 6}; // B(x) = 4 + 5x + 6x^2

    auto expected = multiply_naive(A, B,q);
    auto result = multiply_ntt(A, B, q, root);

    expected.resize(result.size());
    REQUIRE(result == expected);
}

TEST_CASE("NTT handles zeros", "[ntt]")
{
    vector<long long> A = {0, 0, 0};
    vector<long long> B = {1, 2, 3};

    auto result = multiply_ntt(A, B, q, root);
    for (int x : result)
    {
        REQUIRE(x == 0);
    }
}

TEST_CASE("Varying Lengths of Coeff", "[ntt]")
{
    vector<long long> A = {99, 55, 2131, 3123, 12324};
    vector<long long> B = {490, 5131232, 6231231, 123123123, 1231233, 4134324, 2312132};

    auto expected = multiply_naive(A, B,q);
    auto result = multiply_ntt(A, B, q, root);

    expected.resize(result.size());
    for (int i = 0; i < expected.size(); i++)
    {
        REQUIRE(result[i]==expected[i]);
    }
    REQUIRE(result == expected);
}
