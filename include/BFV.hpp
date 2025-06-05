#pragma once

#include "polynomials.hpp"
#include <vector>
#include <cstdint>
#include <cmath>

typedef long long ll;

class BFV
{
public:
    ll n;
    ll q;
    ll t;
    ll mu;
    ll sigma;
    ll root;
    ll T;
    ll l;

    std::vector<ll> np;
    Poly sk;
    std::pair<Poly, Poly> pk;

    std::vector<std::pair<Poly, Poly>> rlk1;
    std::vector<std::pair<Poly, Poly>> rlk2;

    BFV(ll n, ll q, ll t, ll mu, ll sigma, ll root)
        : n(n), q(q), t(t), mu(mu), sigma(sigma), root(root), T(0), l(0) {}

    std::vector<long long> RefPolMul(const std::vector<long long> &A, const std::vector<long long> &B)
    {
        size_t n = A.size();
        std::vector<long long> C(2 * n, 0);
        std::vector<long long> D(n, 0);

        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                __int128_t product = static_cast<__int128_t>(A[i]) * B[j];
                C[i + j] = (C[i + j] + static_cast<long long>(product % q)) % q;
                if (C[i + j] < 0)
                    C[i + j] += q;
            }
        }

        for (size_t i = 0; i < n; ++i)
        {
            D[i] = (C[i] - C[i + n]) % q;
            if (D[i] < 0)
                D[i] += q;
        }

        return D;
    }

    void sk_gen()
    {
        Poly s(n, q, root, np);
        s.randomize(2);
        sk = s;
    }

    void pk_gen()
    {
        Poly a(n, q, root, np);
        Poly e(n, q, root, np);
        a.randomize(q);
        e.randomize(0, false, 1, mu, sigma);
        Poly pk0 = -(a * sk + e);
        Poly pk1 = a;
        pk = {pk0, pk1};
    }
    std::pair<Poly, Poly> Encryption(Poly &encoded_message)
    {
        long long delta = q / t; 

        Poly u(n, q, root, np); 
        Poly e1(n, q, root, np);
        Poly e2(n, q, root, np);

        u.randomize(2);
        e1.randomize(0, false, 1, mu, sigma);
        e2.randomize(0, false, 1, mu, sigma);

        Poly md(n, q, root, np);
        for (long long i = 0; i < n; i++)
        {
            md.F[i] = (delta * encoded_message.F[i]) % q;
            if (md.F[i] < 0)
                md.F[i] += q; 
        }

        Poly c0 = pk.first * u + e1;
        c0 = c0 + md;
        Poly c1 = pk.second * u + e2;

        return {c0, c1};
    }

    Poly Decryption(std::pair<Poly, Poly> encrypted_message)
    {
        Poly m = encrypted_message.second * sk + encrypted_message.first;

        for (long long i = 0; i < n; i++)
        {
            double scaled = (static_cast<double>(t) * m.F[i]) / q;
            long long rounded = static_cast<long long>(std::round(scaled));
            m.F[i] = rounded % t;
            if (m.F[i] < 0)
                m.F[i] += t;
        }

        Poly mr(n, t, root, np); 
        mr.F = m.F;
        mr.inNTT = m.inNTT;

        return mr;
    }

    Poly IntEncode(ll m)
    {
        Poly mr(n, t, root, np);
        if (m > 0)
        {
            ll mt = m;
            for (ll i = 0; i < n; i++)
            {
                mr.F[i] = mt % 2;
                mt = mt / 2;
            }
        }
        else if (m < 0)
        {
            ll mt = -m;
            for (ll i = 0; i < n; i++)
            {
                mr.F[i] = (t - (mt % 2)) % t;
                mt = mt / 2;
            }
        }
        return mr;
    }

    ll IntDecode(const Poly &m)
    {
        ll mr = 0;
        ll thr_ = (t == 2) ? 2 : ((t + 1) >> 1);
        for (ll i = 0; i < n; i++)
        {
            ll c = m.F[i];
            ll c_ = (c >= thr_) ? -(t - c) : c;
            mr += c_ * (1LL << i);
        }
        return mr;
    }


    std::pair<Poly, Poly> HomomorphicAddition(const std::pair<Poly, Poly> &ct0,
                                              const std::pair<Poly, Poly> &ct1)
    {
        Poly c0 = ct0.first + ct1.first;
        Poly c1 = ct0.second + ct1.second;
        return {c0, c1};
    }

    std::pair<Poly, Poly> HomomorphicSubtraction(const std::pair<Poly, Poly> &ct0,
                                                 const std::pair<Poly, Poly> &ct1)
    {
        Poly c0 = ct0.first - ct1.first;
        Poly c1 = ct0.second - ct1.second;
        return {c0, c1};
    }
};