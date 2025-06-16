#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "ntt.h"

typedef int64_t ll;

class Poly
{
public:
    ll n;
    ll q;
    ll root;
    std::vector<ll> np; 
    std::vector<ll> F;
    bool inNTT;

    Poly() : n(0), q(0), root(0), np(4, 0), F(n, 0), inNTT(false) {}

    Poly(ll n, ll q, ll root, std::vector<ll> np = {0, 0, 0, 0})
        : n(n), q(q), root(root), np(np), F(n, 0), inNTT(false) {}

    void randomize(ll B, bool domain = false, int type = 0, double mu = 0.0, double sigma = 0.0)
    {
        static std::mt19937 gen(42); // Fixed seed for reproducible results

        if (type == 0)
        {
            ll bound = B / 2;
            std::uniform_int_distribution<long long> dist(-static_cast<long long>(bound), static_cast<long long>(bound));
            for (ll i = 0; i < n; ++i)
            {
                long long val = dist(gen);
                F[i] = ((val % q) + q) % q;
            }
        }
        else
        {
            std::normal_distribution<double> dist(mu, sigma);
            for (ll i = 0; i < n; ++i)
            {
                long long val = static_cast<long long>(std::round(dist(gen)));
                F[i] = ((val % q) + q) % q;
            }
        }

        inNTT = domain;
    }

    Poly operator+(const Poly &other) const
    {
        if (n != other.n || q != other.q)
            throw std::invalid_argument("Polynomial addition: sizes or moduli do not match");
        if (inNTT != other.inNTT)
            throw std::invalid_argument("Polynomial addition: Not in same domain");

        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = ((F[i] + other.F[i]) % q + q) % q;

        result.inNTT = inNTT;
        return result;
    }

    Poly operator-(const Poly &other) const
    {
        if (n != other.n || q != other.q)
            throw std::invalid_argument("Polynomial subtraction: sizes or moduli do not match");
        if (inNTT != other.inNTT)
            throw std::invalid_argument("Polynomial subtraction: Not in same domain");

        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = ((F[i] - other.F[i]) % q + q) % q;

        result.inNTT = inNTT;
        return result;
    }

    Poly operator*(const Poly &other) const
    {
        if (n != other.n || q != other.q)
            throw std::invalid_argument("Polynomial multiplication: sizes or moduli do not match");
        if (inNTT != other.inNTT)
            throw std::invalid_argument("Polynomial multiplication: Not in same domain");

        Poly result(n, q, root, np);
        result.inNTT = inNTT;

        if (inNTT)
        {
            for (ll i = 0; i < n; i++)
                result.F[i] = (F[i] * other.F[i]) % q;
        }
        else
        {
            if (np.size() >= 4 * n) {
                // Extract tables from np
                std::vector<ll> w_table(np.begin(), np.begin() + n);
                std::vector<ll> wv_table(np.begin() + n, np.begin() + 2*n);
                std::vector<ll> psi_table(np.begin() + 2*n, np.begin() + 3*n);
                std::vector<ll> psiv_table(np.begin() + 3*n, np.begin() + 4*n);
                
                // x1=self*psi, x2=b*psi
                std::vector<ll> s_p(n), b_p(n);
                for (ll i = 0; i < n; i++) {
                    s_p[i] = (F[i] * psi_table[i]) % q;
                    b_p[i] = (other.F[i] * psi_table[i]) % q;
                }
                
                // x1n = NTT(x1,w), x2n = NTT(x2,w)
                std::vector<ll> s_n = ntt_with_table(s_p, w_table, q);
                std::vector<ll> b_n = ntt_with_table(b_p, w_table, q);
                
                // x3n = x1n*x2n
                std::vector<ll> sb_n(n);
                for (ll i = 0; i < n; i++) {
                    sb_n[i] = (s_n[i] * b_n[i]) % q;
                }
                
                // x3 = INTT(x3n,w_inv)
                std::vector<ll> sb_p = intt_with_table(sb_n, wv_table, q);
                
                // c = x3*psi_inv
                for (ll i = 0; i < n; i++) {
                    result.F[i] = (sb_p[i] * psiv_table[i]) % q;
                }
            } else {
                // Fallback to original method
                result.F = multiply_ntt(F, other.F, q, root);
            }
        }

        return result;
    }

    Poly operator%(const ll base) const
    {
        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
        {
            result.F[i] = ((F[i] % base));
            if (result.F[i] < 0)
                result.F[i] += base;
        }

        result.inNTT = inNTT;
        return result;
    }

    Poly pround()
    {
        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = static_cast<ll>(round((long long)F[i]));

        result.inNTT = inNTT;
        return result;
    }

    bool operator==(const Poly &other) const
    {
        if (n != other.n || q != other.q || inNTT != other.inNTT)
            return false;
        for (ll i = 0; i < n; i++)
        {
            if (F[i] != other.F[i])
                return false;
        }
        return true;
    }

    Poly operator-() const
    {
        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = ((-F[i]) % q + q) % q;

        result.inNTT = inNTT;
        return result;
    }
};