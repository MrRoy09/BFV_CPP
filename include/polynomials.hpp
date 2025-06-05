#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "ntt.h"

typedef long long ll;

class Poly
{
public:
    ll n;
    ll q;
    ll root;
    std::vector<ll> np; // Reserved for future root precomputations
    std::vector<ll> F;
    bool inNTT;

    Poly() : n(0), q(0), root(0), np(4, 0), F(n, 0), inNTT(false) {}

    Poly(ll n, ll q, ll root, std::vector<ll> np = {0, 0, 0, 0})
        : n(n), q(q), root(root), np(np), F(n, 0), inNTT(false) {}

    void randomize(ll B, bool domain = false, int type = 0, double mu = 0.0, double sigma = 0.0)
    {
        std::random_device rd;
        std::mt19937 gen(rd());

        if (type == 0)
        {
            ll bound = B / 2;
            std::uniform_int_distribution<ll> dist(-bound, bound - 1);
            for (ll i = 0; i < n; ++i)
            {
                ll val = dist(gen);
                F[i] = ((val % q) + q) % q;
            }
        }
        else
        {
            std::normal_distribution<double> dist(mu, sigma);
            for (ll i = 0; i < n; ++i)
            {
                ll val = static_cast<ll>(std::round(dist(gen)));
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
            result.F = multiply_ntt(F, other.F, q, root);
        }

        return result;
    }

    Poly operator%(const ll base) const
    {
        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = ((F[i] % base) + base) % base;

        result.inNTT = inNTT;
        return result;
    }

    Poly pround()
    {
        Poly result(n, q, root, np);
        for (ll i = 0; i < n; i++)
            result.F[i] = static_cast<ll>(round(F[i]));

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
