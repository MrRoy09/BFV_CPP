#pragma once

#include "polynomials.hpp"
#include <vector>
#include <cstdint>
#include <cmath>
#include <functional>

typedef int64_t ll;

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
        : n(n), q(q), t(t), mu(mu), sigma(sigma), root(root), T(0), l(0) 
    {
        // Generate polynomial arithmetic tables
        // psi = root, w = root^2
        ll psi = root;
        ll w = (root * root) % q;
        
        // Compute modular inverse using extended GCD
        auto modinv = [](ll a, ll m) -> ll {
            std::function<std::tuple<ll, ll, ll>(ll, ll)> egcd = [&](ll a, ll b) -> std::tuple<ll, ll, ll> {
                if (a == 0) return {b, 0, 1};
                auto [g, y, x] = egcd(b % a, a);
                return {g, x - (b / a) * y, y};
            };
            auto [g, x, y] = egcd(a, m);
            if (g != 1) throw std::runtime_error("Modular inverse does not exist");
            return (x % m + m) % m;
        };
        
        ll psiv = modinv(psi, q);
        ll wv = modinv(w, q);
        
        // Generate tables
        std::vector<ll> w_table(n, 1);
        std::vector<ll> wv_table(n, 1);
        std::vector<ll> psi_table(n, 1);
        std::vector<ll> psiv_table(n, 1);
        
        for (int i = 1; i < n; i++) {
            w_table[i] = (w_table[i-1] * w) % q;
            wv_table[i] = (wv_table[i-1] * wv) % q;
            psi_table[i] = (psi_table[i-1] * psi) % q;
            psiv_table[i] = (psiv_table[i-1] * psiv) % q;
        }
        
        // Store in np vector: [w_table, wv_table, psi_table, psiv_table]
        // Note: We'll need to modify this to store the actual tables
        np.resize(4 * n);
        for (int i = 0; i < n; i++) {
            np[i] = w_table[i];
            np[n + i] = wv_table[i];
            np[2*n + i] = psi_table[i];
            np[3*n + i] = psiv_table[i];
        }
    }

    void EvalKeyGen(ll T_param)
    {
        this->T = T_param;
        this->l = (ll)std::log((long long)q) / std::log((long long)T);
        rlk1.clear();
        Poly sk2 = sk * sk;

        for (int i = 0; i < l + 1; i++)
        {
            Poly ai, ei;
            ai = Poly(n, q, root, np);
            ei = Poly(n, q, root, np);
            ai.randomize(q);
            ei.randomize(0, false, 1, mu, sigma);
            Poly Ts2 = Poly(n, q, root, np);
            for (int j = 0; j < sk2.F.size(); j++)
            {
                Ts2.F[j] = (((ll)pow((long long)T, i)) * sk2.F[j]) % q;
            }
            Poly rlki0 = -(ai * sk + ei) + Ts2;
            Poly rlki1 = ai;
            rlk1.push_back({rlki0, rlki1});
        }
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
            long long rounded = static_cast<long long>(std::llround(scaled));
            m.F[i] = rounded % t;
            if (m.F[i] < 0)
                m.F[i] += t;
        }

        Poly mr(n, t, root, np);
        mr.F = m.F;
        mr.inNTT = m.inNTT;

        return mr;
    }

    Poly DecryptionV2(const std::vector<Poly> &ct)
    {
        Poly sk2 = sk * sk;
        Poly m = ct[0] + (ct[1] * sk) + (ct[2] * sk2);


        for (size_t i = 0; i < n; i++)
        {
            long double scaled = ((long double)t * m.F[i]) / q;
            m.F[i] = ((ll)std::llround(scaled)) % t;
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
        Poly mr(n, t, 0, np);
        ll mt = std::abs((long long)m);
        for (int i = 0; i < n; i++)
        {
            if (m > 0)
            {
                mr.F[i] = (mt % 2);
            }
            else if (m < 0)
            {
                mr.F[i] = (t - (mt % 2)) % t;
            }
            mt /= 2;
        }
        return mr;
    }

    ll IntDecode(const Poly &m)
    {
        ll mr = 0;
        ll thr = (t == 2) ? 2 : ((t + 1) >> 1);
        ll power = 1;

        for (size_t i = 0; i < m.F.size(); ++i)
        {
            ll c = m.F[i];
            ll c_ = (c >= thr) ? -(t - c) : c;
            mr += c_ * power;
            power *= 2;
        }
        return mr;
    }

    std::pair<Poly, Poly> Relinearization(std::vector<Poly> &ct)
    {
        if (ct.size() != 3)
        {
            throw std::invalid_argument("Size of ciphertext must be 3 for relinearization");
        }
        Poly c0, c1, c2;
        c0 = ct[0];
        c1 = ct[1];
        c2 = ct[2];

        std::vector<Poly> c2i;
        Poly c2q = Poly(n, q, root, np);
        c2q.F = c2.F;
        for (int i = 0; i < l + 1; i++)
        {
            Poly c2r = Poly(n, q, root, np);
            for (int j = 0; j < n; j++)
            {
                ll qt = (ll)(c2q.F[j] / T);
                ll rt = (ll)(c2q.F[j] - qt * T);
                c2q.F[j] = qt;
                c2r.F[j] = rt;
            }
            c2i.push_back(c2r);
        }
        Poly c0r, c1r;
        c0r = Poly(n, q, root, np);
        c1r = Poly(n, q, root, np);
        c0r.F = c0.F;
        c1r.F = c1.F;

        for (int i = 0; i < l + 1; i++)
        {
            c0r = c0r + (rlk1[i].first * c2i[i]);
            c1r = c1r + (rlk1[i].second * c2i[i]);
        }
        return {c0r, c1r};
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

    std::vector<ll> RefPolyMulv2(const std::vector<ll> &poly_a, const std::vector<ll> &poly_b)
    {
        size_t degree = poly_a.size();

        std::vector<ll> product_coeffs(2 * degree, 0);
        std::vector<ll> reduced_poly(degree, 0);

        for (size_t i = 0; i < degree; i++)
        {
            for (size_t j = 0; j < degree; j++)
            {
                product_coeffs[i + j] = (product_coeffs[i + j] + poly_a[i] * poly_b[j]);
            }
        }

        for (size_t i = 0; i < degree; i++)
        {
            reduced_poly[i] = (product_coeffs[i] - product_coeffs[i + degree]);
        }

        return reduced_poly;
    }

    std::vector<Poly> HomomorphicMultiplication(const std::pair<Poly, Poly> &ct0, const std::pair<Poly, Poly> &ct1)
    {
        std::vector<ll> r0 = RefPolyMulv2(ct0.first.F, ct1.first.F);
        std::vector<ll> r1 = RefPolyMulv2(ct0.first.F, ct1.second.F);
        std::vector<ll> r2 = RefPolyMulv2(ct0.second.F, ct1.first.F);
        std::vector<ll> r3 = RefPolyMulv2(ct0.second.F, ct1.second.F);

        std::vector<ll> c0 = r0;
        std::vector<ll> c1(r1.size());
        std::vector<ll> c2 = r3;

        for (size_t i = 0; i < n; i++)
        {
            // Apply BFV scaling
            c1[i] = r1[i] + r2[i];
            
            double c0_scaled = ((double)t * r0[i]) / q;
            double c1_scaled = ((double)t * c1[i]) / q;
            double c2_scaled = ((double)t * r3[i]) / q;
            
            c0[i] = ((ll)round(c0_scaled)) % q;
            c1[i] = ((ll)round(c1_scaled)) % q;
            c2[i] = ((ll)round(c2_scaled)) % q;

            if (c0[i] < 0) c0[i] += q;
            if (c1[i] < 0) c1[i] += q;
            if (c2[i] < 0) c2[i] += q;
        }
        Poly r_poly0 = Poly(n, q, root, np);
        Poly r_poly1(n, q, root, np);
        Poly r_poly2(n, q, root, np);

        r_poly0.F = c0;
        r_poly1.F = c1;
        r_poly2.F = c2;

        return {r_poly0, r_poly1, r_poly2};
    }
};