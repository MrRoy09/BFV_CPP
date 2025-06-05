#include "include/ntt.h"
#include <algorithm>
#include <vector>
#include <iostream>

typedef long long ll;
using namespace std;



ll mod_pow(ll base, ll exp, ll mod)
{
    ll result = 1;
    while (exp > 0)
    {
        if (exp % 2)
            result = (result * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

void ntt(vector<ll> &a, ll q, ll root, bool invert)
{
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++)
    {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1)
    {
        ll wlen = mod_pow(root, (q - 1) / len, q);
        if (invert)
            wlen = mod_pow(wlen, q - 2, q);

        for (int i = 0; i < n; i += len)
        {
            ll w = 1;
            for (int j = 0; j < len / 2; j++)
            {
                ll u = a[i + j];
                ll v = (a[i + j + len / 2] * w) % q;

                a[i + j] = (u + v < q) ? (u + v) : (u + v - q);
                a[i + j + len / 2] = (u - v >= 0) ? (u - v) : (u - v + q);

                w = (w * wlen) % q;
            }
        }
    }

    if (invert)
    {
        ll inv_n = mod_pow(n, q - 2, q);
        for (ll &x : a)
            x = (x * inv_n) % q;
    }
}

vector<ll> multiply_ntt(vector<ll> a, vector<ll> b, ll q, ll root)
{
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    a.resize(n);
    b.resize(n);

    ntt(a, q, root, false);
    ntt(b, q, root, false);
    for (int i = 0; i < n; i++)
    {
        a[i] = (a[i] * b[i]) % q;
        if (a[i] < 0)
        {
            a[i] += q;
        }
    }

    ntt(a, q, root, true);
    return a;
}

vector<ll> multiply_naive(vector<ll> a, vector<ll> b, ll q)
{
    vector<ll> res(a.size() + b.size() - 1);
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < b.size(); j++)
        {
            res[i + j] = (res[i + j] + a[i] * b[j]) % q;
            if (res[i + j] < 0)
                res[i + j] += q;
        }
    return res;
}
