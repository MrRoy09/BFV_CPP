#include "include/ntt.h"
#include <algorithm>
#include <vector>
#include <iostream>

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
    int n = a.size();
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
    for (int i = 0; i < n; i++)
    {
        a[i] = a[i] % q;
        if (a[i] < 0)
        {
            a[i]+=q;
        }
    }
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

vector<ll> index_reverse(const vector<ll> &a, int r) {
    int n = a.size();
    vector<ll> b(n);
    
    auto int_reverse = [](int a, int n) -> int {
        int result = 0;
        for (int i = 0; i < n; i++) {
            result = (result << 1) | (a & 1);
            a >>= 1;
        }
        return result;
    };
    
    for (int i = 0; i < n; i++) {
        int rev_idx = int_reverse(i, r);
        b[rev_idx] = a[i];
    }
    return b;
}

vector<ll> ntt_with_table(const vector<ll> &A, const vector<ll> &W_table, ll q) {
    int n = A.size();
    vector<ll> B = A;
    
    int v = 0;
    int temp = n;
    while (temp > 1) {
        temp >>= 1;
        v++;
    }
    
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < (1 << i); j++) {
            for (int k = 0; k < (1 << (v - i - 1)); k++) {
                int s = j * (1 << (v - i)) + k;
                int t = s + (1 << (v - i - 1));
                
                ll w = W_table[((1 << i) * k) % n];
                
                ll as_temp = B[s];
                ll at_temp = B[t];
                
                B[s] = (as_temp + at_temp) % q;
                B[t] = ((as_temp - at_temp) * w) % q;
                if (B[t] < 0) B[t] += q;
            }
        }
    }
    
    B = index_reverse(B, v);
    return B;
}

vector<ll> intt_with_table(const vector<ll> &A, const vector<ll> &W_table, ll q) {
    int n = A.size();
    vector<ll> B = A;
    
    int v = 0;
    int temp = n;
    while (temp > 1) {
        temp >>= 1;
        v++;
    }
    
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < (1 << i); j++) {
            for (int k = 0; k < (1 << (v - i - 1)); k++) {
                int s = j * (1 << (v - i)) + k;
                int t = s + (1 << (v - i - 1));
                
                ll w = W_table[((1 << i) * k) % n];
                
                ll as_temp = B[s];
                ll at_temp = B[t];
                
                B[s] = (as_temp + at_temp) % q;
                B[t] = ((as_temp - at_temp) * w) % q;
                if (B[t] < 0) B[t] += q;
            }
        }
    }
    
    B = index_reverse(B, v);
    
    // Apply modular inverse of n
    ll inv_n = mod_pow(n, q - 2, q);
    for (ll &x : B) {
        x = (x * inv_n) % q;
    }
    
    return B;
}