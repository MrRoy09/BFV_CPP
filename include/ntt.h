#pragma once

#include <vector>
#include <cstdint>
using namespace std;

typedef int64_t ll;

const ll MOD = 998244353;
const ll ROOT = 3;

ll mod_pow(ll base, ll exp, ll mod);

void ntt(vector<ll> &a, ll q, ll root, bool invert);

vector<ll> multiply_ntt(vector<ll> a, vector<ll> b, ll q, ll root);

vector<ll> multiply_naive(vector<ll> a, vector<ll> b, ll q);

vector<ll> ntt_with_table(const vector<ll> &A, const vector<ll> &W_table, ll q);
vector<ll> intt_with_table(const vector<ll> &A, const vector<ll> &W_table, ll q);