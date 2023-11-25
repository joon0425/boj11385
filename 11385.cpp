#include <bits/stdc++.h>
using namespace std;
#define fastio ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0)
typedef long long ll;
typedef vector<ll> poly;
const int N = 1 << 21;
const int M1 = (483 << 21) + 1, R1 = 62;
const int M2 = (479 << 21) + 1, R2 = 62;

ll modpow(ll a,ll b,ll m){ // return a^b mod m
    ll res=1;
    while(b){
        if(b&1)res=res*a%m;
        b>>=1;a=a*a%m;
    }
    return res;
}

using uint = unsigned int;
alignas(32) uint a[N];
alignas(32) uint b[N];

alignas(32) uint wdit[2][N];
alignas(32) uint wpdit[2][N];
template <ll M, ll R, bool dit = 0> void ntt(uint a[N]) {
    a = (uint*) __builtin_assume_aligned(a, 32);
    auto w = (uint*) __builtin_assume_aligned(wdit[dit], 32);
    auto wp = (uint*) __builtin_assume_aligned(wpdit[dit], 32);
    w[0] = w[1] = 1;
    wp[0] = wp[1] = (1LL << 32) / M;
    int n = N;
    //static vector<uint> w(2, 1), wp(2, (1LL << 32) / M); int n = a.size();
    for (static int k = 2; k < n; k *= 2) {
        w[k] = 1;
        ll c = modpow(R, M/2/k * (dit ? M-2 : 1), M);
        for (int i = k+1; i < 2*k; i++) w[i] = w[i-1]*c % M;
        for (int i = k; i < 2*k; i++) wp[i] = ((ll) w[i] << 32) / M;
    }
    for (int t = 1; t < n; t *= 2) {
        int k = (dit ? t : n/2/t);
        for (int i = 0; i < n; i += k*2) {
            for (int j = 0; j < k; j++) {
                uint X = a[i+j], Y = a[i+j+k], W = w[j+k], Wp = wp[j+k], T, Q;
                if (!dit) {
                    if ((a[i+j] += Y) >= 2*M) a[i+j] -= 2*M;
                    T = X - Y + 2*M;
                    Q = ((unsigned long long) Wp * T) >> 32;
                    a[i+j+k] = W*T - Q*M;
                } else {
                    if (X >= 2*M) X -= 2*M;
                    Q = ((unsigned long long) Wp * Y) >> 32;
                    T = W*Y - Q*M;
                    a[i+j] = X + T;
                    a[i+j+k] = X - T + 2*M;
                }
            }
        }
    }
    for (int i = 0; i < N; i++) a[i] %= M;
}

template <ll M = (119<<23)+1, ll R = 62>
void conv(uint a[N], uint b[N]) {
    int len = N - 1, n = 1 << (32 - __builtin_clz(len));
    ll inv = modpow(n, M-2, M);
    ntt<M,R>(a); ntt<M,R>(b);
    for(int i=0; i<n; i++) a[i] = 1LL*a[i]*b[i] % M * inv % M;
    ntt<M,R,1>(a);
}

ll ext_gcd(ll a, ll b, ll &x, ll &y) { // ax + by = gcd(a, b)
    ll g = a; x = 1, y = 0;
    if (b) g = ext_gcd(b, a % b, y, x), y -= a / b * x;
    return g;
}
ll CRT(const vector<ll> a, const vector<ll> m){ // a_i mod b_i
    int sz = a.size();
    vector<ll> rmn(sz), lm(sz, 1);
    ll ans = 0, M = 1;
    for(int i=0; i<sz; i++){
        ll k = a[i] - rmn[i]; k %= m[i];
        if(k < 0) k += m[i];
        ll x, y;
        ext_gcd(lm[i], m[i], x, y);
        k *= x; k %= m[i];
        if(k < 0) k += m[i];
        ans += k*M;
        for(int t=i+1; t<sz; t++){
            rmn[t] += lm[t] * k; rmn[t] %= m[t];
            lm[t] *= m[i]; lm[t] %= m[t];
        }
        M *= m[i];
    }
    return ans;
}

int main(){
    fastio;
    ll n,m,r=0;cin>>n>>m;
    vector<ll> u(n+1),v(m+1);
    for(int i=0;i<=n;i++)cin>>u[i];
    for(int q=0;q<=m;q++)cin>>v[q];
    copy(begin(u), end(u), a);
    copy(begin(v), end(v), b);
    conv<M1, R1>(a,b);
    vector<uint> ans1(n+m+1);
    copy(a, a+n+m+1, begin(ans1));
    memset(a, 0, sizeof a);
    memset(b, 0, sizeof b);
    copy(begin(u), end(u), a);
    copy(begin(v), end(v), b);
    conv<M2, R2>(a,b);
    vector<uint> ans2(n+m+1);
    copy(a, a+n+m+1, begin(ans2));
    for(int i=0;i<=n+m;i++)r^=CRT({ans1[i], ans2[i]},{M1, M2});
    cout<<r;
}