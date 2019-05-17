/*
 Shahadat Hossain
 I.C.T Department
 Comilla University
 Session: 2013-2014
 */
#include<bits/stdc++.h>

#include<cstdlib>
#include<cctype>
#include<fstream>
#include<iterator>

#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <bitset>
#include <queue>
#include <stack>
#include <cmath>
#include <ctime>
#include <list>
#include <set>
#include <map>

using namespace std;
/*
 int fx[]= {0, 0, 1, -1};
 int fy[]= {1, -1, 0, 0};
 
 ll dx[]= {-1, 1, 0, 0, 1, 1, -1, -1};
 ll dy[]= {0, 0, 1, -1, 1, -1, -1, 1};
 
 ll cx[] = {1,1,2,2,-1,-1,-2,-2};
 ll cy[] = {2,-2,1,-1,2,-2,1,-1};
 */
//#include <unordered_map>
//#include <ext/pb_ds/assoc_container.hpp> // Common file
//#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update

#define valid(x,y) x>=0 && y>=0 && x<n && y<n
#define all(x) x.begin(),x.end()
#define pb(x) push_back(x)
#define sz(x) (int)x.size()
#define int long long


#define ms(a,b) memset(a, b, sizeof(a))
#define eps 1e-5

#define D(x) cout << #x << " = " << x << endl
#define DD(x,y) cout << #x << " = " << x << "   " << #y << " = " << y << endl
#define DDD(x,y,z) cout << #x << " = " << x << "   " << #y << " = " << y << "   " << #z << " = " << z << endl
#define prntvctr(x) for(int i = 0; i < x.size(); i++) printf("%d%c", x[i], i + 1 == x.size() ? '\n' : ' ' );
#define pc(x) __builtin_popcount(x)

#define input freopen("/Users/shahadat/Desktop/input.txt", "r", stdin)
#define output freopen("/Users/shahadat/Desktop/output.txt", "w", stdout)
//fprintf(stdout,"\nTIME: %.3lf sec\n", (double)clock()/(CLOCKS_PER_SEC));

#define sz(x) (int)x.size()
#define ms(a,b) memset(a, b, sizeof(a))
#define FOR(i, a, n) for (ll i = (ll)a; i < (ll)n; ++i)
#define R(i, n) FOR(i, 0, n)
#define int long long
#define each(a,x)for(auto&a:(x))
#define bit(n,k) (((n)>>(k))&1)
#define pb push_back
#define sz(x) (int)x.size()
#define rep(i,n) for(int i=0;i<n;i++)
#define all(a) (a.begin()),(a.end())
#define D(x) cout << #x << " = " << x << endl
int gcd(int n,int m){return m==0?n:gcd(m,n%m);}
#define input freopen("/Users/shahadat/Desktop/input.txt", "r", stdin)
#define output freopen("/Users/shahadat/Desktop/output.txt", "w", stdout)
#define valid(x,y) x>=0 && y>=0 && x<n && y<n
#define x  first
#define y  second
#define UNIQUE(X) (X).erase(unique(all(X)),(X).end())
#define SORT_UNIQUE(c) (sort(c.begin(),c.end()), c.resize(distance(c.begin(),unique(c.begin(),c.end()))))
#define pr(x) {cout << (x) << '\n';  }
#define GET_POS(c,x) (lower_bound(c.begin(),c.end(),x)-c.begin())
#define mask(i) ((1LL<<(i)))
#define endl '\n'

#define debug(x) cout<<#x<<" :: "<<x<<endl;
#define debug2(x,y) cout<<#x<<" :: "<<x<<"\t"<<#y<<" :: "<<y<<endl;
#define debug3(x,y,z) cout<<#x<<" :: "<<x<<"\t"<<#y<<" :: "<<y<<"\t"<<#z<<" :: "<<z<<endl;

using namespace std;
typedef long long ll;
typedef pair<int, int> ii;
typedef vector<ii> vii;
typedef vector<int> vi;

const int INF = 1e18;
int const N = 1e6 + 10;
int const M = 210;
double const PI = acos(-1);
int const MOD = 998244353;



// sieve design for divisor from 1 to N
vector<int> factor[N + 10];
int sum = 0;
void sieve(){
    for(int i = 1; i<=N; i++){
        for(int j = i; j<=N; j += i){
            factor[j].push_back(i);
            sum++;
        }
    }
}


// prime sieve
int store[N];
bool prime[N];
void sieve(){
    for(int i = 4; i<=N; i+= 2)
        prime[i] = true;
    for(int i = 3; i*i<=N; i+=2){
        if(!prime[i]){
            for(int j = i*i; j<=N; j+=i){
                if(!prime[j]) prime[j] = true;
            }
        }
    }
    for(int i = 1; i<=N; i++)
        if(prime[i]) store.push_back(i);
}

// prime sieve for factorization
int prime[N + 10];
void sieve(){
    for(int i = 1; i <= N; i++) prime[i] = i;
    for(int i = 2; i<=N; i+=2) prime[i] = 2;
    for(long long i = 3; i*i <= N; i += 2){
        if(prime[i] != i) continue;
        for(long long j = i*i; j <= N; j += i)
            prime[j] = i;
    }
}

while ( prime[temp] != 1) {
    store.push_back(prime[temp]);
    temp /= prime[temp];
}




// prime linear sieve
std::vector <int> prime;
bool is_composite[N + 10];

void sieve (int n) {
    std::fill (is_composite, is_composite + n, false);
    for (int i = 2; i < n; ++i) {
        if (!is_composite[i]) prime.push_back (i);
        for (int j = 0; j < prime.size () && i * prime[j] < n; ++j) {
            is_composite[i * prime[j]] = true;
            if (i % prime[j] == 0) break;
        }
    }
}

// factorial factorization calculation
int calc(int div, int n){
    int sum = 0;
    while (n != 0) {
        sum += n/div;
        n /= div;
    }
    return sum;
}

// bitSize of a number
int bitS(int n){
    int sum = 0;
    while (n) { n >>= 1; sum++;}
    return sum;
}


// divisor & sum of divisor up to 1 to n
int divs[N + 10];
long long sig[N + 10];

void sieve(){
    for(int i = 1; i<=N; i++){
        for(int j = i; j<=N; j += i){
            divs[j]++;
            sig[j] += i;
        }
    }
}

//

string hextobin(char a)
{
    int nw = 0;
    if(a >= 'A') nw = a - 'A' + 10;
    else nw = a - '0';
    string ret = "";
    for(int i = 0; i<4; i++)
        ret += (nw&1) + '0', nw >>= 1;
    reverse(all(ret));
    return ret;
}

// random number generator & learn please
int get_rand()
{
    //return rand();
    return ( ( rand() << 15 ) ^ rand() );
} // int idx = 1 + get_rand() % n;

auto rnd = std::bind(std::uniform_int_distribution<int>(1,n), mt19937(seed)); rnd();
//https://codeforces.com/contest/1114/problem/E


// Z algo
vector<int> zfunc(string s){
    int n = s.length();
    vector<int> z(n);
    for(int i = 1, l = 0, r = 0; i<n; i++){
        if(i <= r)
            z[i] = min(r - i + 1, z[i - l]);
        while(i + z[i] < n && s[z[i]] == s[i + z[i]])
            ++z[i];
        if(i + z[i] - 1 > r)
            l = i, r = i + z[i] - 1;
    }
    return z;
}


// Articulation Point & bridge
int parent[N], dlow[N], visit[N], dnow[N];
int dcount = 0;
vector<int> store[N];

void bridge(int pos){
    visit[pos] = 1;
    dnow[pos] = dlow[pos] =  ++dcount;
    for(auto c: store[pos]){
        parent[c] = pos;
        if(visit[c] == 0){
            dfs(c);
            if(dlow[c] > dnow[pos]){
                cout<<pos<<" "<<c<<endl;
            }

            /* Articulation point & you have to care about parent nood (must have 2 or more child)
            if(dlow[c] >= dnow[pos]){
                cout<<pos<<" "<<c<<endl;
            }*/
            dlow[pos] = min(dlow[pos], dlow[c]);
        }
        else if(c != parent[pos]){
            dlow[pos] = min(dlow[pos], dnow[c]);
        }
    }
}


// failure table KMP
int lps[N];
void ftable(string& pattern)
{
    int L=0,R,len=pattern.size();
    for(R=1;R<len;){
        if(pattern[R]==pattern[L]){
            lps[R]=L+1;
            L++;R++;
        }
        else{
            if(L!=0) L=lps[L-1];
            else {lps[R]=0; R++;}
        }
    }
}


// DSU disjoint Set Union
void make_set(int x){
    for(int i = 1; i<=x; i++)
        parent[i] = i;
}
int find(int r){
    if(parent[r] == r) return r;
    return parent[r] = find(parent[r]);
}
void Union(int a, int b){
    int u = find(a);
    int v = find(b);
    if(u != v)
        parent[u] = v;
}


// range divide sum
int divideSum(int n, int g){
    int m = (n / g);
    int t1 = g * (m * (m - 1))/2;
    int t2 = (1 + (n % g))*m;
    return t1 + t2;
}

// gcd
int gcd(int a, int b){
    if(b == 0) return a;
    return gcd(b, a % b);
}
int gcd(int n,int m){return m==0?n:gcd(m,n%m);}


// big multiplication
int bigMul(int a, int b){
    int ret = 0;
    while (b)
    {   if(b & 1) ret = (ret + a) % MOD;
        a = (a + a) % MOD;
        b >>= 1;
    }
    return ret;
}

// big Mod

int bigMod(int b, int p){
    int res = 1;
    while(p != 0){
        if(p & 1)
            res = (res * b) % MOD;
        b = (b*b)%MOD;
        p >>= 1;
    }
    return res;
}
int modInv(int n){ return bigMod(n, MOD - 2)}


// suffix array
vector<int> suffArray(string const & str){
    int len = (int)str.length();
    int alpha = 256;
    vector<int> p(len), c(len), cnt(max(len, alpha), 0);
    
    for(int i = 0; i<len; i++){
        cnt[str[i]]++;
    }
    for(int i = 1; i<alpha; i++){
        cnt[i] += cnt[i-1];
    }
    for(int i = 0; i<len; i++){
        p[--cnt[str[i]]] = i;
    }
    c[p[0]] = 0;
    int cls = 1;
    for(int i = 1; i<len; i++){
        if(str[p[i]] != str[p[i-1]])
            cls++;
        c[p[i]] = cls - 1;
    }
    
    vector<int> pn(len), cn(len);
    
    for(int h = 0; (1 << h) < len; h++){
        
        for(int i = 0; i<len; i++){
            pn[i] = p[i] - (1 << h);
            if(pn[i] < 0) pn[i] += len;
        }
        fill(cnt.begin(), cnt.end(), 0);
        for(int i = 0; i<len; i++){
            cnt[c[pn[i]]]++;
        }
        for(int i = 1; i<cls; i++){
            cnt[i] += cnt[i-1];
        }
        for(int i = len-1; i>=0; i--){
            p[--cnt[c[pn[i]]]] = pn[i];
        }
        
        cn[p[0]] = 0;
        cls = 1;
        for(int i = 1; i<len; i++){
            pair<int,int> cur ={c[p[i]], c[(p[i] + (1 << h)) % len] };
            pair<int,int> prev = {c[p[i-1]], c[(p[i-1] + (1 << h)) % len]};
            if(cur != prev)
                cls++;
            cn[p[i]] = cls - 1;
        }
        c.swap(cn);
    }
    p.erase(p.begin(), p.begin() + 1);
    return p;
}

//modular Inverse 1 to n
void pre(){
    inv[1] = 1;
    int mod = MOD;
    for(int i = 2; i < 5010; ++i)
        inv[i] = (mod - mod / i) * inv[mod % i] % mod;
}

// hash
long long compute_hash(string const& s) {
    const int p = 31;
    const int m = 1e9 + 9;
    long long hash_value = 0;
    long long p_pow = 1;
    for (char c : s) {
        hash_value = (hash_value + (c - 'a' + 1) * p_pow) % m;
        p_pow = (p_pow * p) % m;
    }
    return hash_value;
}


// Longest Increasing Subsequence
int lis(vector<int>& ar){
    vector<int> store;
    store.push_back(ar[0]);
    
    for(int i = 1; i < ar.size(); i++){
        if(store.back() <= ar[i])
            store.push_back(ar[i]);
        else{
            int pos = GET_POS(store, ar[i]);
            store[pos] = ar[i];
        }
    }
    return (int)store.size();
}


