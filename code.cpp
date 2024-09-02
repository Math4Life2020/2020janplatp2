#include <bits/stdc++.h> //AC with runtime ~1300ms on the worst case
using namespace std;
#pragma GCC optimize("O3,unroll-loops,strict-overflow,Ofast") //pragma build is not necessary to AC and only shaves ~200ms on the worst cases (~1500ms -> ~1300ms)
#pragma GCC target("avx2,abm,bmi,bmi2,sse4") //to be completely honest, I have no clue what half of these do

using ll = long long; using pii = pair<ll,ll>;
const ll p = 1e9+7;
const ll Nm = 65536;
const ll E = 16; //log2(Nm)
const ll K = 20;

ll v2(ll x) {
	if (x==0) {
		return 1000;
	}
	return __builtin_ctz(x);
}

struct mtr { //matrix for each component
	bool emp; //is empty?
	int mv[K][K]; //matrix values
};

struct vtr { //vector (for vector-matrix multiplications)
	ll vv[K];
};

mtr mtr0(ll x) { //matrix initializer
	mtr M;
	if (x==-1) {
		M.emp=1;
		return M;
	}
	M.emp=0;
	for (ll j=0;j<K;j++) {
		for (ll i=0;i<j;i++) {
			M.mv[i][j]=0;
		}
		M.mv[j][j]=1;
		for (ll i=j+1;i<K;i++) {
			M.mv[i][j]=0;
		}
	}
	for (ll i=0;i<=x;i++) {
		M.mv[i][x]++;
	}
	return M;
}

vtr vtr0() { //vector initializer
	vtr V;
	for (ll i=1;i<K;i++) {
		V.vv[i]=0;
	}
	V.vv[0]=1;
	return V;
}

mtr prd(mtr m1, mtr m2) { //matrix-matrix product
	if (m2.emp) {
		return m1;
	}
	mtr M;
	ll prdarr[K][K];
	M.emp=0;
	for (ll i=0;i<K;i++) {
		for (ll j=0;j<K;j++) {
			prdarr[i][j]=0;
		}
	}
	for (ll k=0;k<K;k++) {
		for (ll j=0;j<=k;j++) {
			for (ll i=0;i<=j;i++) {
				prdarr[i][k] += (ll)m1.mv[i][j]*(ll)m2.mv[j][k];
			}
		}
	}
	for (ll i=0;i<K;i++) {
		for (ll j=0;j<K;j++) {
			M.mv[i][j]=(prdarr[i][j])%p;
		}
	}
	return M;
}

vtr prdV(vtr v1, mtr m1) { //vector-matrix product; this is the bottleneck of the code
	vtr V;
	for (ll i=0;i<K;i++) {
		V.vv[i]=0;
	}
	for (ll j=0;j<K;j++) {
		for (ll i=0;i<=j;i++) {
			V.vv[j] = (V.vv[j]+v1.vv[i]*(ll)m1.mv[i][j]);
			if ((i&7)==0) { //taking modulo is the slowest operation here; doing it on each i-value gives TLE
				V.vv[j]%=p; //this uses ~1/8 the number of modulo operations (not including it in the loop here (aka only at the end) gives WA)
			}
		}
	}
	for (ll j=0;j<K;j++) {
		V.vv[j]%=p;
	}
	return V;
}

ll sumv(vtr v1) {
	ll val = 0;
	for (ll i=0;i<K;i++) {
		val = (val+v1.vv[i])%p;
	}
	return val;
}

mtr st[2*Nm]; //segment tree

int main() {
	ios_base::sync_with_stdio(false); cin.tie(0);
	freopen("nondec.in","r",stdin);
	freopen("nondec.out","w",stdout);
	ll N,K1; cin >> N >> K1;
	ll A[N];
	vtr vinit=vtr0(); 
	for (ll i=0;i<N;i++) {
		cin >> A[i]; A[i]--;
		st[Nm+i]=mtr0(A[i]);
	}
	for (ll i=N;i<Nm;i++) {
		st[Nm+i]=mtr0(-1);
	}
	for (ll i=(Nm-1);i>=1;i--) {
		st[i]=prd(st[2*i],st[2*i+1]); //constructs segtree of matrices
	}
	ll Q; cin >> Q;
	for (ll q=0;q<Q;q++) {
		ll l,r; cin >> l >> r;
		l--; r--;
		vtr v1 = vinit; //starts with the vector [1,0,0,0,...,0] and then multiplies it by the matrices (in the segtree)
		deque<ll> e1,e2;
		while (l<=r) {
			if (v2(l)<=v2(r+1)) {
				ll V = v2(l);
				e1.push_back((l>>V)+(1LL<<(E-V)));
				l += (1LL<<V);
			} else {
				ll V = v2(r+1);
				e2.push_front((r>>V)+(1LL<<(E-V)));
				r -= (1LL<<V);
			}
		}
		while (!e1.empty()) {
			ll x = e1.front(); e1.pop_front();
			v1 = prdV(v1,st[x]);
		}
		while (!e2.empty()) {
			ll x = e2.front(); e2.pop_front();
			v1 = prdV(v1,st[x]);
		}
		cout << sumv(v1) <<"\n";
	}
}
