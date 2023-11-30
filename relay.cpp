#include <iostream>

struct line { // y = b - ax
	line& operator=(const line &m) { b = m.b; a = m.a; return *this; }
	long long eval(const long long x) const { return b - a*x; }
	long long b, a; 
};

long double intersectX(const line &m, const line &n) {
	return (long double) (m.b - n.b)/(m.a - n.a);
}

class convex_hull_trick {
public:
	convex_hull_trick() : j(0) {}
	convex_hull_trick(const convex_hull_trick &C) : j(C.j) {
		for(int i = 0; i < j; ++i) L[i] = C.L[i];
	}
	
	// insertions must be sorted by a (we need decreasing values of -a) and no parallel lines
	void insert(const long long b, const long long a) {
		const line m = { .b = b, .a = a };
		if (j > 1) while (intersectX(m, L[j-1]) <= intersectX(L[j-1], L[j-2]) && --j > 1);
		L[j] = m;
		j++;

		return;
	}

	long long min(const long long x) const {
		if (j < 2) return L[0].eval(x); // one element so far
		// else binary search the intersection points
		int low = 0, up = j-1, mid;
		long double x1, x2;
		while (low <= up) {
			mid = (low + up)/2;
			x1 = x2 = x;
			if (mid > 0)   x1 = intersectX(L[mid], L[mid-1]);
			if (mid < j-1) x2 = intersectX(L[mid], L[mid+1]);

			if (x1 <= x && x <= x2) return L[mid].eval(x);
			if (x < x1) up  = mid - 1;
			else        low = mid + 1;
		}

		// std::cerr << "Failed to find min x = " << x << std::endl;
		// print(std::cerr);
		return -1; // should not reach here
	}

private:
	line L[100005]; // Convex Hull Data Structure
	int j;          // number of lines currently included
	void print(std::ostream &out) const { // for debug
		out << '{' << std::endl;
		for (int i = 0; i < j; ++i)
			out << "\ty = " << L[i].b << " - " << L[i].a << 'x' << std::endl;
		out << '}' << std::endl;
	}
};

struct node {
	node(const long long v, const int d, const node* const n): V(v), Dij(d), next(n) {}
	long long V; // name of neighbour
	int Dij;
	const node *next;
};

long long P[100005], S[100005], DP[100005]; // 10^5 + 5
const node *L[100005];
bool explored[100005];


int get_n_unexplored_children(const long long v) {
	int n_children = 0;
	const node *n  = L[v];
	while(n != nullptr) {
		if(!explored[n->V]) ++n_children;
		n = n->next;
	}
	return n_children;
}

void DFS(const long long v, convex_hull_trick* const C, const long long D) { // D = D1v, the total distance from 1 to v
	explored[v] = true;

	if (v != 1) { // v = 1 is the base case where DP[1] = D[1] = 0 (that's why we inserted (0, 0) before calling this)
		DP[v] = P[v] + S[v]*D + C->min(S[v]); // DP[v] = Pv + Sv Dv + min_{u ancestor of v}(DP[u] - Du Sv)
		C->insert(DP[v], D);
	}

	int i = 0, n_children = get_n_unexplored_children(v);
	const node *n = L[v];
	while(n != nullptr) {
		if(!explored[n->V]) {
			// don't copy C for the last child (so if there is only one child don't copy C at all -> drops complexity in list cases from O(n^2) to O(n logn))
			// (generally it needs to be copied because any child may modify it but it shouldn't affect other children)
			convex_hull_trick* const H = ++i != n_children ? new convex_hull_trick(*C) : C;
			DFS(n->V, H, D + n->Dij);
			if (C != H) delete H;
		}
		n = n->next;
	}
}

int main() {
	int N, Vi, Vj, Dij;
	std::cin >> N;
	for (int i = 1; i < N; ++i) {
		std::cin >> Vi >> Vj >> Dij;
		L[Vi] = new node(Vj, Dij, L[Vi]);
		L[Vj] = new node(Vi, Dij, L[Vj]);
	}

	for (int i = 2; i <= N; ++i)
		std::cin >> P[i] >> S[i];

	convex_hull_trick* const C = new convex_hull_trick;
	C->insert(0, 0);
	DFS(1, C, 0);

	for (int i = 2; i < N; ++i)
		std::cout << DP[i] << ' ';
	std::cout << DP[N] << std::endl;

	return 0;
}
