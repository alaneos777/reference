#include <bits/stdc++.h>
using namespace std;
using ld = long double;
ld eps = 1e-9, inf = numeric_limits<ld>::max();
// For use with integers, just set eps=0 and everything remains the same
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b

struct point{
	ld x, y;
	point(): x(0), y(0){}
	point(ld x, ld y): x(x), y(y){}

	point operator+(const point & p) const{return point(x + p.x, y + p.y);}
	point operator-(const point & p) const{return point(x - p.x, y - p.y);}
	point operator*(const ld & k) const{return point(x * k, y * k);}
	point operator/(const ld & k) const{return point(x / k, y / k);}

	point operator+=(const point & p){*this = *this + p; return *this;}
	point operator-=(const point & p){*this = *this - p; return *this;}
	point operator*=(const ld & p){*this = *this * p; return *this;}
	point operator/=(const ld & p){*this = *this / p; return *this;}

	point rotate(const ld & angle) const{
		return point(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle));
	}
	point perp() const{return point(-y, x);}

	ld dot(const point & p) const{return x * p.x + y * p.y;}
	ld cross(const point & p) const{return x * p.y - y * p.x;}
	ld norm() const{return x * x + y * y;}
	ld length() const{return sqrtl(x * x + y * y);}
	point unit() const{return (*this) / length();}

	bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
	bool operator!=(const point & p) const{return !(*this == p);}
	bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
	bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}
	bool half(const point & p) const{return le(p.cross(*this), 0) || (eq(p.cross(*this), 0) && le(p.dot(*this), 0));}
};

istream &operator>>(istream &is, point & p){return is >> p.x >> p.y;}
ostream &operator<<(ostream &os, const point & p){return os << "(" << p.x << ", " << p.y << ")";}

int sgn(ld x){
	if(ge(x, 0)) return 1;
	if(le(x, 0)) return -1;
	return 0;
}

void polarSort(vector<point> & P, const point & o, const point & v){
	//sort points in P around o, taking the direction of v as first angle
	sort(P.begin(), P.end(), [&](const point & a, const point & b){
		return point((a - o).half(v), 0) < point((b - o).half(v), (a - o).cross(b - o));
	});
}

bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}

int intersectLinesInfo(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1 and a2+tv2
	ld det = v1.cross(v2);
	if(eq(det, 0)){
		if(eq((a2 - a1).cross(v1), 0)){
			return -1; //infinity points
		}else{
			return 0; //no points
		}
	}else{
		return 1; //single point
	}
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	//lines a1+tv1, a2+tv2
	//assuming that they intersect
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

int intersectLineSegmentInfo(const point & a, const point & v, const point & c, const point & d){
	//line a+tv, segment cd
	point v2 = d - c;
	ld det = v.cross(v2);
	if(eq(det, 0)){
		if(eq((c - a).cross(v), 0)){
			return -1; //infinity points
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v.cross(c - a)) != sgn(v.cross(d - a)); //1: single point, 0: no point
	}
}

int intersectSegmentsInfo(const point & a, const point & b, const point & c, const point & d){
	//segment ab, segment cd
	point v1 = b - a, v2 = d - c;
	int t = sgn(v1.cross(c - a)), u = sgn(v1.cross(d - a));
	if(t == u){
		if(t == 0){
			if(pointInSegment(a, b, c) || pointInSegment(a, b, d) || pointInSegment(c, d, a) || pointInSegment(c, d, b)){
				return -1; //infinity points
			}else{
				return 0; //no point
			}
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v2.cross(a - c)) != sgn(v2.cross(b - c)); //1: single point, 0: no point
	}
}

ld distancePointLine(const point & a, const point & v, const point & p){
	//line: a + tv, point p
	return abs(v.cross(p - a)) / v.length();
}

ld perimeter(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += (P[i] - P[(i + 1) % n]).length();
	}
	return ans;
}

ld area(vector<point> & P){
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return abs(ans / 2);
}

vector<point> convexHull(vector<point> P){
	sort(P.begin(), P.end());
	vector<point> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && leq((L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && leq((U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

bool pointInPerimeter(const vector<point> & P, const point & p){
	int n = P.size();
	for(int i = 0; i < n; i++){
		if(pointInSegment(P[i], P[(i + 1) % n], p)){
			return true;
		}
	}
	return false;
}

bool crossesRay(const point & a, const point & b, const point & p){
	return (geq(b.y, p.y) - geq(a.y, p.y)) * sgn((a - p).cross(b - p)) > 0;
}

int pointInPolygon(const vector<point> & P, const point & p){
	if(pointInPerimeter(P, p)){
		return -1; //point in the perimeter
	}
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += crossesRay(P[i], P[(i + 1) % n], p);
	}
	return rays & 1; //0: point outside, 1: point inside
}

//point in convex polygon in O(log n)
//make sure that P is convex and in ccw
//before the queries, do the preprocess on P:
// rotate(P.begin(), min_element(P.begin(), P.end()), P.end());
// int right = max_element(P.begin(), P.end()) - P.begin();
//returns 0 if p is outside, 1 if p is inside, -1 if p is in the perimeter
int pointInConvexPolygon(const vector<point> & P, const point & p, int right){
	if(p < P[0] || P[right] < p) return 0;
	int orientation = sgn((P[right] - P[0]).cross(p - P[0]));
	if(orientation == 0){
		if(p == P[0] || p == P[right]) return -1;
		return (right == 1 || right + 1 == P.size()) ? -1 : 1;
	}else if(orientation < 0){
		auto r = lower_bound(P.begin() + 1, P.begin() + right, p);
		int det = sgn((p - r[-1]).cross(r[0] - r[-1])) - 1;
		if(det == -2) det = 1;
		return det;
	}else{
		auto l = upper_bound(P.rbegin(), P.rend() - right - 1, p);
		int det = sgn((p - l[0]).cross((l == P.rbegin() ? P[0] : l[-1]) - l[0])) - 1;
		if(det == -2) det = 1;
		return det;
	}
}





vector<point> cutPolygon(const vector<point> & P, const point & a, const point & v){
	//returns the part of the convex polygon P on the left side of line a+tv
	int n = P.size();
	vector<point> lhs;
	for(int i = 0; i < n; ++i){
		if(geq(v.cross(P[i] - a), 0)){
			lhs.push_back(P[i]);
		}
		if(intersectLineSegmentInfo(a, v, P[i], P[(i+1)%n]) == 1){
			point p = intersectLines(a, v, P[i], P[(i+1)%n] - P[i]);
			if(p != P[i] && p != P[(i+1)%n]){
				lhs.push_back(p);
			}
		}
	}
	return lhs;
}
















point centroid(vector<point> & P){
	point num;
	ld den = 0;
	int n = P.size();
	for(int i = 0; i < n; ++i){
		ld cross = P[i].cross(P[(i + 1) % n]);
		num += (P[i] + P[(i + 1) % n]) * cross;
		den += cross;
	}
	return num / (3 * den);
}

vector<pair<int, int>> antipodalPairs(vector<point> & P){
	vector<pair<int, int>> ans;
	int n = P.size(), k = 1;
	auto f = [&](int u, int v, int w){return abs((P[v%n]-P[u%n]).cross(P[w%n]-P[u%n]));};
	while(ge(f(n-1, 0, k+1), f(n-1, 0, k))) ++k;
	for(int i = 0, j = k; i <= k && j < n; ++i){
		ans.emplace_back(i, j);
		while(j < n-1 && ge(f(i, i+1, j+1), f(i, i+1, j)))
			ans.emplace_back(i, ++j);
	}
	return ans;
}

pair<ld, ld> diameterAndWidth(vector<point> & P){
	int n = P.size(), k = 0;
	auto dot = [&](int a, int b){return (P[(a+1)%n]-P[a]).dot(P[(b+1)%n]-P[b]);};
	auto cross = [&](int a, int b){return (P[(a+1)%n]-P[a]).cross(P[(b+1)%n]-P[b]);};
	ld diameter = 0;
	ld width = inf;
	while(ge(dot(0, k), 0)) k = (k+1) % n;
	for(int i = 0; i < n; ++i){
		while(ge(cross(i, k), 0)) k = (k+1) % n;
		//pair: (i, k)
		diameter = max(diameter, (P[k] - P[i]).length());
		width = min(width, distancePointLine(P[i], P[(i+1)%n] - P[i], P[k]));
	}
	return {diameter, width};
}

pair<ld, ld> smallestEnclosingRectangle(vector<point> & P){
	int n = P.size();
	auto dot = [&](int a, int b){return (P[(a+1)%n]-P[a]).dot(P[(b+1)%n]-P[b]);};
	auto cross = [&](int a, int b){return (P[(a+1)%n]-P[a]).cross(P[(b+1)%n]-P[b]);};
	ld perimeter = inf, area = inf;
	for(int i = 0, j = 0, k = 0, m = 0; i < n; ++i){
		while(ge(dot(i, j), 0)) j = (j+1) % n;
		if(!i) k = j;
		while(ge(cross(i, k), 0)) k = (k+1) % n;
		if(!i) m = k;
		while(le(dot(i, m), 0)) m = (m+1) % n;
		//pairs: (i, k) , (j, m)
		point v = P[(i+1)%n] - P[i];
		ld h = distancePointLine(P[i], v, P[k]);
		ld w = distancePointLine(P[j], v.perp(), P[m]);
		perimeter = min(perimeter, 2 * (h + w));
		area = min(area, h * w);
	}
	return {area, perimeter};
}

ld distancePointCircle(const point & c, ld r, const point & p){
	//point p, circle with center c and radius r
	return max((ld)0, (p - c).length() - r);
}

point projectionPointCircle(const point & c, ld r, const point & p){
	//point p (outside the circle), circle with center c and radius r
	return c + (p - c).unit() * r;
}

pair<point, point> pointsOfTangency(const point & c, ld r, const point & p){
	//point p (outside the circle), circle with center c and radius r
	point v = (p - c).unit() * r;
	ld d2 = (p - c).norm(), d = sqrt(d2);
	point v1 = v * (r / d), v2 = v.perp() * (sqrt(d2 - r*r) / d);
	return {c + v1 - v2, c + v1 + v2};
}

vector<point> intersectLineCircle(const point & a, const point & v, const point & c, ld r){
	//line a+tv, circle with center c and radius r
	ld h2 = r*r - v.cross(c - a) * v.cross(c - a) / v.norm();
	point p = a + v * v.dot(c - a) / v.norm();
	if(eq(h2, 0)) return {p}; //line tangent to circle
	else if(le(h2, 0)) return {}; //no intersection
	else{
		point u = v.unit() * sqrt(h2);
		return {p - u, p + u}; //two points of intersection (chord)
	}
}

vector<point> intersectSegmentCircle(const point & a, const point & b, const point & c, ld r){
	//segment ab, circle with center c and radius r
	vector<point> P = intersectLineCircle(a, b - a, c, r), ans;
	for(const point & p : P){
		if(pointInSegment(a, b, p)) ans.push_back(p);
	}
	return ans;
}

pair<point, ld> getCircle(const point & m, const point & n, const point & p){
	//find circle that passes through points p, q, r
	point c = intersectLines((n + m) / 2, (n - m).perp(), (p + n) / 2, (p - n).perp());
	ld r = (c - m).length();
	return {c, r};
}

vector<point> intersectionCircles(const point & c1, ld r1, const point & c2, ld r2){
	//circle 1 with center c1 and radius r1
	//circle 2 with center c2 and radius r2
	point d = c2 - c1;
	ld d2 = d.norm();
	if(eq(d2, 0)) return {}; //concentric circles
	ld pd = (d2 + r1*r1 - r2*r2) / 2;
	ld h2 = r1*r1 - pd*pd/d2;
	point p = c1 + d*pd/d2;
	if(eq(h2, 0)) return {p}; //circles touch at one point
	else if(le(h2, 0)) return {}; //circles don't intersect
	else{
		point u = d.perp() * sqrt(h2/d2);
		return {p - u, p + u};
	}
}

int circleInsideCircle(const point & c1, ld r1, const point & c2, ld r2){
	//test if circle 2 is inside circle 1
	//returns "-1" if 2 touches internally 1, "1" if 2 is inside 1, "0" if they overlap
	ld l = r1 - r2 - (c1 - c2).length();
	return (ge(l, 0) ? 1 : (eq(l, 0) ? -1 : 0));
}

int circleOutsideCircle(const point & c1, ld r1, const point & c2, ld r2){
	//test if circle 2 is outside circle 1
	//returns "-1" if they touch externally, "1" if 2 is outside 1, "0" if they overlap
	ld l = (c1 - c2).length() - (r1 + r2);
	return (ge(l, 0) ? 1 : (eq(l, 0) ? -1 : 0));
}

int pointInCircle(const point & c, ld r, const point & p){
	//test if point p is inside the circle with center c and radius r
	//returns "0" if it's outside, "-1" if it's in the perimeter, "1" if it's inside
	ld l = (p - c).length() - r;
	return (le(l, 0) ? 1 : (eq(l, 0) ? -1 : 0));
}

vector<vector<point>> tangents(const point & c1, ld r1, const point & c2, ld r2, bool inner){
	//returns a vector of segments or a single point
	if(inner) r2 = -r2;
	point d = c2 - c1;
	ld dr = r1 - r2, d2 = d.norm(), h2 = d2 - dr*dr;
	if(eq(d2, 0) || le(h2, 0)) return {};
	point v = d*dr/d2;
	if(eq(h2, 0)) return {{c1 + v*r1}};
	else{
		point u = d.perp()*sqrt(h2)/d2;
		return {{c1 + (v - u)*r1, c2 + (v - u)*r2}, {c1 + (v + u)*r1, c2 + (v + u)*r2}};
	}
}

ld signed_angle(const point & a, const point & b){
	return sgn(a.cross(b)) * acosl(a.dot(b) / (a.length() * b.length()));
}

ld intersectPolygonCircle(const vector<point> & P, const point & c, ld r){
	//Gets the area of the intersection of the polygon with the circle
	int n = P.size();
	ld ans = 0;
	for(int i = 0; i < n; ++i){
		point p = P[i], q = P[(i+1)%n];
		bool p_inside = (pointInCircle(c, r, p) != 0);
		bool q_inside = (pointInCircle(c, r, q) != 0);
		if(p_inside && q_inside){
			ans += (p - c).cross(q - c);
		}else if(p_inside && !q_inside){
			point s1 = intersectSegmentCircle(p, q, c, r)[0];
			point s2 = intersectSegmentCircle(c, q, c, r)[0];
			ans += (p - c).cross(s1 - c) + r*r * signed_angle(s1 - c, s2 - c);
		}else if(!p_inside && q_inside){
			point s1 = intersectSegmentCircle(c, p, c, r)[0];
			point s2 = intersectSegmentCircle(p, q, c, r)[0];
			ans += (s2 - c).cross(q - c) + r*r * signed_angle(s1 - c, s2 - c);
		}else{
			auto info = intersectSegmentCircle(p, q, c, r);
			if(info.size() <= 1){
				ans += r*r * signed_angle(p - c, q - c);
			}else{
				point s2 = info[0], s3 = info[1];
				point s1 = intersectSegmentCircle(c, p, c, r)[0];
				point s4 = intersectSegmentCircle(c, q, c, r)[0];
				ans += (s2 - c).cross(s3 - c) + r*r * (signed_angle(s1 - c, s2 - c) + signed_angle(s3 - c, s4 - c));
			}
		}
	}
	return abs(ans)/2;
}

pair<point, ld> mec2(vector<point> & S, const point & a, const point & b, int n){
	ld hi = inf, lo = -hi;
	for(int i = 0; i < n; ++i){
		ld si = (b - a).cross(S[i] - a);
		if(eq(si, 0)) continue;
		point m = getCircle(a, b, S[i]).first;
		ld cr = (b - a).cross(m - a);
		if(le(si, 0)) hi = min(hi, cr);
		else lo = max(lo, cr);
	}
	ld v = (ge(lo, 0) ? lo : le(hi, 0) ? hi : 0);
	point c = (a + b) / 2 + (b - a).perp() * v / (b - a).norm();
	return {c, (a - c).norm()};
}

pair<point, ld> mec(vector<point> & S, const point & a, int n){
	random_shuffle(S.begin(), S.begin() + n);
	point b = S[0], c = (a + b) / 2;
	ld r = (a - c).norm();
	for(int i = 1; i < n; ++i){
		if(ge((S[i] - c).norm(), r)){
			tie(c, r) = (n == S.size() ? mec(S, S[i], i) : mec2(S, a, S[i], i));
		}
	}
	return {c, r};
}

pair<point, ld> smallestEnclosingCircle(vector<point> S){
	assert(!S.empty());
	auto r = mec(S, S[0], S.size());
	return {r.first, sqrt(r.second)};
}

bool comp1(const point & a, const point & b){
	return le(a.y, b.y);
}
pair<point, point> closestPairOfPoints(vector<point> P){
	sort(P.begin(), P.end(), comp1);
	set<point> S;
	ld ans = inf;
	point p, q;
	int pos = 0;
	for(int i = 0; i < P.size(); ++i){
		while(pos < i && geq(P[i].y - P[pos].y, ans)){
			S.erase(P[pos++]);
		}
		auto lower = S.lower_bound({P[i].x - ans - eps, -inf});
		auto upper = S.upper_bound({P[i].x + ans + eps, -inf});
		for(auto it = lower; it != upper; ++it){
			ld d = (P[i] - *it).length();
			if(le(d, ans)){
				ans = d;
				p = P[i];
				q = *it;
			}
		}
		S.insert(P[i]);
	}
	return {p, q};
}

struct vantage_point_tree{
	struct node
	{
		point p;
		ld th;
		node *l, *r;
	}*root;

	vector<pair<ld, point>> aux;

	vantage_point_tree(vector<point> &ps){
		for(int i = 0; i < ps.size(); ++i)
			aux.push_back({ 0, ps[i] });
		root = build(0, ps.size());
	}

	node *build(int l, int r){
		if(l == r)
			return 0;
		swap(aux[l], aux[l + rand() % (r - l)]);
		point p = aux[l++].second;
		if(l == r)
			return new node({ p });
		for(int i = l; i < r; ++i)
			aux[i].first = (p - aux[i].second).dot(p - aux[i].second);
		int m = (l + r) / 2;
		nth_element(aux.begin() + l, aux.begin() + m, aux.begin() + r);
		return new node({ p, sqrt(aux[m].first), build(l, m), build(m, r) });
	}

	priority_queue<pair<ld, node*>> que;

	void k_nn(node *t, point p, int k){
		if(!t)
			return;
		ld d = (p - t->p).length();
		if(que.size() < k)
			que.push({ d, t });
		else if(ge(que.top().first, d)){
			que.pop();
			que.push({ d, t });
		}
		if(!t->l && !t->r)
			return;
		if(le(d, t->th)){
			k_nn(t->l, p, k);
			if(leq(t->th - d, que.top().first))
				k_nn(t->r, p, k);
		}else{
			k_nn(t->r, p, k);
			if(leq(d - t->th, que.top().first))
				k_nn(t->l, p, k);
		}
	}

	vector<point> k_nn(point p, int k){
		k_nn(root, p, k);
		vector<point> ans;
		for(; !que.empty(); que.pop())
			ans.push_back(que.top().second->p);
		reverse(ans.begin(), ans.end());
		return ans;
	}
};

vector<point> minkowskiSum(vector<point> A, vector<point> B){
	int na = (int)A.size(), nb = (int)B.size();
	if(A.empty() || B.empty()) return {};

	rotate(A.begin(), min_element(A.begin(), A.end()), A.end());
	rotate(B.begin(), min_element(B.begin(), B.end()), B.end());

	int pa = 0, pb = 0;
	vector<point> M;

	while(pa < na && pb < nb){
		M.push_back(A[pa] + B[pb]);
		ld x = (A[(pa + 1) % na] - A[pa]).cross(B[(pb + 1) % nb] - B[pb]);
		if(leq(x, 0)) pb++;
		if(geq(x, 0)) pa++;
	}

	while(pa < na) M.push_back(A[pa++] + B[0]);
	while(pb < nb) M.push_back(B[pb++] + A[0]);

	return M;
}

//Delaunay triangulation in O(n log n)
const point inf_pt(inf, inf);

struct QuadEdge{
	point origin;
	QuadEdge* rot = nullptr;
	QuadEdge* onext = nullptr;
	bool used = false;
	QuadEdge* rev() const{return rot->rot;}
	QuadEdge* lnext() const{return rot->rev()->onext->rot;}
	QuadEdge* oprev() const{return rot->onext->rot;}
	point dest() const{return rev()->origin;}
};

QuadEdge* make_edge(const point & from, const point & to){
	QuadEdge* e1 = new QuadEdge;
	QuadEdge* e2 = new QuadEdge;
	QuadEdge* e3 = new QuadEdge;
	QuadEdge* e4 = new QuadEdge;
	e1->origin = from;
	e2->origin = to;
	e3->origin = e4->origin = inf_pt;
	e1->rot = e3;
	e2->rot = e4;
	e3->rot = e2;
	e4->rot = e1;
	e1->onext = e1;
	e2->onext = e2;
	e3->onext = e4;
	e4->onext = e3;
	return e1;
}

void splice(QuadEdge* a, QuadEdge* b){
	swap(a->onext->rot->onext, b->onext->rot->onext);
	swap(a->onext, b->onext);
}

void delete_edge(QuadEdge* e){
	splice(e, e->oprev());
	splice(e->rev(), e->rev()->oprev());
	delete e->rot;
	delete e->rev()->rot;
	delete e;
	delete e->rev();
}

QuadEdge* connect(QuadEdge* a, QuadEdge* b){
	QuadEdge* e = make_edge(a->dest(), b->origin);
	splice(e, a->lnext());
	splice(e->rev(), b);
	return e;
}

bool left_of(const point & p, QuadEdge* e){
	return ge((e->origin - p).cross(e->dest() - p), 0);
}

bool right_of(const point & p, QuadEdge* e){
	return le((e->origin - p).cross(e->dest() - p), 0);
}

ld det3(ld a1, ld a2, ld a3, ld b1, ld b2, ld b3, ld c1, ld c2, ld c3) {
	return a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) + a3 * (b1 * c2 - c1 * b2);
}

bool in_circle(const point & a, const point & b, const point & c, const point & d) {
	ld det = -det3(b.x, b.y, b.norm(), c.x, c.y, c.norm(), d.x, d.y, d.norm());
	det += det3(a.x, a.y, a.norm(), c.x, c.y, c.norm(), d.x, d.y, d.norm());
	det -= det3(a.x, a.y, a.norm(), b.x, b.y, b.norm(), d.x, d.y, d.norm());
	det += det3(a.x, a.y, a.norm(), b.x, b.y, b.norm(), c.x, c.y, c.norm());
	return ge(det, 0);
}

pair<QuadEdge*, QuadEdge*> build_tr(int l, int r, vector<point> & P){
	if(r - l + 1 == 2){
		QuadEdge* res = make_edge(P[l], P[r]);
		return {res, res->rev()};
	}
	if(r - l + 1 == 3){
		QuadEdge *a = make_edge(P[l], P[l + 1]), *b = make_edge(P[l + 1], P[r]);
		splice(a->rev(), b);
		int sg = sgn((P[l + 1] - P[l]).cross(P[r] - P[l]));
		if(sg == 0)
			return {a, b->rev()};
		QuadEdge* c = connect(b, a);
		if(sg == 1)
			return {a, b->rev()};
		else
			return {c->rev(), c};
	}
	int mid = (l + r) / 2;
	QuadEdge *ldo, *ldi, *rdo, *rdi;
	tie(ldo, ldi) = build_tr(l, mid, P);
	tie(rdi, rdo) = build_tr(mid + 1, r, P);
	while(true){
		if(left_of(rdi->origin, ldi)){
			ldi = ldi->lnext();
			continue;
		}
		if(right_of(ldi->origin, rdi)){
			rdi = rdi->rev()->onext;
			continue;
		}
		break;
	}
	QuadEdge* basel = connect(rdi->rev(), ldi);
	auto valid = [&basel](QuadEdge* e){return right_of(e->dest(), basel);};
	if(ldi->origin == ldo->origin)
		ldo = basel->rev();
	if(rdi->origin == rdo->origin)
		rdo = basel;
	while(true){
		QuadEdge* lcand = basel->rev()->onext;
		if(valid(lcand)){
			while(in_circle(basel->dest(), basel->origin, lcand->dest(), lcand->onext->dest())){
				QuadEdge* t = lcand->onext;
				delete_edge(lcand);
				lcand = t;
			}
		}
		QuadEdge* rcand = basel->oprev();
		if(valid(rcand)){
			while(in_circle(basel->dest(), basel->origin, rcand->dest(), rcand->oprev()->dest())){
				QuadEdge* t = rcand->oprev();
				delete_edge(rcand);
				rcand = t;
			}
		}
		if(!valid(lcand) && !valid(rcand))
			break;
		if(!valid(lcand) || (valid(rcand) && in_circle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())))
			basel = connect(rcand, basel->rev());
		else
			basel = connect(basel->rev(), lcand->rev());
	}
	return {ldo, rdo};
}

vector<tuple<point, point, point>> delaunay(vector<point> & P){
	sort(P.begin(), P.end());
	auto res = build_tr(0, (int)P.size() - 1, P);
	QuadEdge* e = res.first;
	vector<QuadEdge*> edges = {e};
	while(le((e->dest() - e->onext->dest()).cross(e->origin - e->onext->dest()), 0))
		e = e->onext;
	auto add = [&P, &e, &edges](){
		QuadEdge* curr = e;
		do{
			curr->used = true;
			P.push_back(curr->origin);
			edges.push_back(curr->rev());
			curr = curr->lnext();
		}while(curr != e);
	};
	add();
	P.clear();
	int kek = 0;
	while(kek < (int)edges.size())
		if(!(e = edges[kek++])->used)
			add();
	vector<tuple<point, point, point>> ans;
	for(int i = 0; i < (int)P.size(); i += 3){
		ans.emplace_back(P[i], P[i + 1], P[i + 2]);
	}
	return ans;
}

int main(){
	/*vector<pair<point, point>> centers = {{point(-2, 5), point(-8, -7)}, {point(14, 4), point(18, 6)}, {point(9, 20), point(9, 28)},
										  {point(21, 20), point(21, 29)}, {point(8, -10), point(14, -10)}, {point(24, -6), point(34, -6)},
										  {point(34, 8), point(36, 9)}, {point(50, 20), point(56, 24.5)}};
	vector<pair<ld, ld>> radii = {{7, 4}, {3, 5}, {4, 4}, {4, 5}, {3, 3}, {4, 6}, {5, 1}, {10, 2.5}};
	int n = centers.size();
	for(int i = 0; i < n; ++i){
		cout << "\n" << centers[i].first << " " << radii[i].first << " " << centers[i].second << " " << radii[i].second << "\n";
		auto extLines = tangents(centers[i].first, radii[i].first, centers[i].second, radii[i].second, false);
		cout << "Exterior tangents:\n";
		for(auto par : extLines){
			for(auto p : par){
				cout << p << " ";
			}
			cout << "\n";
		}
		auto intLines = tangents(centers[i].first, radii[i].first, centers[i].second, radii[i].second, true);
		cout << "Interior tangents:\n";
		for(auto par : intLines){
			for(auto p : par){
				cout << p << " ";
			}
			cout << "\n";
		}
	}*/

	/*int n;
	cin >> n;
	vector<point> P(n);
	for(auto & p : P) cin >> p;
	auto triangulation = delaunay(P);
	for(auto triangle : triangulation){
		cout << get<0>(triangle) << " " << get<1>(triangle) << " " << get<2>(triangle) << "\n";
	}*/

	/*int n;
	cin >> n;
	vector<point> P(n);
	for(auto & p : P) cin >> p;
	auto ans = smallestEnclosingCircle(P);
	cout << ans.first << " " << ans.second << "\n";*/

	/*vector<point> P;
	srand(time(0));
	for(int i = 0; i < 1000; ++i){
		P.emplace_back(rand() % 1000000000, rand() % 1000000000);
	}
	point o(rand() % 1000000000, rand() % 1000000000), v(rand() % 1000000000, rand() % 1000000000);
	polarSort(P, o, v);
	auto ang = [&](point p){
		ld th = atan2(p.y, p.x);
		if(th < 0) th += acosl(-1)*2;
		ld t = atan2(v.y, v.x);
		if(t < 0) t += acosl(-1)*2;
		if(th < t) th += acosl(-1)*2;
		return th;
	};
	for(int i = 0; i < P.size()-1; ++i){
		assert(leq(ang(P[i] - o), ang(P[i+1] - o)));
	}*/
	return 0;
}