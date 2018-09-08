#include <bits/stdc++.h>
using namespace std;

double eps = 1e-8;
# define M_PI 3.14159265358979323846
# define M_E 2.71828182845904523536

struct point{
	double x, y;

	point(): x(0), y(0){}

	point(double x, double y): x(x), y(y){}

	void read(){cin >> x >> y;}

	point operator+(const point & p) const{return point(x + p.x, y + p.y);}
	
	point operator-(const point & p) const{return point(x - p.x, y - p.y);}
	
	point operator*(const double & k) const{return point(x * k, y * k);}

	point operator/(const double & k) const{return point(x / k, y / k);}

	point operator+=(const point & p){*this = *this + p; return *this;}

	point operator-=(const point & p){*this = *this - p; return *this;}

	point operator*=(const double & p){*this = *this * p; return *this;}

	point operator/=(const double & p){*this = *this / p; return *this;}

	point rotate(const double angle) const{
		return point(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle));
	}
	point rotate(const double angle, const point & p){
		return p + ((*this) - p).rotate(angle);
	}
	point perpendicular() const{
		return point(-y, x);
	}

	double dot(const point & p) const{
		return x * p.x + y * p.y;
	}
	double length() const{
		return hypot(x, y);
	}
	double cross(const point & p) const{
		return x * p.y - y * p.x;
	}

	point normalize() const{
		return (*this) / length();
	}

	point projection(const point & p) const{
		return (*this) * p.dot(*this) / dot(*this);
	}
	point normal(const point & p) const{
		return p - projection(p);
	}

	bool operator==(const point & p) const{
		return abs(x - p.x) < eps && abs(y - p.y) < eps;
	}
	bool operator!=(const point & p) const{
		return !(*this == p);
	}
	bool operator<(const point & p) const{
		if(abs(x - p.x) < eps) return y < p.y;
		return x < p.x;
	}
	bool operator>(const point & p) const{
		if(abs(x - p.x) < eps) return y > p.y;
		return x > p.x;
	}
};

istream &operator>>(istream &is, point & P){
	point p;
    is >> p.x >> p.y;
    P = p;
    return is;
}

ostream &operator<<(ostream &os, const point & p) { 
    return os << "(" << p.x << ", " << p.y << ")";
}

int sgn(double x){
	if(abs(x) < eps){
		return 0;
	}else if(x > 0){
		return 1;
	}else{
		return -1;
	}
}


bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return abs((p - a).cross(v)) < eps;
}

bool pointInSegment(point a, point b, const point & p){
	//segment ab, point p
	if(a > b) swap(a, b);
	return pointInLine(a, b - a, p) && !(p < a || p > b);
}

int intersectLinesInfo(const point & a1, const point & v1, const point & a2, const point & v2){
	//line a1+tv1
	//line a2+tv2
	double det = v1.cross(v2);
	if(abs(det) < eps){
		if(abs((a2 - a1).cross(v1)) < eps){
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
	double det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

int intersectLineSegmentInfo(const point & a, const point & v, const point & c, const point & d){
	//line a+tv, segment cd
	point v2 = d - c;
	double det = v.cross(v2);
	if(abs(det) < eps){
		if(abs((c - a).cross(v)) < eps){
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

double distancePointLine(const point & a, const point & v, const point & p){
	//line: a + tv, point p
	return abs(v.cross(p - a)) / v.length();
}

double perimeter(vector<point> & P){
	int n = P.size();
	double ans = 0;
	for(int i = 0; i < n; i++){
		ans += (P[i] - P[(i + 1) % n]).length();
	}
	return ans;
}

double area(vector<point> & P){
	int n = P.size();
	double ans = 0;
	for(int i = 0; i < n; i++){
		ans += P[i].cross(P[(i + 1) % n]);
	}
	return abs(ans / 2);
}

vector<point> convexHull(vector<point> P){
	sort(P.begin(), P.end());
	vector<point> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && (L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]) <= 0){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && (U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]) <= 0){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

bool pointInPerimeter(vector<point> & P, const point & p){
	int n = P.size();
	for(int i = 0; i < n; i++){
		if(pointInSegment(P[i], P[(i + 1) % n], p)){
			return true;
		}
	}
	return false;
}

int pointInPolygon(vector<point> & P, const point & p){
	if(pointInPerimeter(P, p)){
		return -1; //point in the perimeter
	}
	point bottomLeft = (*min_element(P.begin(), P.end())) - point(M_E, M_PI);
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += (intersectSegmentsInfo(p, bottomLeft, P[i], P[(i + 1) % n]) == 1 ? 1 : 0);
	}
	return rays & 1; //0: point outside, 1: point inside
}

bool comp1(const point & a, const point & b){
	return a.y < b.y;
}
pair<point, point> closestPairOfPoints(vector<point> P){
	sort(P.begin(), P.end(), comp1);
	set<point> S;
	double ans = 1e9;
	point p, q;
	int pos = 0;
	for(int i = 0; i < P.size(); ++i){
		while(pos < i && abs(P[i].y - P[pos].y) >= ans){
			S.erase(P[pos++]);
		}
		auto lower = S.lower_bound({P[i].x - ans - eps, -1e9});
		auto upper = S.upper_bound({P[i].x + ans + eps, -1e9});
		for(auto it = lower; it != upper; ++it){
			double d = (P[i] - *it).length();
			if(d < ans){
				ans = d;
				p = P[i];
				q = *it;
			}
		}
		S.insert(P[i]);
	}
	return {p, q};
}

point centroid(vector<point> & P){
	point num;
	double den = 0;
	int n = P.size();
	for(int i = 0; i < n; ++i){
		double cross = P[i].cross(P[(i + 1) % n]);
		num += (P[i] + P[(i + 1) % n]) * cross;
		den += cross;
	}
	return num / (3.0 * den);
}

struct vantage_point_tree{
	struct node
	{
		point p;
		double th;
		node *l, *r;
	}*root;

	vector<pair<double, point>> aux;

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
			aux[i].first = (p - aux[i].second).length() * (p - aux[i].second).length();
		int m = (l + r) / 2;
		nth_element(aux.begin() + l, aux.begin() + m, aux.begin() + r);
		return new node({ p, sqrt(aux[m].first), build(l, m), build(m, r) });
	}

	priority_queue<pair<double, node*>> que;

	void k_nn(node *t, point p, int k){
		if(!t)
			return;
		double d = (p - t->p).length();
		if(que.size() < k)
			que.push({ d, t });
		else if(que.top().first > d){
			que.pop();
			que.push({ d, t });
		}
		if(!t->l && !t->r)
			return;
		if(d < t->th){
			k_nn(t->l, p, k);
			if(t->th - d <= que.top().first)
				k_nn(t->r, p, k);
		}else{
			k_nn(t->r, p, k);
			if(d - t->th <= que.top().first)
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

vector<pair<int, int>> antipodalPairs(vector<point> & P){
	vector<pair<int, int>> ans;
	int n = P.size(), k = 1;
	auto f = [&](int u, int v, int w){return abs((P[v%n]-P[u%n]).cross(P[w%n]-P[u%n]));};
	while(f(n-1, 0, k+1) > f(n-1, 0, k)) ++k;
	for(int i = 0, j = k; i <= k && j < n; ++i){
		ans.emplace_back(i, j);
		while(j < n-1 && f(i, i+1, j+1) > f(i, i+1, j))
			ans.emplace_back(i, ++j);
	}
	return ans;
}

pair<double, double> diameterAndWidth(vector<point> & P){
	int n = P.size(), k = 0;
	auto dot = [&](int a, int b){return (P[(a+1)%n]-P[a]).dot(P[(b+1)%n]-P[b]);};
	auto cross = [&](int a, int b){return (P[(a+1)%n]-P[a]).cross(P[(b+1)%n]-P[b]);};
	double diameter = numeric_limits<double>::min();
	double width = numeric_limits<double>::max();
	while(dot(0, k) > eps) k = (k+1) % n;
	for(int i = 0; i < n; ++i){
		while(cross(i, k) > eps) k = (k+1) % n;
		//pair: (i, k)
		diameter = max(diameter, (P[k] - P[i]).length());
		width = min(width, distancePointLine(P[i], P[(i+1)%n] - P[i], P[k]));
	}
	return make_pair(diameter, width);
}

pair<double, double> smallestEnclosingRectangle(vector<point> & P){
	int n = P.size();
	auto dot = [&](int a, int b){return (P[(a+1)%n]-P[a]).dot(P[(b+1)%n]-P[b]);};
	auto cross = [&](int a, int b){return (P[(a+1)%n]-P[a]).cross(P[(b+1)%n]-P[b]);};
	double perimeter = numeric_limits<double>::max(), area = perimeter;
	for(int i = 0, j = 0, k = 0, m = 0; i < n; ++i){
		while(dot(i, j) > eps) j = (j+1) % n;
		if(!i) k = j;
		while(cross(i, k) > eps) k = (k+1) % n;
		if(!i) m = k;
		while(dot(i, m) < -eps) m = (m+1) % n;
		//pairs: (i, k) , (j, m)
		point v = P[(i+1)%n] - P[i];
		double h = distancePointLine(P[i], v, P[k]);
		double w = distancePointLine(P[j], v.perpendicular(), P[m]);
		perimeter = min(perimeter, 2 * (h + w));
		area = min(area, h * w);
	}
	return make_pair(area, perimeter);
}

double distancePointCircle(const point & p, const point & c, double r){
	//point p, center c, radius r
	return max(0.0, (p - c).length() - r);
}

point projectionPointCircle(const point & p, const point & c, double r){
	//point p (outside the circle), center c, radius r
	return c + (p - c) / (p - c).length() * r;
}

pair<point, point> pointsOfTangency(const point & p, const point & c, double r){
	//point p (outside the circle), center c, radius r
	point v = (p - c).normalize() * r;
	double theta = acos(r / (p - c).length());
	return {c + v.rotate(-theta), c + v.rotate(theta)};
}

vector<point> intersectLineCircle(const point & a, const point & v, const point & c, double r){
	//line a+tv, center c, radius r
	double A = v.dot(v);
	double B = (a - c).dot(v);
	double C = (a - c).dot(a - c) - r * r;
	double D = B*B - A*C;
	if(abs(D) < eps) return {a + v * (-B/A)}; //line tangent to circle
	else if(D < 0) return {}; //no intersection
	else{ //two points of intersection (chord)
		D = sqrt(D);
		double t1 = (-B + D) / A;
		double t2 = (-B - D) / A;
		return {a + v * t1, a + v * t2};
	}
}

pair<point, double> getCircle(const point & m, const point & n, const point & p){
	//find circle that passes through points p, q, r
	point c = intersectLines((n + m) / 2, (n - m).perpendicular(), (p + n) / 2, (p - n).perpendicular());
	double r = (c - m).length();
	return {c, r};
}

vector<point> intersectionCircles(const point & c1, double r1, const point & c2, double r2){
	//circle 1 with center c1 and radius r1
	//circle 2 with center c2 and radius r2
	double A = 2*r1*(c2.y - c1.y);
	double B = 2*r1*(c2.x - c1.x);
	double C = (c1 - c2).dot(c1 - c2) + r1*r1 - r2*r2;
	double D = A*A + B*B - C*C;
	if(abs(D) < eps) return {c1 + point(B, A) * r1 / C};
	else if(D < 0) return {};
	else{
		D = sqrt(D);
		double cos1 = (B*C + A*D) / (A*A + B*B);
		double sin1 = (A*C - B*D) / (A*A + B*B);
		double cos2 = (B*C - A*D) / (A*A + B*B);
		double sin2 = (A*C + B*D) / (A*A + B*B);
		return {c1 + point(cos1, sin1) * r1, c1 + point(cos2, sin2) * r1};
	}
}

int circleInsideCircle(const point & c1, double r1, const point & c2, double r2){
	//test if circle 2 is inside circle 1
	//returns "-1" if 2 touches internally 1, "1" if 2 is inside 1, "0" if they overlap
	double l = r1 - r2 - (c1 - c2).length();
	return (l > eps ? 1 : (abs(l) < eps ? -1 : 0));
}

int circleOutsideCircle(const point & c1, double r1, const point & c2, double r2){
	//test if circle 2 is outside circle 1
	//returns "-1" if they touch externally, "1" if 2 is outside 1, "0" if they overlap
	double l = (c1 - c2).length() - (r1 + r2);
	return (l > eps ? 1 : (abs(l) < eps ? -1 : 0));
}

vector<vector<point>> commonExteriorTangents(const point & c1, double r1, const point & c2, double r2){
	//returns a vector of segments or a single point
	if(r1 < r2) return commonExteriorTangents(c2, r2, c1, r1);
	if(c1 == c2 && abs(r1-r2) < 0) return {};
	int in = circleInsideCircle(c1, r1, c2, r2);
	if(in == 1) return {};
	else if(in == -1) return {{c1 + (c2 - c1).normalize() * r1}};
	else{
		pair<point, point> t;
		if(abs(r1-r2) < eps)
			t = {c1 - (c2 - c1).perpendicular(), c1 + (c2 - c1).perpendicular()};
		else
			t = pointsOfTangency(c2, c1, r1 - r2);
		t.first = (t.first - c1).normalize();
		t.second = (t.second - c1).normalize();
		return {{c1 + t.first * r1, c2 + t.first * r2}, {c1 + t.second * r1, c2 + t.second * r2}};
	}
}

vector<vector<point>> commonInteriorTangents(const point & c1, double r1, const point & c2, double r2){
	if(c1 == c2 && abs(r1-r2) < 0) return {};
	int out = circleOutsideCircle(c1, r1, c2, r2);
	if(out == 0) return {};
	else if(out == -1) return {{c1 + (c2 - c1).normalize() * r1}};
	else{
		auto t = pointsOfTangency(c2, c1, r1 + r2);
		t.first = (t.first - c1).normalize();
		t.second = (t.second - c1).normalize();
		return {{c1 + t.first * r1, c2 - t.first * r2}, {c1 + t.second * r1, c2 - t.second * r2}};
	}
}

vector<point> minkowskiSum(vector<point> A, vector<point> B){
	int na = (int)A.size(), nb = (int)B.size();
	if(A.empty() || B.empty()) return {};

	rotate(A.begin(), min_element(A.begin(), A.end()), A.end());
	rotate(B.begin(), min_element(B.begin(), B.end()), B.end());

	int pa = 0, pb = 0;
	vector<point> M;

	while(pa < na && pb < nb){
		M.push_back(A[pa] + B[pb]);
		double x = (A[(pa + 1) % na] - A[pa]).cross(B[(pb + 1) % nb] - B[pb]);
		if(x <= eps) pb++;
		if(-eps <= x) pa++;
	}

	while(pa < na) M.push_back(A[pa++] + B[0]);
	while(pb < nb) M.push_back(B[pb++] + A[0]);

	return M;
}

int main(){
	vector<pair<point, point>> centers = {{point(-2, 5), point(-8, -7)}, {point(14, 4), point(18, 6)}, {point(9, 20), point(9, 28)},
										  {point(21, 20), point(21, 29)}, {point(8, -10), point(14, -10)}, {point(24, -6), point(34, -6)},
										  {point(34, 8), point(36, 9)}, {point(50, 20), point(56, 24.5)}};
	vector<pair<double, double>> radii = {{7, 4}, {3, 5}, {4, 4}, {4, 5}, {3, 3}, {4, 6}, {5, 1}, {10, 2.5}};
	int n = centers.size();
	for(int i = 0; i < n; ++i){
		cout << "\n" << centers[i].first << " " << radii[i].first << " " << centers[i].second << " " << radii[i].second << "\n";
		auto extLines = commonExteriorTangents(centers[i].first, radii[i].first, centers[i].second, radii[i].second);
		cout << "Exterior tangents:\n";
		for(auto par : extLines){
			for(auto p : par){
				cout << p << " ";
			}
			cout << "\n";
		}
		auto intLines = commonInteriorTangents(centers[i].first, radii[i].first, centers[i].second, radii[i].second);
		cout << "Interior tangents:\n";
		for(auto par : intLines){
			for(auto p : par){
				cout << p << " ";
			}
			cout << "\n";
		}
	}
	return 0;
}