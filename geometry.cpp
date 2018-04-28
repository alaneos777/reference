#include <bits/stdc++.h>
using namespace std;

double eps = 1e-8;
# define M_PI 3.14159265358979323846
# define M_E 2.71828182845904523536

struct point{
	double x, y;

	point(){
		x = y = 0;
	}
	point(double x, double y){
		this->x = x, this->y = y;
	}

	point operator+(const point & p) const{
		return point(x + p.x, y + p.y);
	}
	point operator-(const point & p) const{
		return point(x - p.x, y - p.y);
	}
	point operator*(const double & k) const{
		return point(x * k, y * k);
	}
	point operator/(const double & k) const{
		return point(x / k, y / k);
	}

	point rotate(const double angle) const{
		return point(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle));
	}
	point rotate(const double angle, const point & p){
		return p + ((*this) - p).rotate(angle);
	}

	double dot(const point & p) const{
		return x * p.x + y * p.y;
	}
	double length() const{
		return sqrt(dot(*this));
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
		if(abs(x - p.x) < eps){
			return y < p.y;
		}else{
			return x < p.x;
		}
	}
	bool operator>(const point & p) const{
		if(abs(x - p.x) < eps){
			return y > p.y;
		}else{
			return x > p.x;
		}
	}
};

istream &operator>>(istream &is, point & P){
	point p;
    is >> p.x >> p.y;
    P = p;
    return is;
}

ostream &operator<<(ostream &os, const point & p) { 
    return os << fixed << setprecision(8) << p.x << " " << p.y;
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

bool pointInLine(point & a, point & b, point & p){
	//line ab, point p
	return abs((p - a).cross(b - a)) < eps;
}

bool pointInSegment(point a, point b, point & p){
	//segment ab, point p
	if(a > b) swap(a, b);
	return pointInLine(a, b, p) && !(p < a || p > b);
}

int intersectLinesInfo(point & a, point & b, point & c, point & d){
	//line ab, line cd
	point v1 = b - a, v2 = d - c;
	double det = v1.cross(v2);
	if(abs(det) < eps){
		if(abs((c - a).cross(v1)) < eps){
			return -1; //infinity points
		}else{
			return 0; //no points
		}
	}else{
		return 1; //single point
	}
}

point intersectLines(point & a, point & b, point & c, point & d){
	//assuming that they intersect
	point v1 = b - a, v2 = d - c;
	double det = v1.cross(v2);
	return a + v1 * ((c - a).cross(v2) / det);
}

int intersectLineSegmentInfo(point & a, point & b, point & c, point & d){
	//line ab, segment cd
	point v1 = b - a, v2 = d - c;
	double det = v1.cross(v2);
	if(abs(det) < eps){
		if(abs((c - a).cross(v1)) < eps){
			return -1; //infinity points
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v1.cross(c - a)) != sgn(v1.cross(d - a)); //1: single point, 0: no point
	}
}

int intersectSegmentsInfo(point & a, point & b, point & c, point & d){
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

double distancePointLine(point & a, point & v, point & p){
	//line: a + tv, point p
	return abs(v.cross(p - a)) / v.length();
}

double perimeter(vector<point> & points){
	int n = points.size();
	double ans = 0;
	for(int i = 0; i < n; i++){
		ans += (points[i] - points[(i + 1) % n]).length();
	}
	return ans;
}

double area(vector<point> & points){
	int n = points.size();
	double ans = 0;
	for(int i = 0; i < n; i++){
		ans += points[i].cross(points[(i + 1) % n]);
	}
	return abs(ans / 2);
}

vector<point> convexHull(vector<point> points){
	sort(points.begin(), points.end());
	vector<point> L, U;
	for(int i = 0; i < points.size(); i++){
		while(L.size() >= 2 && (L[L.size() - 2] - points[i]).cross(L[L.size() - 1] - points[i]) <= 0){
			L.pop_back();
		}
		L.push_back(points[i]);
	}
	for(int i = points.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && (U[U.size() - 2] - points[i]).cross(U[U.size() - 1] - points[i]) <= 0){
			U.pop_back();
		}
		U.push_back(points[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

bool pointInPerimeter(vector<point> & points, point & p){
	int n = points.size();
	for(int i = 0; i < n; i++){
		if(pointInSegment(points[i], points[(i + 1) % n], p)){
			return true;
		}
	}
	return false;
}

int pointInPolygon(vector<point> & points, point & p){
	if(pointInPerimeter(points, p)){
		return -1; //point in the perimeter
	}
	point bottomLeft = (*min_element(points.begin(), points.end())) - point(M_E, M_PI);
	int n = points.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += (intersectSegmentsInfo(p, bottomLeft, points[i], points[(i + 1) % n]) == 1 ? 1 : 0);
	}
	return rays & 1; //0: point outside, 1: point inside
}