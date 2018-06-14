#pragma once
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#include<utility>

#define A_ARG 3//A>1
#define B_ARG 1.2//0.5<B<2
#define P_ARG 0.8//0<p<1
using namespace std;
class Pixel;
class Line;
extern vector<vector<Pixel>> picture_des;
extern vector<vector<Pixel>> picture_orign;//两个图像矩阵
extern vector<pair<Line, Line>> Lines;
extern int pic_1_W, pic_1_H, pic_2_W, pic_2_H, P_vectors_num;
class RGB_NEW
{
public:
	int R;
	int G;
	int B;
public:
	RGB_NEW(){}
	RGB_NEW(int s1, int s2, int s3) :R(s1), G(s2), B(s3){}
	RGB_NEW(RGB_NEW& item):R(item.R), G(item.G), B(item.B) {}
	explicit RGB_NEW(const RGB_NEW& item) :R(item.R), G(item.G), B(item.B) {}
	RGB_NEW& operator=(const RGB_NEW& item) { R = item.R; G = item.G; B = item.B; return *this; }
	friend const RGB_NEW operator*(const RGB_NEW& item1, const double item2) { RGB_NEW new_one; new_one.R = item1.R*item2, new_one.G = item1.G*item2, new_one.B = item1.B*item2; return new_one; }
	friend const RGB_NEW operator/(const RGB_NEW& item1, const double item2) { RGB_NEW new_one; new_one.R = item1.R/item2, new_one.G = item1.G/item2, new_one.B = item1.B/item2; return new_one; }
	friend const RGB_NEW operator+(const RGB_NEW& item1, const RGB_NEW& item2) { RGB_NEW new_one; new_one.R = item1.R + item2.R, new_one.G = item1.G + item2.G, new_one.B = item1.B + item2.B; return new_one; }
	void SetColor(RGB_NEW item) { R = item.R; G = item.G; B = item.B; }
};

class Point
{
public:
	int x;
	int y;
public:
	Point(){}
	Point(int x_pos, int y_pos) :x(x_pos), y(y_pos){}
	Point(Point& item):x(item.x),y(item.y){}
	Point(const Point& item) :x(item.x), y(item.y) {}
	Point& operator=(const Point& item) { x = item.x; y = item.y; return *this; }
	friend const Point operator-(const Point& item1, const Point& item2) { Point new_one; new_one.x = item1.x - item2.x; new_one.y = item1.y - item2.y; return new_one; }
	friend const Point operator+(const Point& item1, const Point& item2) { Point new_one; new_one.x = item1.x + item2.x; new_one.y = item1.y + item2.y; return new_one; }
	friend const Point operator*(const Point& item1, const double item2) { Point new_one; new_one.x = item1.x * item2; new_one.y = item1.y * item2; return new_one; }
	friend const Point operator/(const Point& item1, const double item2) { Point new_one; new_one.x = item1.x / item2; new_one.y = item1.y / item2; return new_one; }
	friend int Dot(Point item1, Point item2) { return item1.x*item2.x + item1.y*item2.y; }
	bool operator==(const Point& item) { return (x == item.x&&y == item.y); }
	bool operator!=(const Point& item) { return (x != item.x || y != item.y); }
	friend double Norm(Point item1, Point item2) { return sqrt(pow(item1.x - item2.x, 2) + pow(item1.y - item2.y, 2)); }
	void SetPoint(int x_pos, int y_pos) { x = x_pos; y = y_pos; }
};

class Line
{
public:
	Point P;
	Point Q;
public:
	Line(){}
	Line(int s1, int s2, int s3, int s4) :P(s1, s2), Q(s3, s4) {}
};

class Pixel
{
public:
	Point Pos;
	RGB_NEW Color;
public:
	Pixel() {}
	Pixel(int s1, int s2, int s3, int s4, int s5) :Pos(s1,s2),Color(s3,s4,s5){}
	Pixel(Pixel& item):Pos(item.Pos),Color(item.Color){}
	explicit Pixel(const Pixel& item) :Pos(item.Pos), Color(item.Color) {}
	void SetPos(int s1, int s2) { Pos.SetPoint(s1, s2); }
	void SetColor(RGB_NEW item) { Color.SetColor(item); }

};

void Tranfer_To_Line(vector<pair<Line, Line>>& lines, vector<string> s)
{
	Line l1(atoi(s[0].c_str()), atoi(s[1].c_str()),atoi(s[2].c_str()), atoi(s[3].c_str()));
	Line l2(atoi(s[4].c_str()), atoi(s[5].c_str()), atoi(s[6].c_str()), atoi(s[7].c_str()));
	lines.push_back(pair<Line, Line>(l1, l2));
}

void Tranfer_To_Pixel(vector<vector<Pixel>>& pic, vector<string> s,int n)
{
	vector<Pixel> pic_tmp;
	for (int i = 0; i < s.size(); ++i)
	{
		if (i != 0 && (i + 1) % 3 == 0)
			pic_tmp.push_back(Pixel((i + 1) / 3 - 1, n, atoi(s[i - 2].c_str()), atoi(s[i - 1].c_str()), atoi(s[i].c_str())));
	}
	pic.push_back(pic_tmp);
}

void Tranfer_To_Pixel_2(vector<vector<Pixel>>& pic, vector<string> s, int n)
{
	vector<Pixel> pic_tmp;
	for (int i = 0; i < s.size(); ++i)
	{
		if (i != 0 && (i + 1) % 3 == 0)
			pic_tmp.push_back(Pixel((i + 1) / 3 - 1, n - pic_1_H - 2, atoi(s[i - 2].c_str()), atoi(s[i - 1].c_str()), atoi(s[i].c_str())));
	}
	pic.push_back(pic_tmp);
}

void GetData(string& s)
{
	ifstream fin(s);
	string line;
	int flag = 0;//确保不执行第一行
	while (getline(fin, line))
	{
		if (flag == 0)
		{
			pic_1_W = atoi(line.c_str());
			++flag;
			continue;
		}
		if (flag == 1)
		{
			pic_1_H = atoi(line.c_str());
			++flag;
			continue;
		}
		if (flag >= 2 && flag <= pic_1_H + 1)
		{
			istringstream sin(line);
			vector<string> fields;
			string field;
			while (getline(sin, field, ' ')) //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
				fields.push_back(field); //将刚刚读取的字符串添加到向量fields中
			Tranfer_To_Pixel(picture_orign, fields,flag-2);
			++flag;
			continue;
		}
		if(flag==pic_1_H+2)
		{
			pic_2_W = atoi(line.c_str());
			++flag;
			continue;
		}
		if (flag == pic_1_H + 3)
		{
			pic_2_H = atoi(line.c_str());
			++flag;
			continue;
		}
		if (flag >= pic_1_H + 4 && flag <= pic_1_H + 4 + pic_2_H - 1)
		{
			istringstream sin(line);
			vector<string> fields;
			string field;
			while (getline(sin, field, ' ')) //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
				fields.push_back(field); //将刚刚读取的字符串添加到向量fields中
			Tranfer_To_Pixel_2(picture_des, fields, flag - 2);
			++flag;
			continue;
		}
		if (flag == pic_1_H + pic_2_H + 4)
		{
			P_vectors_num = atoi(line.c_str());
			++flag;
			continue;
		}
		if (flag > pic_1_H + pic_2_H + 4)
		{
			istringstream sin(line);
			vector<string> fields;
			string field;
			while (getline(sin, field, ' ')) //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
				fields.push_back(field); //将刚刚读取的字符串添加到向量fields中
			Tranfer_To_Line(Lines, fields);
			++flag;
			continue;
		}
	}
}

Point tranf_vertical(Point p1)
{
	int temp;
	temp = p1.x;
	p1.x = p1.y;
	p1.y = -temp;
	return p1;
}

vector<double> get_u(Pixel X_Des, vector<pair<Line, Line>>& line)
{
	vector<double> u_result;
	for (auto i = line.begin(); i != line.end(); ++i)
		u_result.push_back(Dot(Point(X_Des.Pos - (*i).second.P), (*i).second.Q - (*i).second.P) / pow(Norm((*i).second.Q, (*i).second.P), 2));
	return u_result;
}

vector<double> get_v(Pixel& X_Des, vector<pair<Line, Line>>& line)
{
	vector<double> v_result;
	for (auto i = line.begin(); i != line.end(); ++i)
	{
		Point Q_P = (*i).second.Q - (*i).second.P;
		v_result.push_back(Dot(X_Des.Pos - (*i).second.P, tranf_vertical(Q_P)) / Norm((*i).second.Q, (*i).second.P));
	}
	return v_result;
}

vector<Pixel> get_x_org(vector<double>& u, vector<double>& v,vector<pair<Line,Line>>& line)
{
	vector<Pixel> X_org;
	Point p;
	int index = 0;
	for (auto i = line.begin(); i != line.end(); ++i)
	{
		Point Q_P = (*i).first.Q - (*i).first.P;
		Point o = Q_P * u[index];
		p = (*i).first.P + Q_P * u[index] + tranf_vertical(Q_P)*v[index] / Norm((*i).first.Q, (*i).first.P);
		if (p.x >= pic_1_W)
			p.x = pic_1_W - 1;
		if (p.x < 0)
			p.x = 0;
		if (p.y >= pic_1_H)
			p.y = pic_1_H - 1;
		if (p.y < 0)
			p.y = 0;
		X_org.push_back(picture_orign[p.y][p.x]);
		++index;
	}
	return X_org;
}

vector<double> get_w(vector<double> u, vector<double> v, vector<pair<Line, Line>>& line,vector<Pixel> X_org)
{
	vector<double> pw;
	int index = 0;
	for (auto i = line.begin(); i != line.end(); ++i)
	{
		double dist;
		if (u[index] > 0 && u[index] < 1)
			dist = abs(v[index]);
		else if (u[index] < 0)
			dist = Norm(X_org[index].Pos, line[index].first.P);
		else
			dist = Norm(X_org[index].Pos, line[index].first.Q);
		pw.push_back(pow((pow(Norm((*i).first.Q, (*i).first.P), P_ARG) / (A_ARG + dist)), B_ARG));
	}
	return pw;
}

void Get_X_Orign_RGB(Pixel& X_Des,RGB_NEW& R_N)
{
	vector<double> u = get_u(X_Des, Lines);
	vector<double> v = get_v(X_Des, Lines);
	vector<Pixel> X_org = get_x_org(u, v, Lines);
	vector<double> w = get_w(u, v, Lines, X_org);
	RGB_NEW Sum_RGB(0,0,0);
	double Sum_w = 0;
	for (int i = 0; i != X_org.size(); ++i)
	{
		Sum_w += w[i];
		Sum_RGB =Sum_RGB+X_org[i].Color*w[i];
	}
	RGB_NEW R_N_1;
	if (X_org.size() == 1)
		R_N_1 = X_org[0].Color;
	else
		R_N_1 = Sum_RGB / Sum_w;
	R_N = (R_N_1 + X_Des.Color) / 2;
}

void Display_Des()
{
	string s = "output.txt";
	ofstream fout(s);
	if (fout.is_open())
	{
		fout << pic_1_W << endl;
		fout << pic_1_H << endl;
		for (int i = 0; i < pic_1_H; ++i)
		{
			for (int j = 0; j < pic_1_W; ++j)
			{
				fout << picture_des[i][j].Color.R << " " << picture_des[i][j].Color.G << " " << picture_des[i][j].Color.B << " ";
			}
			fout << endl;
		}
	}
	fout.close();
}