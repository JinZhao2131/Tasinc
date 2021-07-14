#ifndef PIPE_H
#define PIPE_H 1

#include <fstream>
#include "Material.h"
#include <iomanip>
#include <sstream>
#include <vector>
#include <numeric>

using namespace std;

typedef long double Tdouble;
typedef int Tint;

//传热方式：导热、对流传热、辐射传热
enum class HeatTransType{ Cond, Conv, Rad};
enum class FluType{Solid, Fluid, Gas};

static const Tdouble Gravity = 9.8;
static const Tdouble Pi = 3.14;

FluType FluFlag(string flag);


HeatTransType HeatTransFlag(string flag);


struct Boundary {
	const Tdouble* Tin = nullptr;
	Tdouble* Tout = nullptr;
	const Tdouble* Pin = nullptr;
	Tdouble* Pout = nullptr;
	const Tdouble* W = nullptr;
	Tdouble* dW = nullptr;
};

class Pipe {
protected:
	HeatTransType* Htype;
	FluType* Ftype;
	Tint N;
	Tint* RCv;
	Tint* NCv;
	Tdouble* T, * P, * H, * W, *P_interface; 
	Tdouble* dH, * dW, *dT;
	Tdouble *Length, *Area, *KCoe, *De, *FI, *Thick;
	Tdouble Angle;
	Material* mat;

public:
	Boundary boundary;
	
public:
	//Pipe();
	Pipe(Tint N0 = 2);
	~Pipe();
	void input(ifstream& FILE);
	void set_boundary();
	void output(ofstream& FILE);
	void steady();
	void transnt();
	Tdouble Fric(Tdouble Re);
};

#endif