#include "../include/Pipe.h"

//Pipe::Pipe():RCv(0), NCv(0), T(nullptr), P(nullptr), H(nullptr), W(nullptr), dT(nullptr), dP(nullptr), dH(nullptr), dW(nullptr),
//Length(0.0), Area(0.0), KCoe(0.0), De(0.0), Angle(0.0), FI(0.0), Thick(0.0)
//{}
FluType FluFlag(string flag) {
	if (flag == "Solid") return FluType::Solid;
	if (flag == "Fluid") return FluType::Fluid;
	if (flag == "gas")	return FluType::Gas;
}

HeatTransType HeatTransFlag(string flag)
{
	if (flag == "Cond") return HeatTransType::Cond;
	if (flag == "Conv") return HeatTransType::Conv;
	if (flag == "Rad") return HeatTransType::Rad;
}

Pipe::Pipe(Tint N0) :Htype(nullptr), Ftype(nullptr), N(2), RCv(nullptr), NCv(nullptr),
T(nullptr), P(nullptr), H(nullptr), W(nullptr), P_interface(nullptr),
dH(nullptr), dW(nullptr),dT(nullptr),
Length(nullptr), Area(nullptr), KCoe(nullptr), De(nullptr), Angle(0.0), FI(nullptr), Thick(nullptr),mat(nullptr)
{
	Htype = new HeatTransType[N];
	Ftype = new FluType[N];
	RCv = new Tint[N];
	NCv = new Tint[N];
	mat = new Material[N];
	Length = new Tdouble[N];
	Area = new Tdouble[N];
	Thick = new Tdouble[N];
}




Pipe::~Pipe()
{
	delete[]T;
	delete[]P;
	delete[]H;
	//delete[]W;
	delete[]P_interface;
	delete[]dH;
	delete[]dW;
	delete[]dT;
	delete[]Length;
	delete[]Area;
	delete[]KCoe;
	delete[]De;
	delete[]FI;
	delete[]Thick;
	delete[]mat;
	//delete boundary.Tin;
	//delete boundary.Tout;
	//delete boundary.Pin;
	//delete boundary.Pout;
	T = NULL;
	P = NULL;
	H = NULL;
	W = NULL;
	P_interface = NULL;
	dH = NULL;
	dW = NULL;
	dT = NULL;
	Length = NULL;
	Area = NULL;
	KCoe = NULL;
	De = NULL;
	FI = NULL;
	Thick = NULL;
	mat = NULL;
	boundary.Tin = NULL;
	boundary.Tout = NULL;
	boundary.Pin = NULL;
	boundary.Pout = NULL;
}




void Pipe::input(ifstream& FILE)
{	
	//-------------边界划分-----------------
	Tint CountSolid = 0;	//固体区域数目
	string flag;			//区域标识关键字
	string st;
	getline(FILE, st);
	stringstream ss(st);
	int i = 0;
	while (ss >> flag)
	{
		Ftype[i] = FluFlag(flag);
		switch (Ftype[i]) {
		case(FluType::Fluid):
			Htype[i] = HeatTransType::Conv;
			break;
		case(FluType::Solid):
			CountSolid++, Htype[i] = HeatTransType::Cond;
			break;
		case(FluType::Gas):
			Htype[i] = HeatTransType::Conv;
			break;
		}
		i++;
	}
	ss.clear();
	getline(FILE, st);
	ss.str(st);
	i = 0;
	while (ss >> flag)
	{

		mat[i] = MaterialFlag(flag);
		i++;
	}
	ss.str(""), ss.clear();
	KCoe = new Tdouble[N - CountSolid];
	De = new Tdouble[N - CountSolid];
	FI = new Tdouble[N - CountSolid];
	
	//-------------几何参数-----------------
	//除固体外径向节点均为1
	//各层材料轴向节点数一致,添加check函数
	Tint count = 0;
	Tint RCv_count = 0; //径向节点总数
	Tint NCv_count = 0; //轴向节点总数
	Tint RCv_count_fluid = 0;
	Tint NCv_count_fluid = 0;
	for (int j = 0; j < N; j++) //文件流输入跳单个空格
	{
		getline(FILE, st);
		ss << st;
		ss >> RCv[j] >> NCv[j];
		//cout << RCv[j] << NCv[j];
		ss.str(""), ss.clear();
		getline(FILE, st);
		ss << st;
		ss >> Length[j] >> Area[j] >> Thick[j];
		ss.str(""), ss.clear();
		if (Ftype[j] == FluType::Fluid) {
			getline(FILE, st);
			ss << st;
			ss >> KCoe[count] >> De[count] >> FI[count];
			RCv_count_fluid += RCv[j];
			NCv_count_fluid += NCv[j];
			ss.str(""), ss.clear();
		}
		RCv_count += RCv[j];
		NCv_count += NCv[j];
	}
	NCv_count = NCv_count / N;
	getline(FILE, st);
	ss << st;
	ss >> Angle;
	ss.str(""), ss.clear();
	//-------------热工参数-----------------
	Tdouble *t, *p, *w;		//热力学基本参数
	t = new Tdouble[N];
	p = new Tdouble[N-CountSolid];
	w = new Tdouble[N-CountSolid];
	for (int j = 0; j < N; j++) {
		getline(FILE, st);
		ss << st;
		if (Ftype[j] == FluType::Fluid)
			ss >> t[j] >> w[j] >> p[j];
		else
			ss >> t[j];
		ss.str(""), ss.clear();
	}

	T = new Tdouble[RCv_count * NCv_count];							//流体固体温度
	P = new Tdouble[RCv_count_fluid * NCv_count];					//流体系统压力
	P_interface = new Tdouble[RCv_count_fluid * (NCv_count -1)];	//流体界面压力
	H = new Tdouble[RCv_count_fluid * NCv_count ];					//流体焓
	W = new Tdouble[N - CountSolid];								//流量
	dH = new Tdouble[RCv_count_fluid * NCv_count]{ 0. };			
	dW = new Tdouble[N - CountSolid]{ 0. };
	dT = new Tdouble[CountSolid]{ 0. };

	Tint Count_T = 0;
	Tint Count_P = 0;
	Tint Count_H = 0;
	Tint Count_Pi = 0;
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < RCv[j] * NCv_count; k++)
		{
			T[Count_T + k] = t[j];
			if (Ftype[j] == FluType::Fluid) {
				P[Count_P + k] = p[j];
				if (k < RCv[j] * (NCv_count - 1))
					P_interface[Count_Pi + k] = p[j];
				H[Count_H + k] = mat[j].H(T[Count_H + k], P[Count_P + k]);
			}
		}
		Count_T += RCv[j] * NCv_count;
		if (Ftype[j] == FluType::Fluid)
		{
			Count_P += RCv[j] * NCv_count;
			Count_H += RCv[j] * NCv_count;
			Count_Pi += RCv[j] * (NCv_count - 1);
		}
	}
	delete[]t;
	delete[]w;
	delete[]p;
	t = NULL;
	w = NULL;
	p = NULL;
}
 

void Pipe::output(ofstream& FILE)
{
	FILE << "Pipe:" << setw(5) << N << " regions" << endl;
	for (int i = 0; i < N; i++) {
		if (Ftype[i] == FluType::Solid)
			FILE << "Region " << i + 1 << ": Solid" << endl;
		else
			FILE << "Region " << i + 1 << ": Solid" << endl;
		FILE << "RCv：" << setw(5) << RCv[i] << endl;
		FILE << "NCv: " << setw(5) << NCv[i] << endl;
		FILE << "Length: " << setw(5) << Length[i] << setw(10) << "Area: " << setw(5) << Area[i] << setw(10) << "Width: " << setw(5) << Thick[i] << endl;
		FILE << endl;
	}
}

void Pipe::steady()
{
	Tdouble *dPFric, *dPGra, *dPLoc, *dP;
	Tdouble *Dens, *Vel, *Re, dPsum;
	Tdouble NCv0;
	NCv0 = NCv[0];
	Tdouble PIn, POut;
	Tdouble TIn, Tout;
	Tdouble Win;
	dPsum = 0;
	dPFric = new Tdouble[NCv0]{ 0. };
	dPGra = new Tdouble[NCv0]{ 0. };
	dPLoc = new Tdouble[NCv0]{ 0. };
	dP = new Tdouble[NCv0]{ 0. };
	Dens = new Tdouble[NCv0]{ 0. };
	Vel = new Tdouble[NCv0]{ 0. };
	Re = new Tdouble[NCv0]{ 0. };
	Win = *(boundary.W);
	PIn = *(boundary.Pin);
	
	for (int j = 0; j < N; j++)
	{
		if (Win < 0.001 && Win > 0)
			Ftype[j] = FluType::Solid;
		else 
			if (Win < 0) Win = abs(Win), cout << "Massflow less than 0!" << endl;
	}

	for (int j = 0; j < N; j++)					//每个区域计算压降
	{
		if (Ftype[j] == FluType::Fluid)
		{
			for (int i = 0; i < NCv0; i++) {
				Dens[i] = mat->Dens(H[i], P[i]);
				Vel[i] = Win / Dens[i] / Area[j];
				Re[i] = Dens[i] * Vel[i] * De[j];
				dPFric[i] = 0.5 * Fric(Re[i]) * 0.5 * Length[j] / NCv[j] * Dens[i] * Vel[i]*Vel[i];
				dPGra[i] = Dens[i] * Gravity * Length[j] / NCv[j] * sin(Angle * Pi / 180.0);
				dP[i] = dPFric[i] + dPGra[i];
				dPsum += dP[i];
				if (i == 0)
				{
					P_interface[i] = PIn - dP[i];
					P[i] = PIn - dP[i] / 2;
				}
				else
				{
					if (i < NCv0 - 1) P_interface[i] = P_interface[i - 1] - dP[i];
					P[i] = P_interface[i - 1] - dP[i] / 2;
				}
					
			}
		}
	}
	
	POut = PIn - dPsum;
	*(boundary.Pout) = POut;
	delete[]dPFric;
	delete[]dPGra;
	delete[]dPLoc;
	delete[]dP;
	delete[]Dens;
	delete[]Vel;
	delete[]Re;

}




void Pipe::transnt(){}
void Pipe::set_boundary(){}

Tdouble Pipe::Fric(Tdouble Re) {
	Tdouble Fric0, Fric1;
	if (Re <= 1000.0)
		return 64 / Re;
	else if (Re > 2300 && Re < 1.e5)
		return 0.3164 / pow(Re, 0.25);
	else
		Fric0 = 0.3164 / pow(Re, 0.25);
	Fric1 = pow((1.0 / (-2.0 * log(2.51 / (Re * sqrt(Fric0))))), 2.0);
	while (abs(Fric1 - Fric0) > 1.0e-8)
	{
		Fric0 = Fric1;
		Fric1 = pow((1.0 / (-2.0 * log(2.51 / (Re * sqrt(Fric0))))), 2.0);
	}
	return Fric0;
}
