#include "../include/Material.h"

MaterialTab MaterialFlag(string flag)
{
if (flag == "NaK")	return MaterialTab::NaK;
if (flag == "AsFe") return MaterialTab::AsFe;
if (flag == "UO2")	return MaterialTab::UO2;
if (flag == "B4C")	return MaterialTab::B4C;
if (flag == "MOX")	return MaterialTab::MOX;
if (flag == "MoNb")	return MaterialTab::MoNb;
if (flag == "Mo")	return MaterialTab::Mo;
if (flag == "ZrH")	return MaterialTab::ZrH;
if (flag == "BeO")	return MaterialTab::BeO;
if (flag == "Be")	return MaterialTab::Be;
if (flag == "CO2")	return MaterialTab::CO2;
if (flag == "He")	return MaterialTab::He;
if (flag == "Xe")	return MaterialTab::Xe;
if (flag == "Cs")	return MaterialTab::Cs;
if (flag == "Wick")	return MaterialTab::Wick;
if (flag == "K")    return MaterialTab::K;
}

Material::Material( MaterialTab name0) {
	name = name0;
}

Material::Material(const Material& mat) {
	name = mat.name;
}


Material::~Material() {}

Tdouble Material::PfSat(Tdouble T)
{
	Tdouble P;
	switch (name) {
	case(MaterialTab::NaK):
		{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		P = exp(11.9463 - 12633.73/T-0.4672*log(T));
		return P;
		}
	case(MaterialTab::K):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		P = 4.0168e11 * pow(10, (-4625.3 / T)) * 1. / pow(T, 0.7);
		return P;
	}
	default:
		cout << "Error: Can not get the material!" << endl;

	}
}

Tdouble Material::PgSat(Tdouble T, Tdouble D)
{
	Tdouble P;
	switch (name) {
	case(MaterialTab::K):
	{
		Tdouble VV, CON1, CON2, CON3;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		VV = 1000 * 39.098 / D;
		CON1 = -0.3812E-5 * T * exp(7007 / T);
		CON2 = 6.193E-2 * exp(8168 / T);
		CON3 = -0.4615 * exp(10059 / T);
		P = (1 + CON1 / VV + CON2 / pow(VV, 2) + CON3 / pow(VV, 3)) / VV * T * 8.314;
		return P;
	}
	default:
		cout << "Error: Can not get the material!" << endl;

	}
}

Tdouble Material::TSat(Tdouble P)
{
	Tdouble T;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble A = 7.813;
		Tdouble B = 11209.0;
		Tdouble C = 5.249 * 100000.0;
		T = (B + sqrt(pow(B, 2) + 4.0 * C * (A - log(P)))) / (2.0 * (A - log(P)));
		Tdouble D = T + 1.0;
		while (abs(D-T)>0.001) {
			D = T;
			A = log(P);
			T = T - (11.9463 - log(P) - 12633.73 / T - 0.4672 * log(T)) / (12633.0 / pow(T, 2) - 0.4672 / T);
		}
		return T;
	}
	case(MaterialTab::K):
	{
		cout << " Waiting for setting! " << endl;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::VfSat(Tdouble T)
{
	Tdouble V;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		V = 393.37 * (1.0 - T / T0) + 4398.6 * pow((1.0 - T / T0), 0.29302);

		return V;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Hfg(Tdouble T)
{
	Tdouble H;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		H = 393.37 * (1.0 - T / T0) + 4398.6 * pow((1.0 - T / T0) ,0.29302);

		return H;
	}
	case(MaterialTab::K):
	{
		if (T < 336)
			H = 14.17 * 4.1858518;
		else
			H = 473.9 * 4.1858518;
		return H;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::VgSat(Tdouble T)
{
	Tdouble V;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble B;
		B = (12633.73 / (T * T) - 0.4672 / T) * exp(11.9463 - 12633.73 / T - 0.4672 * log(T));
		Tdouble T0 = 2503.7;
		V = VfSat(T) + Hfg(T) / (1000.0 * T * B);
		return V;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Vol(Tdouble T, Tdouble P)
{
	Tdouble V;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		Tdouble R = 8.3144;
		Tdouble A;
		A = 219.0 + 275.32 * (1.0 - T / T0) + 511.58 * sqrt(1.0 - T / T0);
		V = 1.0 / A;
		if (T < TSat(P))
			return V;
		else
		{
			Tdouble AN1,AN2,AN4, AP, AK2, AK4, AM, error;
			AN1 = 0.8;
			AP = P / 0.1013125;
			AK2 = exp(-9.95845 + 9117.50383 / T);
			AK4 = exp(-24.59115 + 20660.83322 / T);
			AN1 = 0.8;
			error = 0.1;
			while (error > 1.0e-12) {
				AN1 = AN1 - (AK4 * pow(AP, 3.) * pow(AN1, 4) + AK2 * AP * pow(AN1, 2.0) + AN1 - 1.0) / (4.0 * AK4 * pow(AP, 3.0) * pow(AN1, 3.0) + 2.0 * AK2 * AP * AN1 + 1.0);
				error = (AK4 * pow(AP, 3) * pow(AN1 ,4.0) + AK2 * AP * pow(AN1, 2.0) + AN1 - 1.0) / (4.0 * AK4 * pow(AP, 3.0) * pow(AN1, 3.0) + 2.0 * AK2 * AP * AN1 + 1.0);
			}
			AN2 = pow(AN1, 2.0) * AP * AK2;
			AN4 = pow(AN1, 4.0) * pow(AP, 3.0) * AK4;
			AM = 22.98977 * AN1 + 22.98977 * 2.0 * AN2 + 22.98977 * 4.0 * AN4;
			V = T * R / (1000.0 * P * AM);
			return V;
		}
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}

}

Tdouble Material::HfSat(Tdouble T) 
{
	Tdouble H;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		H = abs(1.6582 * T - 365.7 - 4.2395 * pow(T, 2.0) / 10000.0 + 1.487 * pow((T / 100.0), 3.0) / 10.0 + 2992.6 / T);
		return H;
	}
	case(MaterialTab::K):
	{
		Tdouble T0 = T - 273.15;
		if (T0 > 63.2 && T0 < 800)
			H = 56.179 + 0.84074 * T0 - 1.5844E-4 * pow(T0, 2) + 1.0499E-7 * pow(T0, 3);
		else
			H = 0.7104 * T0 + 1.0385E-3 * pow(T0, 2);
		return H;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::HgSat(Tdouble T)
{
	Tdouble H;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		H = Hfg(T)+HfSat(T);
		return H;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}

}

Tdouble Material::H(Tdouble T, Tdouble P)
{
	Tdouble H;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		H = -126.03 + 0.8844 * T;
		return H;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::DensdT(Tdouble T)
{
	Tdouble DensdT;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		DensdT = -275.32 / T0 - 511.58 * 0.5 /(sqrt(1.0 - T / T0) * T0);
		return DensdT;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Cond(Tdouble T)
{
	Tdouble Cond;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 13.9723 + 0.0331088 * T - 0.0000222398 * pow(T, 2);
		return Cond;
	}
	case(MaterialTab::AsFe):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 9.20 + 0.01750 * T - 2.0e-6 * T * T;
		return Cond;
	}
	case(MaterialTab::UO2):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 4575.0 / (1500 - 273.15 + 2360.0 - 2120.0 * pow(0.95, 3.0)) + 5.67e-10 * pow(0.95, 4.0) * pow((1500 - 273.15 - 860.0), 3.0);
		return Cond;
	}
	case(MaterialTab::B4C):
	{
		Tdouble TCB[13] = { 16.1 ,14.67 ,13.69 ,13.62 ,12.81 ,10.91 ,12.17 ,12.32 ,11.48 ,9.196 ,9.903,10.21,8.434 };
		Tdouble T00[13] = { 400. ,450. ,500. ,550. ,600. ,650. ,700. ,750. ,800. ,850. ,900. ,950. ,1000. };
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		for (int i = 0; i < 13; i++) {
			if ((T > T00[i]) && (T < T00[i + 1]))
				Cond = TCB[i] + (T - T00[i]) * (TCB[i + 1] - TCB[i]) / (T00[i + 1] - T00[i]);
		}
		if (T > 1000.0) Cond = TCB[12];
		return Cond;
	}
	case(MaterialTab::MoNb):
	{
		Tdouble A1, B1;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		A1 = 3.9e-3;
		B1 = 0.645;
		Cond = 120.0 - 1.4e-2 * T;
		return Cond;
	}
	case(MaterialTab::Mo):
	{
		Tdouble A1, B1;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = (110.0 - 0.015 * T) * 0.9;
		return Cond;
	}
	case(MaterialTab::ZrH):
	{
		Tdouble A1, B1;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 20.0;
		return Cond;
	}
	case(MaterialTab::BeO):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 4.0e7/pow(T, 2);
		return Cond;
	}
	case(MaterialTab::Be):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T <= 198.1)
			Cond = 304.362 - 0.472 * 198.1 + 3.38035e-4 * pow(198.1, 2.0) - 9.0729e-8 * pow(198.1, 3);
		else if (T >= 1556.0)
			Cond = 304.362 - 0.472 * 1556.0 + 3.38035e-4 * pow(1556.0, 2.0) - 9.0729e-8 * pow(1556.0, 3);
		else
			Cond = 304.362 - 0.472 * T + 3.38035e-4 * pow(T, 2.0) - 9.0729e-8 * pow(T, 3.0);
		return Cond;
	}
	case(MaterialTab::CO2):
	{
		Tdouble A1, B1;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		A1 = 3.9e-3;
		B1 = 0.645;
		Cond = 0.0064 + 14.4e-5 * T;
		return Cond;
	}
	case(MaterialTab::He):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 0.00031 * T + 0.059;
		return Cond;
	}
	case(MaterialTab::Xe):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 0.3;
		return Cond;
	}
	case(MaterialTab::Cs):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 0.15;
		return Cond;
	}
	case(MaterialTab::Wick):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = 0.15;
		return Cond;
	}
	case(MaterialTab::K):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T > 2.3 && T < 336.35)
			Cond = (0.301 - 1.44e-4 * T) * 418.68;
		else
			Cond = 43.8 - 2.22E-2 * T + 3.95E3 / T;
		return Cond;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Cond(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM)
{
	Tdouble Cond;
	switch (name) {
	case(MaterialTab::MOX):
	{
		Tdouble X = 2.0 - OM;
		Tdouble A = T / 1.0e3;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cond = (1.0 - P) / (1. + 2.0 * P) * (1.158 * (1.0 / (0.035 + 2.85 * X + (-0.715 * X + 0.286) * A) + 6400. * exp(-16.35 / A) / (pow(A, 2.5))));
		return Cond;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::CondSat(Tdouble T)
{
	Tdouble Cond;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Cond = (7.80 * pow(T,0.3) + 4.75 / T - 30.4 + 1.14e-2 * T) * 1.0e-3;
		return Cond;
	}
	case(MaterialTab::K):
	{
		Tdouble TT[7] = { 700, 800, 900,  1000, 1100, 1200, 1300 };
		Tdouble TH[7] = { 118, 148, 175,   200,	 225,  248,  269 };
		for (int i = 0; i < 6; i++)
		{
			if (T > TT[i] && T < TT[i + 1]) Cond = (T - TT[i]) / (TT[i + 1] - TT[i]) * (TH[i + 1] - TH[i]) + TH[i];
		}
		return Cond;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::CpfSat(Tdouble T)
{
	Tdouble Cp;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.7) T = 2503.7, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		Tdouble YXGM, CXGM, PL, AERFAX, AERFAP, BS, BT;
		YXGM = (12633.7 / pow(T, 2.0) - 0.4672 / T) * exp(11.9463 - 12633.7 / T - 0.4672 * log(T));
		PL = 219.0 + 275.32 * (1.0 - T / T0) + 511.58 * sqrt(1.0 - T / T0);
		if ((T >= 371.0) && (T < 2000.0))
			CXGM = 1.65820 - 8.47900 * T / 1.0e4 + 4.4541 * pow(T, 2) / 1.0e7 - 2992.6 / pow(T, 2) - 1000.0 * YXGM / PL;
		else
			CXGM = 0.86496 + 0.5 * (393.37 + 0.29302 * 4398.6 * pow((1.0 - T / 2503.70) ,(-0.70698))) / 2503.7 - 1000.0 * YXGM / PL;
		AERFAX = (275.32 / T0 + 511.58 * 0.5 / (sqrt(1.0 - T / T0) * T0)) / PL;
		BS = 1.717 * 0.0001 * (1 + (T - 371.0) / (T0 - 371.0) / 3.2682) / (1.0 - (T - 371.0) / (T0 - 371.0));
		BT = (BS * CXGM + 1000.0 * T / PL * AERFAX * (AERFAX + BS * YXGM)) / (CXGM - 1000.0 * T / PL * YXGM * (AERFAX + BS * YXGM));
		AERFAP = AERFAX + BT * YXGM;
		Cp = CXGM + 1000.0 * T * AERFAP * YXGM / PL;
		return Cp;
	}
	case(MaterialTab::K):
	{
		if (T > 336.35 && T < 1523.15)
			Cp = (0.8389 - 3.6741e-4 * (T - 273.15) + 4.592E-7 * pow((T - 273.15), 2)) * 1.E3;
		else
		{
			Cp = (0.1285 + 1.9116E-4 * T) * 4.18 * 1.E3;
			cout << "Properties need to be extended! " << endl;
		}
		return Cp;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Cp(Tdouble T)
{
	Tdouble Cp;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Tdouble T0 = 2503.7;
		Cp = (1097.73 - 0.556577 * T + 3.43167e-4 * pow(T, 2)) / 1000.0;
		return Cp;
	}
	case(MaterialTab::AsFe):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = (472.0 + 0.1360 * T - 2.82e6 / T / T) / 1000.0;
		return Cp;
	}
	case(MaterialTab::UO2):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T < 700)
			Cp = 297.0 + 2.5e-2 * (700 - 273.15) - 6.12e6 / (pow((700 - 273.15), 2.0));
		else
			Cp = 297.0 + 2.5e-2 * (T - 273.15) - 6.12e6 / (pow((T - 273.15), 2.0));
		Cp = 1.0 * Cp * 1.0e-3;
		return Cp;
	}
	case(MaterialTab::B4C):
	{
		Tdouble TCB[22] = { 0.959 ,1.206 ,1.341 ,1.468 ,1.589 ,1.704 ,1.812 ,1.913 ,2.008 ,2.096 ,2.178 ,2.253 ,2.321 ,2.383 ,2.438 ,2.487 ,2.529 ,2.564 ,2.593 ,2.615 ,2.631 ,2.640  };
		Tdouble T00[22] = { 300.0,500.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0 };
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		for (int i = 0; i < 13; i++) {
			if ((T > T00[i]) && (T < T00[i + 1]))
				Cp = TCB[i] + (T - T00[i]) * (TCB[i + 1] - TCB[i]) / (T00[i + 1] - T00[i]);
		}
		if (T > 1000.0) Cp = TCB[21];
		return Cp;
	}
	case(MaterialTab::MoNb):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 65.05e-3 + 0.0133e-3 * T;
		return Cp;
	}
	case(MaterialTab::Mo):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 101.8e-3 + 0.025e-3 * T;
		return Cp;
	}
	case(MaterialTab::ZrH):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = (187.0 + 0.745 * T) / 1000.0;
		return Cp;
	}
	case(MaterialTab::BeO):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = (1343.0 + 1.8469 * (T - 0.0)) / 1000.0;
		return Cp;
	}
	case(MaterialTab::Be):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T <= 223.1)
			Cp = 142.3;
		else if (T >= 1556.0)
			Cp = 357.7;
		else
			Cp = 256.23 + 6.58 * T - 5.262e-3 * pow(T, 2.0) + 1.54667e-6 * pow(T, 3.0);
		Cp = Cp / 1000.0;
		return Cp;
	}
	case(MaterialTab::CO2):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 100;
		return Cp;
	}
	case(MaterialTab::He):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 10000;
		return Cp;
	}
	case(MaterialTab::Xe):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 100;
		return Cp;
	}
	case(MaterialTab::Cs):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 100;
		return Cp;
	}
	case(MaterialTab::Wick):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Cp = 100;
		return Cp;
	}
	case(MaterialTab::K):
	{
		if (T > 100 && T < 336)
			Cp = 537.9 + 0.8002;
		else if (T > 336.0 && T < 1523.15)
			Cp = (0.8389 - 3.6741e-4 * (T - 273.15) + 4.592E-7 * pow((T - 273.15), 2)) * 1.E3;
		else
			cout << "Properties need to be fixed!";
		return Cp;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Cp(Tdouble T, Tdouble Y)
{
	Tdouble Cp;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble C0, C1, C2, C3, C4, C5, C6, T1, T2, CPT1, CPT2, CPMOX1, T0;
		Tdouble THEITA, EA, CPU, CPPU, C7;
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		T0 = T / 1.0e3;
		C1 = -78.4303;
		C2 = 193.238;
		C3 = 162.8647;
		C4 = -104.0014;
		C5 = 29.2056;
		C6 = -1.9507;
		C7 = 2.6441;
		CPU = C2 + 2.0 * C3 * T0 + 3.0 * C4 * pow(T0, 2.0) + 4.0 * C5 * pow(T0, 3.0)+ 5.0 * C6 * pow(T0, 4.0) - C7 * pow(T0, -2.0);
		C1 = -118.2062;
		C2 = 311.7866;
		C3 = 19.629;
		C4 = -0.752;
		C5 = 0.;
		C6 = 0.;
		C7 = 7.0131;
		CPPU = C2 + 2.0 * C3 * T0 + 3.0 * C4 * pow(T0, 2.0) + 4.0 * C5 * pow(T0, 3.0) + 5.0 * C6 * pow(T0, 4.0) - C7 * pow(T0, -2.0);
		Cp = (1 - Y) * CPU + Y * CPPU;
		Cp = Cp / 1.0e3;
		return Cp;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}


Tdouble Material::VCp(Tdouble T)
{
	Tdouble VCp;
	switch (name) {
	case(MaterialTab::UO2):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T < 2373.15)
			VCp = 0.147248498e4 + 0.4770160190e1 * T - 0.4261192860e-2 * T * T + 0.14571965e-5 * pow(T, 3) - 0.1090418860e-9 * pow(T, 4.0);
		else
			VCp = -0.882635e5 + 0.9059268180e2 * T - 0.2918222540e-1 * T * T + 0.31788104e-5 * pow(T, 3.0);
		return VCp;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}



Tdouble Material::T(Tdouble H, Tdouble P)
{
	Tdouble T;
	switch (name) {
	case(MaterialTab::NaK):
	{
		T = 142.54 + 1.13064 * H;
		return T;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Dens(Tdouble H, Tdouble P)
{
	Tdouble Dens;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble T0;
		T0 = T(H, P);
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 938.399 - 0.215438 * T0 - 1.71386e-5 * pow(T0, 2);
		return Dens;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Dens(Tdouble T0)
{
	Tdouble Dens;
	switch (name) {
	case(MaterialTab::AsFe):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 8084.0 - 0.42090 * T0 - 3.8940 * pow(T0,2.0) / 1.0e5;
		return Dens;
	}
	case(MaterialTab::UO2):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 10400.0;
		return Dens;
	}
	case(MaterialTab::B4C):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 2300.0;
		return Dens;
	}
	case(MaterialTab::MoNb):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 11506.0;
		return Dens;
	}
	case(MaterialTab::Mo):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 10200.0;
		return Dens;
	}
	case(MaterialTab::ZrH):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 5615.0;
		return Dens;
	}
	case(MaterialTab::BeO):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 2800.0;
		return Dens;
	}
	case(MaterialTab::Be):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 1830.0;
		return Dens;
	}
	case(MaterialTab::CO2):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 1.1;
		return Dens;
	}
	case(MaterialTab::He):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 4.0 / 22.4;
		return Dens;
	}
	case(MaterialTab::Xe):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 1.0;
		return Dens;
	}
	case(MaterialTab::Cs):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 1.0;
		return Dens;
	}
	case(MaterialTab::Wick):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 1.0;
		return Dens;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Dens(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM)
{
	Tdouble Dens;
	switch (name) {
	case(MaterialTab::AsFe):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 2503.07) T = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Dens = 10970.0 + Y * 490.0;
		Dens = Dens * (1 - P) * pow((1 + Expand(T, P, Y, OM)), -3);
		return Dens;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}


Tdouble Material::DensdH(Tdouble H0, Tdouble P)
{
	Tdouble DensdH;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble PP, X2, TS, A, B, X, HF, HG, X1, HFG, T0;
		T0 = T(H0, P);
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		TS = TSat(P);
		if (H0 < HfSat(TS))
			if (H0 <= (HfSat(TS) - 1.1))
			{
				A = T(H0 + 1.0, P);
				B = T(H0 - 1.0, P);
				DensdH = -(Vol(A, P) - Vol(B, P)) / (2.0 * pow(Vol(T0, P), 2));
			}
			else
			{
				A = T(H0, P);
				B = T(H0 - 0.001, P);
				DensdH = -(VfSat(A) - VfSat(B)) / (0.001 * pow(Vol(T0, P), 2));
			}
		else if (H0 >= HfSat(TS) && H0 <= HgSat(TS))
		{
			HF = HfSat(TS);
			HG = HgSat(TS);
			HFG = HG - HF;
			X = (H0 - HF) / HFG;
			if (H0 < (HfSat(TS) + 0.001))
			{
				X1 = (H0 + 0.001 - HF) / HFG;
				X2 = (H0 - HF) / HFG;
				DensdH = (1.0 / (VfSat(TS) * (1.0 - X1) + VgSat(TS) * X1) - 1.0 / (VfSat(TS) * (1.0 - X2) + VgSat(TS) * X2)) / 0.001;
			}
			else if (H0 >= (HfSat(TS) + 0.001) && H0 <= (HgSat(TS) - 0.001))
			{
				X1 = (H0 + 0.001 - HF) / HFG;
				X2 = (H0 - 0.001 - HF) / HFG;
				DensdH = (1.0 / (VfSat(TS) * (1.0 - X1) + VgSat(TS) * X1) - 1.0 / (VfSat(TS) * (1.0 - X2) + VgSat(TS) * X2)) / 0.002;
			}
			else
			{
				X1 = (H0 - HF) / HFG;
				X2 = (H0 - 0.001 - HF) / HFG;
				DensdH = (1.0 / (VfSat(TS) * (1.0 - X1) + VgSat(TS) * X1) - 1.0 / (VfSat(TS) * (1.0 - X2) + VgSat(TS) * X2)) / 0.001;
			}
		}
		else // if (H0 > HgSat(TS))
		{
			if (H0 <= (HgSat(TS) + 1.0))
			{
				A = T(H0 + 1.0, P);
				B = T(H0, P);
				DensdH = -(Vol(A, P) - Vol(B, P)) / (1.0 * pow(Vol(T0, P), 2));
			}
			else
			{
				A = T(H0 + 1.0, P);
				B = T(H0 - 1.0, P);
				DensdH = -(Vol(A, P) - Vol(B, P)) / (2.0 * pow(Vol(T0, P), 2));
			}
		}
		return DensdH;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::DensdP(Tdouble H, Tdouble P)
{
	Tdouble DensdP;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble T0;
		Tdouble T1, X2, A, B, TS, X, HF, HG, X1, xe, TS2, PM, PF, PG, DHGS, DHFS, DPFS, DPGS, ERROR, APP, K2, AN1, K4, A2, A4, F, DFDY1, DHAS, Y1, Y2, Y4, AK2, AK4, DY1DP, DENDP, DK2DT, DY1DT, DADT, DDHAST, EA2, EA1, EB2, EB1, DDHASDT, DADP, KTV, DENDT, HFG, CXGM, YXGM, PL, AERFAX, AERFAP, BS, BT, SDENPHDPSO1, TS1;
		T0 = T(H, P);
		TS = TSat(P);
		T1 = 2503.7;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (H < HfSat(TS))
		{
			PL = 1 / Vol(T0, P);
			YXGM = (12633.7 / pow(T0, 2.0) - 0.4672 / T0) * exp(11.9463 - 12633.7 / T0 - 0.4672 * log(T0));
			if (T0 >= 371.0 && T0 <= 2000.0)
				CXGM = 1.6582 - 8.4790 * T0 / 1.0e4 + 4.4541 * pow(T0, 2) / 1.0e7 - 2992.6 / pow(T0, 2.0) - 1000.0 * YXGM / PL;
			else //if (T0 > 2000 && T0 < 2503.7)
				CXGM = 0.86496 + 0.5 * (393.37 + 0.29302 * 4398.6 * pow((1.0 - T0 / 2503.7), -0.70698)) / 2503.7 - 1000.0 * YXGM / PL;
			AERFAX = (275.32 / T1 + 511.58 * 0.5 / (sqrt(1.0 - T0 / T1) * T1)) / PL;
			BS = 1.717 * 0.0001 * (1.0 + (T0 - 371.0) / (T1 - 371.0) / 3.2682) / (1.0 - (T0 - 371.0) / (T1 - 371.0));
			BT = (BS * CXGM + 1000.0 * T0 / PL * AERFAX * (AERFAX + BS * YXGM)) / (CXGM - 1000.0 * T0 / PL * YXGM * (AERFAX + BS * YXGM));
			AERFAP = AERFAX + BT * YXGM;
			DensdP = BT * PL + AERFAP * (1.0 - T0 * AERFAP) / CpfSat(TS);
		}
		else if (H >= HfSat(TS) && H <= HgSat(TS))
		{
			A = P + 0.0001;
			B = P + 0.0001;
			if (H < VfSat(TSat(A)))
				A = P;
			else if (H > VgSat(TSat(B)))
				B = P;
			TS1 = TSat(A);
			TS2 = TSat(B);
			HF = HfSat(TS1);
			HG = HgSat(TS1);
			HFG = HG - HF;
			X1 = (H - HF) / HFG;
			HF = HfSat(TS2);
			HG = HgSat(TS2);
			HFG = HG - HF;
			X2 = (H - HF) / HFG;
			DensdP = (1.0 / (VfSat(TS1) * (1.0 - X1) + VgSat(TS1) * X1) - 1.0 / (VfSat(TS2) * (1.0 - X2) + VgSat(TS2) * X2)) / (A - B);
		}
		else
		{
			A = P + 0.00001;
			B = P - 0.00001;
			if (H < VgSat(TSat(A))) A = P;
			DensdP = (1.0 / Vol(T0, A) - 1.0 / Vol(T0, B)) / (A - B);
		}

		return DensdP;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Visg(Tdouble T0)
{
	Tdouble Vis;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Vis = 5.196 * pow(T0, 0.675) * 1e-7;
		return Vis;
	}
	case(MaterialTab::K):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Vis = 4.86372E-6 + 1.15683E-8 * T0;
		return Vis;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Visf(Tdouble T0)
{
	Tdouble Vis;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Vis = 0.0016363 - 4.58269e-6 * T0 + 4.94085e-9 * pow(T0, 2) - 1.84513e-12 * pow(T0, 3);
		return Vis;
	}
	case(MaterialTab::K):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		if (T0 < 653.15 && T0>336.35)
			Vis = 1.131E-5 * Dens(T0) / 3.0 * exp(0.68 * Dens(T0) / T0);
		else if (T0 < 1073.15)
			Vis = 0.0016363 - 4.58269e-6 * T0 + 4.94085e-9 * pow(T0, 2) - 1.84513e-12 * pow(T0, 3);
		else
			Vis = 0, cout << "Properties need to be extended!" << endl;
		return Vis;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Pr(Tdouble T0, Tdouble P0)
{
	Tdouble Pr;
	switch (name) {
	case(MaterialTab::NaK):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Pr = Visf(T0) * 1000.0 * Cp(T0) / Cond(T0);
		return Pr;
	}
	case(MaterialTab::K):
	{
		Pr = 0.8;
		return Pr;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::STen(Tdouble T0)
{
	Tdouble STen;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble T1;
		T1 = 2503.7;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		STen = 240.5*pow((1.0-T0/T1),1.126);
		return STen;
	}
	case(MaterialTab::K):
	{
		STen = 115.7 - 0.064 * (T0 - 273.15);
		return STen;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Compf(Tdouble T0)
{
	Tdouble Compf;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble T1,T2,B;
		B = 3.2682;
		T1 = 2503.7;
		T2 = 371.0;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Compf = 1.717 * (1.0 + ((T0 - T2) / (T1 - T2)) / B) / ((1.0 - (T0 - T2) / (T1 - T2)) * 1.0e4);
		return Compf;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::CpgSat(Tdouble T0)
{
	Tdouble CpgSat;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble A, B, C, D, E, THERMALPRE, VAPORPRE, THERMALEXP, DHDT, AERFAPG, DENDT;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		DENDT = (1.0 / VgSat(T0 + 1.0) - 1.0 / VgSat(T0 - 1.0)) / 2.0;
		THERMALEXP = -VgSat(T0) * DENDT;
		VAPORPRE = (12633.730 / (T0 * T0) - 0.46720 / T0) * exp(11.94630 - 12633.730 / T0 - 0.46720 * log(T0));
		A = 8.353070;
		B = -12905.60;
		C = -0.458240;
		D = 0.00209490;
		E = -0.000000507860;
		THERMALPRE = (-B / pow(T0, 2.0) + C / T0 + D + 2.0 * E * T0) * exp(A + B / T0 + C * log(T0) + D * T0 + E * pow(T0, 2.0));
		AERFAPG = THERMALEXP / (1.0 - VAPORPRE / THERMALPRE);
		DHDT = (HgSat(T0 + 0.10) - HgSat(T0 - 0.10)) / 0.20;
		CpgSat = DHDT + VAPORPRE * VgSat(T0) * (T0 * AERFAPG - 1.0) * 1000.0;
		return CpgSat;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Compg(Tdouble T0)
{
	Tdouble Compg;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble TCR, DENDT, THERMALEXP, VAPORPRE, ERFAPG, DHDT, CP, CV, BEITAT, THERMALPRE, AERFAPG, A, B, C, D, E;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		TCR = 2503.7;
		DENDT = (1.0 / VgSat(T0 + 1.0) - 1.0 / VgSat(T0 - 1.0)) / 2.0;
		THERMALEXP = -VgSat(T0) * DENDT;
		VAPORPRE = (12633.730 / (T0 * T0) - 0.46720 / T0) * exp(11.94630 - 12633.730 / T0 - 0.46720 * log(T0));
		A = 8.353070;
		B = -12905.60;
		C = -0.458240;
		D = 0.00209490;
		E = -0.000000507860;
		THERMALPRE = (-B / pow(T0 ,2) + C / T0 + D + 2 * E * T0) * exp(A + B / T0 + C * log(T0) + D * T0 + E * pow(T0, 2));
		AERFAPG = THERMALEXP / (1 - VAPORPRE / THERMALPRE);
		DHDT = (HgSat(T0 + 1.0) - HgSat(T0 - 1)) / 2.0;
		CP = DHDT + VAPORPRE * VgSat(T0) * (T0 * AERFAPG - 1) * 1000.0;
		CV = CP - T0 * AERFAPG * THERMALPRE * VgSat(T0) * 1000.0;
		BEITAT = AERFAPG / THERMALPRE;
		Compg = BEITAT * CV / CP;
		return Compg;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Soundf(Tdouble T0)
{
	Tdouble Soundf;
	switch (name) {
	case(MaterialTab::NaK):
	{
		Tdouble T1;
		T1 = 2503.7;
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Soundf = 2660.70 - 0.376670 * T0 - 9.03560 * pow(T0, 2.0) / 100000.00;
		return Soundf;
	}
	case(MaterialTab::K):
	{
		if (T0 >= 100 && T0 <= 800)
			Soundf = 1922.3 - 0.539 * T0;
		else
			cout << "Propertied need to be extended!";
		return Soundf;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Expand(Tdouble T0)
{
	Tdouble Expand;
	switch (name) {
	case(MaterialTab::AsFe):
	{
		if (T0 < 371) T0 = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T0 > 2503.07) T0 = 2503.07, cout << "High Temperature in Material. Temperature reset!";
		Expand = 1.7890 / 1.0e5 + 2.3980 * T0 / 1.0e9 + 3.2690 * pow(T0, 2.0) / 1.0e13;
		return Expand;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}

Tdouble Material::Expand(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM)
{
	Tdouble Expand;
	switch (name) {
	case(MaterialTab::MOX):
	{
		if (T < 371) T = 371, cout << "Low Temperature in Material. Temperature reset!";
		if (T > 3120.0) T = 3120.0, cout << "High Temperature in Material. Temperature reset!";
		if (T < 923.0)
			Expand = -2.66e-3 + 9.802e-6 * T - 2.705e-10 * T * T + 4.391e-13 * pow(T, 3.0);
		else if (T < 3120.0)
			Expand = -3.28e-3 + 1.179e-5 * T - 2.429e-9 * T * T + 1.219e-12 * pow(T, 3.0);
		else
			Expand = 0.9285 / (8860.0 - 0.9285 * (T - 3120.0));
		return Expand;
	}
	default:
		cout << "Error: Can not get the material!" << endl;
	}
}