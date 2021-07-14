#ifndef MATERIAL_H
#define MATERIAL_H 1

#include <string>
#include <iostream>

using namespace std;

typedef long double Tdouble;
typedef int Tint;

enum class MaterialTab {NaK, AsFe, UO2, B4C, MOX, MoNb, Mo, ZrH, BeO, Be, CO2, He, Xe, Cs, Wick, K};
MaterialTab MaterialFlag(string flag);

//ƽ������δ���
class Material
{
private:
	Material(const Material& mat);

protected:
	MaterialTab name;

public:
	Material() {};
	Material(MaterialTab name0);
	~Material();

	void SetMat(MaterialTab name0) { name = name0; };

	Tdouble T(Tdouble H, Tdouble P);
	Tdouble TSat(Tdouble P);
	
	Tdouble PfSat(Tdouble T);
	Tdouble PgSat(Tdouble T, Tdouble D);


	//����
	Tdouble Vol(Tdouble T, Tdouble P);
	Tdouble VfSat(Tdouble T);
	Tdouble VgSat(Tdouble T);

	//��
	Tdouble H(Tdouble T, Tdouble P);
	Tdouble HfSat(Tdouble T);
	Tdouble HgSat(Tdouble T);
	//����Ǳ��
	Tdouble Hfg(Tdouble T);


	//�ܶ�
	Tdouble Dens(Tdouble H, Tdouble P);
	Tdouble Dens(Tdouble T);
	Tdouble Dens(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);
	Tdouble DensdT(Tdouble T);
	Tdouble DensdH(Tdouble H, Tdouble P);
	Tdouble DensdP(Tdouble H, Tdouble P);

	//����ϵ��
	Tdouble Cond(Tdouble T);
	Tdouble Cond(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);
	//���履�͵���ϵ��
	Tdouble CondSat(Tdouble T);

	//������
	Tdouble Cp(Tdouble T);
	Tdouble Cp(Tdouble T, Tdouble Y);
	//Һ�履�ͱ�����
	Tdouble CpfSat(Tdouble T);
	//���履�ͱ�����
	Tdouble CpgSat(Tdouble T);
	//�������
	Tdouble VCp(Tdouble T);


	//ճ��
	//���嶯��ճ��
	Tdouble Visg(Tdouble T);
	//Һ�嶯��ճ��
	Tdouble Visf(Tdouble T);
	
	//��������
	Tdouble Pr(Tdouble T, Tdouble P);

	//��������
	Tdouble STen(Tdouble T);

	//����Һ̬��ѹ��ϵ��
	Tdouble Compf(Tdouble T);
	//������̬��ѹ��ϵ��
	Tdouble Compg(Tdouble T);

	//Һ̬����
	Tdouble Soundf(Tdouble T);

	//������ϵ��
	Tdouble Expand(Tdouble T);
	Tdouble Expand(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);

};



#endif
