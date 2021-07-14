#ifndef MATERIAL_H
#define MATERIAL_H 1

#include <string>
#include <iostream>

using namespace std;

typedef long double Tdouble;
typedef int Tint;

enum class MaterialTab {NaK, AsFe, UO2, B4C, MOX, MoNb, Mo, ZrH, BeO, Be, CO2, He, Xe, Cs, Wick, K};
MaterialTab MaterialFlag(string flag);

//平滑导数未添加
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


	//比容
	Tdouble Vol(Tdouble T, Tdouble P);
	Tdouble VfSat(Tdouble T);
	Tdouble VgSat(Tdouble T);

	//焓
	Tdouble H(Tdouble T, Tdouble P);
	Tdouble HfSat(Tdouble T);
	Tdouble HgSat(Tdouble T);
	//蒸发潜热
	Tdouble Hfg(Tdouble T);


	//密度
	Tdouble Dens(Tdouble H, Tdouble P);
	Tdouble Dens(Tdouble T);
	Tdouble Dens(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);
	Tdouble DensdT(Tdouble T);
	Tdouble DensdH(Tdouble H, Tdouble P);
	Tdouble DensdP(Tdouble H, Tdouble P);

	//导热系数
	Tdouble Cond(Tdouble T);
	Tdouble Cond(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);
	//气体饱和导热系数
	Tdouble CondSat(Tdouble T);

	//比热容
	Tdouble Cp(Tdouble T);
	Tdouble Cp(Tdouble T, Tdouble Y);
	//液体饱和比热容
	Tdouble CpfSat(Tdouble T);
	//气体饱和比热容
	Tdouble CpgSat(Tdouble T);
	//体积热容
	Tdouble VCp(Tdouble T);


	//粘度
	//气体动力粘度
	Tdouble Visg(Tdouble T);
	//液体动力粘度
	Tdouble Visf(Tdouble T);
	
	//普朗特数
	Tdouble Pr(Tdouble T, Tdouble P);

	//表面张力
	Tdouble STen(Tdouble T);

	//饱和液态可压缩系数
	Tdouble Compf(Tdouble T);
	//饱和气态可压缩系数
	Tdouble Compg(Tdouble T);

	//液态声速
	Tdouble Soundf(Tdouble T);

	//热膨胀系数
	Tdouble Expand(Tdouble T);
	Tdouble Expand(Tdouble T, Tdouble P, Tdouble Y, Tdouble OM);

};



#endif
