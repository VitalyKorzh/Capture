#ifndef __PHYSIC_VALUES_H__
#define __PHYSIC_VALUES_H__

struct  PhysicValues
{
public:
    static constexpr double ee = 4.8e-10; // элементарный заряд СГСЭ
    static constexpr double ME = 9.11e-28; // масса электрона в г
    static constexpr double MP = 1.672e-24; // масса протона в г
    static constexpr double C = 2.998e10; // скорость света в см/c
    static constexpr double EV_TO_ERG = 1.6e-12; // коэффициент перевода из эВ в эрг
    static constexpr double KEV_TO_ERG = 1.6e-9; // коэффициент перевода из кэВ в эрг
    static constexpr double T_TO_GS = 1e4; // коэффициент перевода из Т в Гс
    static constexpr double A_TO_STAT_AMPER = 3e9; // скольк кл в секунду в одном стат ампере
    static constexpr double E_A_TO_P_TO_S = 6.24e18; // перевод в экв.А. в шт/c
    static constexpr double KWT_TO_ERG_S = 1e10; // перевод мощности из кВт в эрг/c
    static constexpr double BARN_TO_CM2 = 1e-24; // перевод барн в см^2
    virtual ~PhysicValues()=0;
};


#endif