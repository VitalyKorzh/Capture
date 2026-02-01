#ifndef __LINE_H__
#define __LINE_H__

#include <vector>

typedef std::vector<double> darray;
typedef std::vector<unsigned> uiarray;
typedef std::vector<bool> barray;
typedef unsigned uint;


struct LineData
{
    uint iZ;
    uint iR;
    uint iPhi;
    double s;
    LineData(uint iZ=0, uint iR=0, uint iPhi=0, double s=0) : 
    iZ(iZ), iR(iR), iPhi(iPhi), s(s)
    {}

};

class Line
{
private:

    struct Boundary
    {
        enum class IntersectionType {
            R,
            Z,
            Phi
        };

        double t_boundary;
        IntersectionType type;
        uint index;

        Boundary(double t_boundary, IntersectionType type, uint index) : 
                        t_boundary(t_boundary), type(type), index(index) {}

    };

    const double rho; // прицельный параметер
    const double rho2;
    const double theta; // питч-угол
    const double phi0_rho; // угол при r = rho
    const double z0_rho; // растояние по оси z при r = rho

    const double sinTheta;
    const double cosTheta;
    const double sinPhi0;
    const double cosPhi0;

    double d; // sqrt(R0**2-rho**2)

    double x00; // начальная точка луча x
    double y00; // начальная точка луча y
    double z00; // начальная точка луча z

    const double sx; // направление луча x -sin(theta)*sin(phi0_rho)
    const double sy; // направление луча y sin(theta)*cos(phi0_rh0)
    const double sz; // направление луча z cos(theta)


    uint ns;
    barray intersectionPoints; // точки пересечения линии
    std::vector <LineData> data; // положение линии

    uint nz;
    uint nr;
    uint nphi;
    double Rmax;
    double zmin;
    double zmax;

    double t_crit;
    bool crossAxis;

    const double t_epsilon;
    bool lineWork; // удалось ли прорисовать луч

    inline double getXPoint(double t) const { return x00 + t * sx; }
    inline double getYPoint(double t) const { return y00 + t * sy; }
    inline double getZPoint(double t) const { return z00 + t * sz; }
    inline double getRPoint(double t) const; // sqrt(x**2+y**2)
    inline double getPhiPoint(double t) const;
    

    inline double getTR(double r, double tPrevious=-1.) const;
    inline double getTZ(double z) const;
    inline double getTPhi(double phi) const;

    void getNewIndex(const Boundary &boundary, uint &iZ, uint &iR, uint &iPhi);

    LineData traceLine(uint nz, uint nr, uint nphi, const darray &zArray, 
                        const darray &rArray, const darray &phiArray, 
                        double &tPrevious, uint &iZ, uint &iR, uint &iPhi); // сместить луч на одну позицию


    bool createDataArray(uint nz, uint nr, uint nphi, const darray &zArray, 
                        const darray &rArray, const darray &phiArray);

public:
    Line(double rho, double theta, double phi0_rho, double z0_rho, 
        uint nz, uint nr, uint nphi, 
        const darray &zArray, const darray &rArray, const darray &phiArray,
        double t_epsilon=0.      
    );



    bool getLineWork() const { return lineWork; }
    const std::vector<LineData> & getData() const { return data; }
    uint getNs() const { return ns; }
    const barray & getIntersectionPoints() const { return intersectionPoints; }
};

#endif
