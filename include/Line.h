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
            Phi,
            None
        };

        enum class IntersectionDirection {
            MIN,
            MAX
        };

        double t_boundary;
        IntersectionType type;
        IntersectionDirection direction;
        uint index;

        Boundary(double t_boundary=0, IntersectionType type=IntersectionType::None, IntersectionDirection direction=IntersectionDirection::MIN,
                        uint index=0) : 
                        t_boundary(t_boundary), type(type), direction(direction), index(index) {}


        bool operator != (const Boundary& boundary) const { 
            return this->type == IntersectionType::R || this->type != boundary.type || this->direction == boundary.direction; 
        }
        bool operator == (const Boundary& boundary) const { 
            return  this->type != IntersectionType::R && this->type == boundary.type && this->direction != boundary.direction; 
        }


        Boundary& operator=(const Boundary& boundary)
        {
            if(&boundary != this)
            {
                this->t_boundary = boundary.t_boundary;
                this->type = boundary.type;
                this->direction = boundary.direction;
                this->index = boundary.index;
            }
            return *this;
        }

    };

    const double rho; // прицельный параметер
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

    const int sign;
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
    uint iR_crit;
    bool crossAxis;

    const double t_epsilon;
    const double t_epsilon_first;
    const double rhoLarmor;

    bool lineWork; // удалось ли прорисовать луч


    inline double getX0L(int l=0) const { return x00 + l*sign*rhoLarmor*cosPhi0; }
    inline double getY0L(int l=0) const { return y00 + l*sign*rhoLarmor*sinPhi0; }

    inline double getXPoint(double t, int l=0) const { return getX0L(l) + t * sx; }
    inline double getYPoint(double t, int l=0) const { return getY0L(l) + t * sy; }
    inline double getZPoint(double t, int l=0) const { return z00 + t * sz; }
    inline double getRPoint(double t, int l=0) const; // sqrt(x**2+y**2)
    inline double getPhiPoint(double t, int l=0) const;
    

    inline double getTR(double r, double tPrevious=-1., int l=0) const;
    inline double getTZ(double z, int l=0) const;
    inline double getTPhi(double phi, int l=0) const;

    void getNewIndex(const Boundary &boundary, uint &iZ, uint &iR, uint &iPhi);


    bool checkBoundary(const std::vector <Boundary> &boundaryArray, const Boundary & boundary, uint index=0) const;
    bool checkBoundaryR(const std::vector <Boundary> &boundaryArray, Boundary::IntersectionDirection direction, uint index=0) const;

    LineData traceLine(uint nz, uint nr, uint nphi, const darray &zArray, 
                        const darray &rArray, const darray &phiArray, 
                        double &tPrevious, uint &iZ, uint &iR, uint &iPhi, std::vector <Boundary> &bPrevious); // сместить луч на одну позицию


    bool createDataArray(uint nz, uint nr, uint nphi, const darray &zArray, 
                        const darray &rArray, const darray &phiArray);


    uint findIndex(const darray &array, uint n, double value, bool left=true) const;

public:
    Line(double rho, double theta, double phi0_rho, double z0_rho, 
        uint nz, uint nr, uint nphi, 
        const darray &zArray, const darray &rArray, const darray &phiArray,
        double t_epsilon=0., double t_epsilon_first=0., bool plusDirection=true,
        double rhoLarmor=0., double rMax=-1.    
    );



    bool getLineWork() const { return lineWork; }
    const std::vector<LineData> & getData() const { return data; }
    uint getNs() const { return ns; }
    const barray & getIntersectionPoints() const { return intersectionPoints; }
};

#endif
