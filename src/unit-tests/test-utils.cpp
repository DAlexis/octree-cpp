#include "test-utils.hpp"

void PointsGenerator::addGrid(int n, double size, Octree& oct, std::vector<Position>* positions)
{
    CenterMassUpdatingMute m(oct);
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int k=0; k<n; k++)
            {
                Position p(
                    -size/2.0 + size / (n-1) * i,
                    -size/2.0 + size / (n-1) * j,
                    -size/2.0 + size / (n-1) * k
                );
                oct.add(ElementValue(p, 1.0));
                if (positions)
                    positions->push_back(p);
            }
}

Position& PointsGenerator::findNearestBruteForce(const Position& pos, std::vector<Position>& positions)
{
    double distMin = (positions.front() - pos).len();
    Position* p = &positions.front();
    for (auto it = positions.begin(); it != positions.end(); ++it)
    {
        double d = (*it - pos).len();
        if (d < distMin)
        {
            distMin = d;
            p = &(*it);
        }
    }
    return *p;
}
