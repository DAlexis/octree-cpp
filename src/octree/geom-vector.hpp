/*
 * geom-vector.hpp
 *
 *  Created on: 7 июн. 2017 г.
 *      Author: dalexies
 */

#ifndef OCTREE_GEOM_VECTOR_HPP_
#define OCTREE_GEOM_VECTOR_HPP_

#include <string>
#include <initializer_list>
#include <cstring>
#include <cmath>

template<typename T>
T sqr(T x)
{
	return x*x;
}

/**
 * This is a general-purpose vector class. You may use it in your project and may not.
 * If octree library needs GeomVector as argument, you may pass raw double array, so
 * it will not conflict with your linear algebra libraries
 */
template<int dim>
class GeomVector
{
public:

    GeomVector(const double* coords)
    {
        memcpy(x, coords, sizeof(double)*dim);
    }

    GeomVector(std::initializer_list<double> initList)
    {
        double* px = this->x;
        for (const double& coord : initList)
            *(px++) = coord;
    }

    GeomVector(double x, double y, double z)
    	{ this->x[0] = x; this->x[1] = y; this->x[2] = z; }

    GeomVector(double x, double y)
    	{ this->x[0] = x; this->x[1] = y; }

    GeomVector(double x)
    	{ this->x[0] = x; }


    GeomVector()
    	{ memset(x, 0, sizeof(x[0])*dim); }

    double len() const
    {
    	double result = 0;
    	for (int i=0; i<dim; i++)
    	{
    		result += sqr(x[i]);
    	}
    	return std::sqrt(result);
    }

    GeomVector<dim> operator=(std::initializer_list<double> initList)
    {
        double* px = this->x;
        for (const double& coord : initList)
            *(px++) = coord;
        return *this;
    }

    /// @todo optimize this functions. Remove this cycles for concrete dimensions
    GeomVector<dim> operator-() const
    {
        GeomVector<dim> result;
        for (int i=0; i<dim; ++i)
            result.x[i] = -x[i];
        return result;
    }

    GeomVector<dim> operator-(const GeomVector<dim>& right) const
    {
        GeomVector<dim> result;
        for (int i=0; i<dim; ++i)
            result.x[i] = x[i] - right.x[i];
        return result;
    }

    GeomVector<dim> operator+(const GeomVector<dim>& right) const
    {
        GeomVector<dim> result;
        for (int i=0; i<dim; ++i)
            result.x[i] = x[i] + right.x[i];
        return result;
    }

    GeomVector<dim> operator+=(const GeomVector<dim>& right)
    {
        for (int i=0; i<dim; ++i)
            x[i] += right.x[i];
        return *this;
    }

    GeomVector<dim> operator-=(const GeomVector<dim>& right)
    {
        for (int i=0; i<dim; ++i)
            x[i] -= right.x[i];
        return *this;
    }

    GeomVector<dim> operator*=(double right)
    {
        for (int i=0; i<dim; ++i)
            x[i] *= right;
        return *this;
    }

    GeomVector<dim> operator/=(double right)
    {
        for (int i=0; i<dim; ++i)
            x[i] /= right;
        return *this;
    }

    bool operator==(const GeomVector<dim>& right) const
    {
        for (int i = 0; i<dim; ++i)
            if (x[i] != right.x[i])
                return false;
        return true;
    }

    bool operator!=(const GeomVector<dim>& right) const
    {
        return ! (*this == right);
    }

    GeomVector<dim> operator*(double right) const
    {
        GeomVector<dim> result;
        for (int i=0; i<dim; ++i)
            result.x[i] = x[i] * right;
        return result;
    }

    /**
     * @brief Scalar product
     */
    double operator*(const GeomVector<dim>& right) const
    {
    	double result = 0;
    	for (int i=0; i<dim; i++)
    		result += x[i] * right.x[i];

    	return result;
    }

    /**
     * @brief Vector product
     */
    GeomVector<dim> operator%(const GeomVector<dim>& right) const
    {
        static_assert(dim != 3, "Vector product works only for dimension = 3!");
    	GeomVector<dim> result;
    	result[0] =   (*this)[1] * right[2] - (*this)[2] * right[1];
		result[1] = - (*this)[0] * right[2] + (*this)[2] * right[0];
		result[2] =   (*this)[0] * right[1] - (*this)[1] * right[0];
		return result;
    }

    GeomVector<dim> operator/(double right)
    {
        GeomVector<dim> result;
        for (int i=0; i<dim; ++i)
            result.x[i] = x[i] / right;
        return result;
    }

    std::string str() const
    {
    	return std::string("(") +
    			std::to_string(x[0]) + "; " +
				std::to_string(x[1]) + "; " +
				std::to_string(x[2]) + ")";
    }

    double& operator[](unsigned int i) { return x[i]; }

    const double& operator[](unsigned int i) const { return x[i]; }

    void normalize()
    {
    	double l = len();
    	for (int i=0; i<dim; i++)
    		x[i] /= l;
    }

    double distTo(const GeomVector<dim>& target) const
    {
        double result = 0.0;
        const double *t = target.x;
        for (int i=0; i<dim; i++)
        {
            double d = x[i] - t[i];
            result += d*d;
        }
        return sqrt(result);
    }

    double x[dim];
};

using Position = GeomVector<3>;

#endif /* OCTREE_GEOM_VECTOR_HPP_ */
