#ifndef OCTREE_HPP_INCLUDED
#define OCTREE_HPP_INCLUDED

#include "geom-vector.hpp"

#include <ostream>
#include <functional>
#include <vector>
#include <list>

#include <memory>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace octree {

class Node;
class Octree;

/**
 * @brief Octree element with reference to its value
 */
struct Element
{
    virtual ~Element() {}
    Element(const Position& p, double& value) :
		pos(p),
        value(value)
	{ }
    Element(double x, double y, double z, double& value) :
		pos(x, y, z),
		value(value)
	{ }
	Position pos;
    double &value;

    Node* parent = nullptr;
};

/**
 * @brief Octree element storing its value inside
 */
struct ElementValue : public Element
{
    ElementValue(const Position& p, double value = 0.0) :
        Element(p, storedValue),
        storedValue(value)
    {}
    ElementValue(double x, double y, double z, double value = 0.0) :
        Element(x, y, z, storedValue),
        storedValue(value)
    {}
    double storedValue;
};

struct SubdivisionPos
{
	SubdivisionPos();
	SubdivisionPos(Position center, Position point);

	constexpr static unsigned char npos = 4;
	unsigned char s[3] = {npos, npos, npos};
	unsigned char index()
	{
		return s[0]+2*s[1]+4*s[2];
	}
};

struct DistToNode
{
    double nearest = 0.0, farest = 0.0;
};

/**
 * @brief The octree Node class.
 * Node has 3 states:
 *  - Empty =>           Element == nullptr, subnodes[i] == nullptr
 *  - Holds 1 element => Element != nullptr, subnodes[i] == nullptr
 *  - Holds subnodes =>  Element == nullptr, subnodes[i] != nullptr
 */
class Node
{
friend class Octree;
public:
    Node(Octree* octree, SubdivisionPos subdivision, Node* parent);
    Node(Octree* octree, Position center, double size);
    void addElement(std::shared_ptr<Element> e);

    size_t elementsCount() const;

    /**
     * @brief Returns minimal and maximal distance to node (to its corners)
     * @param pos Point that distance should be calculated from
     * @return
     */
    DistToNode getDistsToNode(Position pos) const;

    double getMinDist(const Position& pos) const;

    double getDistToCenter(const Position& pos) const;

	/**
	* @brief Checks if some point is inside this cell
	* @param pos Point to test
	* @return true if inside, false otherwise
	*/
    bool isInside(const Position& pos) const;

    void dbgOutCoords(std::ostream& s) const;

    std::shared_ptr<Element> element;
	
	SubdivisionPos subdivisionPos;

	Position center;
	double size;

    double dia;

    Position massCenter;
    double mass;

    std::unique_ptr<Node> subnodes[8];

    void updateMassCenterReqursiveUp();
    void updateMassCenterReqursiveDown();
    void updateMassCenter();

    /**
     * @brief Put all non-zero const pointers to subnodes into container
     * @param container Any container supporting push_back method
     */
    template<class T>
    void pushBackSubnodes(T& container) const
    {
        for (int i=0; i<8; i++)
        {
            const Node *subnode = subnodes[i].get();
            if (subnode != nullptr)
                container.push_back(subnode);
        }
    }

    template<class T>
    void pushBackAllElements(T& container) const
    {
        if (element != nullptr)
        {
            container.push_back(element.get());
            return;
        }

        for (int i=0; i<8; i++)
        {
            const Node *subnode = subnodes[i].get();
            if (subnode != nullptr)
                subnode->pushBackAllElements(container);
        }
    }

private:
    int subdivisionLevel = 0;
    bool hasSubnodes = false;
    Node* parent = nullptr;

    void giveElementToSubnodes(std::shared_ptr<Element> e);
    void calculateCorners();
    void updateDiameter();

    Octree* m_octree = nullptr;
    Position m_corners[8];
};

class Octree
{
public:
	Octree(double initialSize = 1.0);
	Octree(Position center, double initialSize = 1.0);
    void clear();
    bool empty() const;
    void add(std::shared_ptr<Element> e);
	void update();
	size_t count();
	
	void dbgOutCoords(std::ostream& s);

    const Element& getNearest(Position pos);
    void getClose(std::vector<Element*>& target, const Position& pos, double dist) const;
	
    const Node& root() const;
    double mass();
    const Position& massCenter();

    bool centerMassUpdatingEnabled() const;
    void muteCenterMassCalculation();
    void unmuteCenterMassCalculation();

private:
	void enlargeSpaceIteration(const Position& p);
	bool isPointInsideRoot(const Position& p);

    std::unique_ptr<Node> m_root;
	Position m_center;
	double m_initialSize;
	bool m_centerIsSet;
    bool m_centerMassUpdatingEnabled = true;
};

/**
 * @brief The CenterMassUpdatingMute class
 * RAII object to mute center mass calculation in octree
 */
class CenterMassUpdatingMute
{
public:
    CenterMassUpdatingMute(Octree& octree);
    ~CenterMassUpdatingMute();

    void unmute();

private:
    Octree& m_octree;
};

class IScalesConfig
{
public:
    virtual ~IScalesConfig() {}
    virtual double findScale(double distance) const = 0;
};

class LinearScales : public IScalesConfig
{
public:
    LinearScales(double k=0.5);
    double findScale(double distance) const override;
private:
    double m_k;
};

/**
 * @brief The ScalesConfig class stores averaging scales and
 * corresponding minimal distances
 */
class DiscreteScales : public IScalesConfig
{
public:
    DiscreteScales();
    /**
     * @brief addScale Allow averaging with scale averagingScale when objects are farer than minDistance
     * @param minDistance Minimal distance that alow this averaging
     * @param averagingScale Space size of blocks where objects masses may be averaged
     */
    void addScale(double minDistance, double averagingScale);

    double findScale(double distance) const override;

private:
    void sortDistsScales();
    std::vector<std::pair<double, double>> m_distsScales;
};

template<typename ResultType = double>
class Convolution
{
public:
    using Visitor = std::function<ResultType(const Position& target, const Position& object, double mass)>;

    Convolution(const IScalesConfig& scalesConfig) :
        m_scalesConfig(scalesConfig)
    {
    }

    /**
     * @brief Calculate convolution by whole octree without exclusions
     * Algorythm is upgraded. No isInside check used
     * @param oct       Octree
     * @param target    Point where to calculate
     * @param v         Visitor function
     * @return Result of convolution
     */
    ResultType convolute(const Octree& oct, const Position& target, Visitor v)
    {
        // Vector is used instead of list to prevent new/deletes for single pointers
        std::vector<const Node*> nodesVector;
        nodesVector.reserve(200);
        ResultType result = ResultType();
        if (oct.empty())
            return result;

        nodesVector.push_back(&oct.root());
        for (size_t i=0; i != nodesVector.size(); i++)
        {
            const Node *n = nodesVector[i];

            // This variant approximate a cube by a sphere and it is faster,
            // because it does not contain any ifs and min/max finding
            double dia = n->dia;
            double dist = n->getDistToCenter(target) - dia * 0.5;
            double scale = m_scalesConfig.findScale(dist);
            if (dia <= scale)
            {
                // We can use averaging over this node
                result += v(target, n->massCenter, n->mass);
            } else {
                // Node is too large, so we should devide it
                n->pushBackSubnodes(nodesVector);
            }
        }
        return result;
    }

private:

    const IScalesConfig& m_scalesConfig;
};

}

#endif // LIBHEADER_INCLUDED
