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
	bool isDirty = false;
	bool toRemove = false;

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
public:
    Node(Octree* octree, SubdivisionPos subdivision, Node* parent);
    Node(Octree* octree, Position center, double size);
    void addElement(std::shared_ptr<Element> e);
    size_t elementsCount() const;
    Element& findNearest(Position pos);

    double diameter() const;

    /**
     * @brief Returns minimal and maximal distance to node (to its corners)
     * @param pos Point that distance should be calculated from
     * @return
     */
    DistToNode getDistsToNode(Position pos) const;

	/**
	* @brief Checks if some point is inside this cell
	* @param pos Point to test
	* @return true if inside, false otherwise
	*/
    bool isInside(const Position& pos) const;

	void dbgOutCoords(std::ostream& s);

    std::shared_ptr<Element> element;
    Node* parent = nullptr;
	
	SubdivisionPos subdivisionPos;
	int subdivisionLevel = 0;
	bool hasSubnodes = false;

    //OctreeNode*& getSubnode(const SubdivisionPos& sp);

	Position center;
	double size;

    Position massCenter;
    double mass;

    std::unique_ptr<Node> subnodes[8];

    void updateMassCenterReqursiveUp();
    void updateMassCenterReqursiveDown();
    void updateMassCenter();

private:

    void giveElementToSubnodes(std::shared_ptr<Element> e);
    Octree* m_octree = nullptr;
};

class Octree
{
public:
	Octree(double initialSize = 1.0);
	Octree(Position center, double initialSize = 1.0);
    void clear();
    void add(std::shared_ptr<Element> e);
	void update();
	size_t count();
	
	void dbgOutCoords(std::ostream& s);

    const Element& getNearest(Position pos);
	
    const Node& root() const { return *m_root; }
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

template<typename ResultType = double>
class Convolution
{
public:
    using Visitor = std::function<ResultType(const Position& target, const Position& object, double mass)>;

    Convolution()
    {
        addScale(0.0, 0.0);
    }

    /**
     * @brief addScale Allow averaging with scale averagingScale when objects are farer than minDistance
     * @param minDistance Minimal distance that alow this averaging
     * @param averagingScale Space size of blocks where objects masses may be averaged
     */
    void addScale(double minDistance, double averagingScale)
    {
        m_distsScales.push_back(std::pair<double, double>(minDistance, averagingScale));
    }

    /**
     * @brief Calculate convolution by whole octree without exclusions
     * @param oct       Octree
     * @param target    Point where to calculate
     * @param v         Visitor function
     * @return Result of convolution
     */
    ResultType convolute(const Octree& oct, const Position& target, Visitor v)
    {
        return convoluteExcludingElement(oct, target, v, nullptr);
    }

    /**
     * @brief Calculate convolution by octree excluding element in position of which
     *        convolution should be calculated
     * @param oct            Octree
     * @param excludedTarget Element that should be excluded. In its pos convolution
     *                       will be calculated
     * @param v              Visitor function
     * @return Result of convolution
     */
    ResultType convoluteExcludingElement(const Octree& oct, const Element& excludedTarget, Visitor v)
    {
        return convoluteExcludingElement(oct, excludedTarget.pos, v, &excludedTarget);
    }

private:

    /**
     * @brief Calculate convolution of visitor v by all octree elements. One element may be excluded.
     *        Parameters target and ignoreElement are ionterconnected, see below.
     * @param oct           Octree object
     * @param target        Point where we are calculating convolution
     * @param v             Visitor function
     * @param ignoreElement Element that will be excluded. Important that this element must be
     *                      in one node with point target! This is the reason this function is
     *                      unsafe for user.
     * @return Result of convolution
     */
    ResultType convoluteExcludingElement(const Octree& oct, const Position& target, Visitor v, const Element* ignoreElement = nullptr)
    {
        sortDistsScales();
        m_nodesList.clear();
        ResultType result = ResultType();
        for (
             m_nodesList.push_back(&oct.root());
             m_nodesList.size() != 0;
             m_nodesList.pop_front()
        )
        {
            const Node *n = m_nodesList.front();
            if (n->isInside(target))
            {
                // We are inside this node
                if (n->element != nullptr)
                {
                    // We have only one object in cube with point target, so we can ise directly its mass and center

                    // Check if we must ignore this node
                    if (n->element.get() != ignoreElement)
                    {
                        result += v(target, n->massCenter, n->mass);
                    }
                } else {
                    // Node we located in has subnodes, lets devide it
                    addSubnodesToList(n);
                }
                continue;
            }
            DistToNode d = n->getDistsToNode(target);
            double scale = findScale(d.nearest);
            if (n->diameter() <= scale)
            {
                // We can use averaging over this node
                result += v(target, n->massCenter, n->mass);
            } else {
                // Node is too large, so we should devide it
                addSubnodesToList(n);
            }
        }
        return result;
    }

    void sortDistsScales()
    {
        std::sort(m_distsScales.begin(), m_distsScales.end(),
            [](const std::pair<double, double> p1, std::pair<double, double> p2)
            { return p1.first < p2.first; }
        );
    }

    double findScale(double distance)
    {
        if (distance >= m_distsScales.back().first)
            return m_distsScales.back().second;
        int l = 0, r = m_distsScales.size() - 1;
        int c = (l+r) / 2;
        while (r-l > 1)
        {
            if (m_distsScales[c].first <= distance)
                l = c;
            else
                r = c;
            c = (l+r) / 2;
        }
        return m_distsScales[l].second;
    }

    void addSubnodesToList(const Node* n)
    {
        for (int i=0; i<8; i++)
        {
            const Node *subnode = n->subnodes[i].get();
            if (subnode != nullptr)
                m_nodesList.push_back(subnode);
        }
    }
    std::vector<std::pair<double, double>> m_distsScales;
    std::list<const Node*> m_nodesList;
};

}

#endif // LIBHEADER_INCLUDED
