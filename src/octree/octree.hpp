#ifndef LIBHEADER_INCLUDED
#define LIBHEADER_INCLUDED

#include "geom-vector.hpp"

#include <ostream>
#include <functional>
#include <vector>
#include <list>

#include <memory>
#include <cmath>

class OctreeNode;
class Octree;

struct OctreeElement
{
	OctreeElement(const Position& p, double value = 0.0) :
		pos(p),
		value(value)
	{ }
	OctreeElement(double x = 0.0, double y = 0.0, double z = 0.0, double value = 0.0) :
		pos(x, y, z),
		value(value)
	{ }
	Position pos;
	double value = 0.0;
	bool isDirty = false;
	bool toRemove = false;

	OctreeNode* parent = nullptr;
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

class OctreeNode
{
public:
    OctreeNode(Octree* octree, SubdivisionPos subdivision, OctreeNode* parent);
    OctreeNode(Octree* octree, Position center, double size);
	void addElement(const OctreeElement& e);
    size_t elementsCount() const;
	OctreeElement& findNearest(Position pos);

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

	std::unique_ptr<OctreeElement> element;
	OctreeNode* parent = nullptr;
	
	SubdivisionPos subdivisionPos;
	int subdivisionLevel = 0;
	bool hasSubnodes = false;

    //OctreeNode*& getSubnode(const SubdivisionPos& sp);

	Position center;
	double size;

    Position massCenter;
    double mass;

	std::unique_ptr<OctreeNode> subnodes[8];

    void updateMassCenterReqursiveUp();
    void updateMassCenterReqursiveDown();
    void updateMassCenter();

private:

	void giveElementToSubnodes(const OctreeElement& e);
    Octree* m_octree = nullptr;
};

class Octree
{
public:
	Octree(double initialSize = 1.0);
	Octree(Position center, double initialSize = 1.0);
	void add(const OctreeElement& e);
	void update();
	size_t count();
	
	void dbgOutCoords(std::ostream& s);

	OctreeElement& getNearest(Position pos);
	
    const OctreeNode& root() const { return *m_root; }
    double mass();
    const Position& massCenter();

    bool centerMassUpdatingEnabled() const;
    void muteCenterMassCalculation();
    void unmuteCenterMassCalculation();

private:
	void enlargeSpaceIteration(const Position& p);
	bool isPointInsideRoot(const Position& p);

	std::unique_ptr<OctreeNode> m_root;
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

class Convolution
{
public:
    using Visitor = std::function<double(const Position& target, const Position& object, double mass)>;

    Convolution();
    /**
     * @brief addScale Allow averaging with scale averagingScale when objects are farer than minDistance
     * @param minDistance Minimal distance that alow this averaging
     * @param averagingScale Space size of blocks where objects masses may be averaged
     */
    void addScale(double minDistance, double averagingScale);

    /**
     * @brief Calculate convolution of visitor v by all octree elements
     * @param oct Octree object
     * @param target Point where we are calculating convolution
     * @param v Visitor function
     * @return result of convolution
     */
    double convolute(const Octree& oct, const Position& target, Visitor v);
    void sortDistsScales();
    double findScale(double distance);
private:
    void addSubnodesToList(const OctreeNode* n);
    std::vector<std::pair<double, double>> m_distsScales;
    std::list<const OctreeNode*> m_nodesList;
};

#endif // LIBHEADER_INCLUDED
