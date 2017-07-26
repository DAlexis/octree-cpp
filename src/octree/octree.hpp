#ifndef LIBHEADER_INCLUDED
#define LIBHEADER_INCLUDED

#include "geom-vector.hpp"

#include <ostream>

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
	size_t elementsCount();
	OctreeElement& findNearest(Position pos);

    /**
     * @brief Returns minimal and maximal distance to node (to its corners)
     * @param pos Point that distance should be calculated from
     * @return
     */
	DistToNode getDistsToNode(Position pos);

	/**
	* @brief Checks if some point is inside this cell
	* @param pos Point to test
	* @return true if inside, false otherwise
	*/
	bool isInside(const Position& pos);

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
    Position mass;

	std::unique_ptr<OctreeNode> subnodes[8];

    void updateMassCenterReqursiveUp();
    void updateMassCenterReqursiveDown();

private:
    void updateMassCenter();

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
	
	const OctreeNode& root() { return *m_root; }

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

private:
};

#endif // LIBHEADER_INCLUDED
