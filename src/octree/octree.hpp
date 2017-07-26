#ifndef LIBHEADER_INCLUDED
#define LIBHEADER_INCLUDED

#include "geom-vector.hpp"

#include <ostream>

#include <memory>
#include <cmath>

class OctreeNode;

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
	OctreeNode(SubdivisionPos subdivision, OctreeNode* parent);
	OctreeNode(Position center, double size);
	void addElement(const OctreeElement& e);
	size_t elementsCount();
	OctreeElement& findNearest(Position pos);

	/**
	 * Returns minimal and maximal distance to node (to its corners)
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

	std::unique_ptr<OctreeNode> subnodes[8];

private:
	void giveElementToSubnodes(const OctreeElement& e);
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

private:
	void enlargeSpaceIteration(const Position& p);
	bool isPointInsideRoot(const Position& p);

	std::unique_ptr<OctreeNode> m_root;
	Position m_center;
	double m_initialSize;
	bool m_centerIsSet;
};

#endif // LIBHEADER_INCLUDED
