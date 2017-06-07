#ifndef LIBHEADER_INCLUDED
#define LIBHEADER_INCLUDED

#include <memory>

struct Position
{
	Position(double x = 0.0, double y = 0.0, double z = 0.0)
		{ this->x[0] = x; this->x[1] = y; this->x[2] = z; }
	double x[3] = {0.0, 0.0, 0.0};
};

struct OctreeElement
{
	OctreeElement(double x = 0.0, double y = 0.0, double z = 0.0, double value = 0.0) :
		pos(x, y, z),
		value(value)
	{ }
	Position pos;
	double value = 0.0;
	bool isDirty = false;
	bool toRemove = false;
};

struct SubdivisionPos
{
	constexpr static unsigned char npos = 4;
	unsigned char s[3] = {npos, npos, npos};
	unsigned char index()
	{
		return s[0]+2*s[1]+4*s[2];
	}
};

class OctreeNode
{
public:
	OctreeNode(SubdivisionPos subdivision, OctreeNode* parent);
	OctreeNode(Position cornerMin, Position cornerMax);

	void addElement(const OctreeElement& e);
	size_t elementsCount();

	std::unique_ptr<OctreeElement> element;
	OctreeNode* parent = nullptr;
	
	SubdivisionPos subdivisionPos;
	unsigned int subdivisionLevel = 0;
	bool hasSubnodes = false;
	
	//OctreeNode*& getSubnode(const SubdivisionPos& sp);

	//OctreeNode* getNeighbour();

	Position cornerMin;
	Position cornerMax;
	Position center;

	std::unique_ptr<OctreeNode> subnodes[8];

private:
	void calculateCenter();
	void giveElementToSubnodes(const OctreeElement& e);

};

class Octree
{
public:
	Octree(Position minCorner, Position maxCorner);
	Octree(double initialSize = 1.0);
	void add(const OctreeElement& e);
	void update();
	size_t count();
	
	OctreeElement& getNearest(Position pos);
	
private:
	std::unique_ptr<OctreeNode> m_root;
	bool m_initializedByCorners = false;
	Position m_minCorner, m_maxCorner;
	double m_initialSize;
};

#endif // LIBHEADER_INCLUDED
