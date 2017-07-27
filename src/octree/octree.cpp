#include "octree.hpp"
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <list>

SubdivisionPos::SubdivisionPos()
{
}

SubdivisionPos::SubdivisionPos(Position center, Position point)
{
	for(int i=0; i<3; i++)
		s[i] = point.x[i] < center.x[i] ? 0 : 1;
}


OctreeNode::OctreeNode(Octree* octree, SubdivisionPos subdivision, OctreeNode* parent) :
	parent(parent),
	subdivisionPos(subdivision),
	subdivisionLevel(parent->subdivisionLevel + 1),
    size(parent->size / 2.0),
    m_octree(octree)
{
	for (int i=0; i<3; i++)
	{
		if (subdivision.s[i] == 0)
			center.x[i] = parent->center.x[i] - size / 2.0;
		else
			center.x[i] = parent->center.x[i] + size / 2.0;
	}
}

OctreeNode::OctreeNode(Octree* octree, Position center, double size) :
        subdivisionLevel(0), center(center), size(size), m_octree(octree)
{
}

void OctreeNode::addElement(const OctreeElement& e)
{
	if (!hasSubnodes)
	{
		if (element == nullptr)
		{
			element.reset(new OctreeElement(e));
			element->parent = this;
            if (m_octree->centerMassUpdatingEnabled())
                updateMassCenterReqursiveUp();
			return;
		}

		if (element->pos == e.pos)
		{
			throw std::runtime_error("Cannot work with 2 elements at one place");
		}
		element->parent = nullptr;
		giveElementToSubnodes(*element);
		element.reset();
	}
	giveElementToSubnodes(e);
}


size_t OctreeNode::elementsCount()
{
	if (!hasSubnodes && element != nullptr)
		return 1;

	size_t count = 0;
	for (int i=0; i<8; i++)
	{
		if (subnodes[i] != nullptr)
			count += subnodes[i]->elementsCount();
	}
	return count;
}

DistToNode OctreeNode::getDistsToNode(Position pos)
{
	DistToNode result;
	if (element != nullptr)
	{
		result.nearest = result.farest = (element->pos - pos).len();
		return result;
	}
	result.nearest = -1;
	result.farest = -1;
	Position corner;
	corner.x[0] = center.x[0] + size/2.0;
	corner.x[1] = center.x[1] + size/2.0;
	corner.x[2] = center.x[2] + size/2.0;
	result.nearest = result.farest = (corner - pos).len();

	for (int x = -1; x <=1; x += 2)
		for (int y = -1; y <=1; y += 2)
			for (int z = -1; z <=1; z += 2)
			{
				corner.x[0] = center.x[0] + x*size/2.0;
				corner.x[1] = center.x[1] + y*size/2.0;
				corner.x[2] = center.x[2] + z*size/2.0;
				double dist = (corner - pos).len();
				if (result.farest < dist)
					result.farest = dist;
				if (result.nearest > dist)
					result.nearest = dist;
			}
	return result;
}

bool OctreeNode::isInside(const Position& pos)
{
    double hs = size / 2.0;
    for (int i=0; i<3; i++)
    {
        if (
            pos.x[i] >= center.x[i] + hs
            || pos.x[i] < center.x[i] - hs
        )
            return false;
    }
    return true;
}

void OctreeNode::dbgOutCoords(std::ostream& s)
{
	for (int x = -1; x <=1; x += 2)
		for (int y = -1; y <=1; y += 2)
			for (int z = -1; z <=1; z += 2)
				s << center[0] + x*size/2.0 << ","
				  << center[1] + y*size/2.0 << ","
				  << center[2] + z*size/2.0
				  << std::endl;

	for (int i=0; i<8; i++)
	{
		if (subnodes[i] != nullptr)
			subnodes[i]->dbgOutCoords(s);
	}
}

void OctreeNode::updateMassCenter()
{
    if (element != nullptr)
    {
        massCenter = element->pos;
        mass = element->value;
        return;
    }
    massCenter = {0.0, 0.0, 0.0};
    mass = 0.0;
    int count = 0;
    for (int i=0; i<8; i++)
    {
        if (subnodes[i] != nullptr)
        {
            double nodeMass = subnodes[i]->mass;
            massCenter += subnodes[i]->massCenter * nodeMass;
            mass += nodeMass;
            count++;
        }
    }
    massCenter /= mass;
}

void OctreeNode::updateMassCenterReqursiveUp()
{
    updateMassCenter();
    if (parent != nullptr)
        parent->updateMassCenterReqursiveUp();
}

void OctreeNode::updateMassCenterReqursiveDown()
{
    /// @todo May be optimized by using one cycle for calls and sum calculations, but this will duplicate code
    for (int i=0; i<8; i++)
        if (subnodes[i] != nullptr)
            subnodes[i]->updateMassCenterReqursiveDown();
    updateMassCenter();
}

void OctreeNode::giveElementToSubnodes(const OctreeElement& e)
{
	SubdivisionPos targerSubdivision(center, e.pos);

	int index = targerSubdivision.index();
	if (subnodes[index] == nullptr)
	{
        subnodes[index].reset(new OctreeNode(m_octree, targerSubdivision, this));
		hasSubnodes = true;
	}
	subnodes[index]->addElement(e);
}

/////////////////////////////////
// Octree
Octree::Octree(Position center, double initialSize) :
	m_center(center),
	m_initialSize(initialSize),
	m_centerIsSet(true)
{
}

Octree::Octree(double initialSize) :
	m_initialSize(initialSize),
	m_centerIsSet(false)
{
}

void Octree::add(const OctreeElement& e)
{
    // Creating foor if no
	if (m_root == nullptr)
	{
		if (!m_centerIsSet)
		{
			m_center = e.pos;
			m_centerIsSet = true;
		}

		m_root.reset(
            new OctreeNode(this, m_center, m_initialSize)
		);
	}
    // Enlarging root cell
    while (!m_root->isInside(e.pos))
    {
        enlargeSpaceIteration(e.pos);
    }
	m_root->addElement(e);
}

size_t Octree::count()
{
	if (m_root != nullptr)
		return m_root->elementsCount();
	else
		return 0;
}

OctreeElement& Octree::getNearest(Position pos)
{
	using namespace std;
	if (m_root == nullptr)
		throw(std::runtime_error("Octree is empty"));
	// ndp = node-distance pair
	using ndp = pair<const OctreeNode*, DistToNode>;
	list<ndp> nodes;
	list<ndp> nodesNext;
	nodes.push_back(ndp(m_root.get(), m_root->getDistsToNode(pos)));

	double minFarest = nodes.front().second.farest;

	do {
		// Finding closes
		for (auto it=nodes.begin(); it!=nodes.end(); it++)
		{
			double farest = it->second.farest;
			if (farest < minFarest)
				minFarest = farest;
		}

		// Removing nodes that are too far
		for (auto it=nodes.begin(); it != nodes.end(); )
		{
			if (it->second.nearest > minFarest)
				it = nodes.erase(it);
			else
				it++;
		}

		nodesNext.clear();
		// Subdivision
		for (auto it=nodes.begin(); it!=nodes.end(); it++)
		{
			const OctreeNode& n = *(it->first);
			if (n.element != nullptr)
			{
				nodesNext.push_back(*it);
				continue;
			}
			for (int i=0; i<8; i++)
			{
				if (n.subnodes[i] == nullptr)
					continue;

				nodesNext.push_back(ndp(n.subnodes[i].get(), n.subnodes[i]->getDistsToNode(pos)));
			}
		}
		swap(nodes, nodesNext);
	} while (!(nodes.size() == 1 && nodes.front().first->element != nullptr));

	// const_cast is not bad, because const modifier used only for code above
	// in this function, and its job is done
	return const_cast<OctreeElement&>(*(nodes.front().first->element));
}

double Octree::mass()
{
    if (m_root == nullptr)
        return 0.0;
    return m_root->mass;
}

const Position& Octree::massCenter()
{
    return m_root->massCenter;
}

void Octree::dbgOutCoords(std::ostream& s)
{
	m_root->dbgOutCoords(s);
}

bool Octree::centerMassUpdatingEnabled() const
{
    return m_centerMassUpdatingEnabled;
}

void Octree::muteCenterMassCalculation()
{
    m_centerMassUpdatingEnabled = false;
}

void Octree::unmuteCenterMassCalculation()
{
    m_centerMassUpdatingEnabled = true;
    /// @todo Optimization: update only if changed
    m_root->updateMassCenterReqursiveDown();
}

void Octree::enlargeSpaceIteration(const Position& p)
{
	Position newRootCenter;
	for (int i=0; i<3; i++)
	{
		double cx = m_root->center.x[i];
		double dcx = m_root->size / 2.0;
        	newRootCenter.x[i] = cx + (p.x[i] > cx ? dcx : -dcx);
	}
	SubdivisionPos subPos(newRootCenter, m_root->center);
    std::unique_ptr<OctreeNode> n(new OctreeNode(this, newRootCenter, m_root->size * 2));
	n->hasSubnodes = true;
	n->subdivisionLevel = m_root->subdivisionLevel - 1;
	n->subdivisionPos = subPos;
	n->subnodes[subPos.index()] = std::move(m_root);
	m_root = std::move(n);
    if (centerMassUpdatingEnabled())
        m_root->updateMassCenter();
}

///////////////////////////
/// CenterMassUpdatingMute

CenterMassUpdatingMute::CenterMassUpdatingMute(Octree& octree) :
    m_octree(octree)
{
    m_octree.muteCenterMassCalculation();
}

CenterMassUpdatingMute::~CenterMassUpdatingMute()
{
    unmute();
}

void CenterMassUpdatingMute::unmute()
{
    if (!m_octree.centerMassUpdatingEnabled())
        m_octree.unmuteCenterMassCalculation();
}
