#include "octree.hpp"
#include <iostream>
#include <cstring>
#include <stdexcept>

OctreeNode::OctreeNode(SubdivisionPos subdivision, OctreeNode* parent) :
	parent(parent),
	subdivisionPos(subdivision)
{
	if (parent != nullptr)
	{
		cornerMin = parent->cornerMin;
		cornerMax = parent->cornerMax;
		for (int i=0; i<3; i++)
		{
			if (subdivision.s[0] == 0)
				cornerMax.x[0] = parent->center.x[0];
			else
				cornerMin.x[0] = parent->center.x[0];
		}
	}
	calculateCenter();
}

OctreeNode::OctreeNode(Position cornerMin, Position cornerMax) :
	parent(nullptr),
	cornerMin(cornerMin),
	cornerMax(cornerMax)
{
	calculateCenter();
}

void OctreeNode::addElement(const OctreeElement& e)
{
	if (!hasSubnodes)
	{
		if (element != nullptr)
		{
			giveElementToSubnodes(*element);
			element.reset();
		} else {
			element.reset(new OctreeElement(e));
			return;
		}
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

void OctreeNode::giveElementToSubnodes(const OctreeElement& e)
{
	SubdivisionPos targerSubdivision;
	for (int i=0; i<3; i++)
	{
		if (e.pos.x[i] < center.x[i])
			targerSubdivision.s[i] = 0;
		else
			targerSubdivision.s[i] = 1;
	}
	int index = targerSubdivision.index();
	if (subnodes[index] == nullptr)
	{
		subnodes[index].reset(new OctreeNode(targerSubdivision, this));
		hasSubnodes = true;
	}
	subnodes[index]->addElement(e);
}

void OctreeNode::calculateCenter()
{
	for (int i=0; i<3; i++)
		center.x[i] = (cornerMin.x[i] + cornerMax.x[i]) / 2.0;
}

/////////////////////////////////
// Octree
Octree::Octree(Position minCorner, Position maxCorner) :
	m_initializedByCorners(true),
	m_minCorner(minCorner),
	m_maxCorner(maxCorner),
	m_initialSize(0.0)
{
}

Octree::Octree(double initialSize) :
	m_initializedByCorners(false),
	m_initialSize(0.0)
{
}

void Octree::add(const OctreeElement& e)
{
	if (m_root == nullptr)
	{
		if (!m_initializedByCorners)
		{
			for (int i=0; i<3; i++)
			{
				m_minCorner.x[i] = e.pos.x[i] - m_initialSize / 2.0;
				m_maxCorner.x[i] = e.pos.x[i] + m_initialSize / 2.0;
			}
		} else {
			for (int i=0; i<3; i++)
			{
				if (m_minCorner.x[i] > e.pos.x[i] || m_maxCorner.x[i] < e.pos.x[i])
					throw std::runtime_error("Space extension not supported yet");
			}
		}

		m_root.reset(
			new OctreeNode(
				m_minCorner, m_maxCorner
			)
		);
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
