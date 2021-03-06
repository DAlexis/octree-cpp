#include "octree.hpp"
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <list>
#include <algorithm>

using namespace octree;

SubdivisionPos::SubdivisionPos()
{
}

SubdivisionPos::SubdivisionPos(Position center, Position point)
{
    for(int i=0; i<3; i++)
        s[i] = point.x[i] < center.x[i] ? 0 : 1;
}


Node::Node(Octree* octree, SubdivisionPos subdivision, Node* parent) :
    subdivisionPos(subdivision),
    size(parent->size * 0.5),
    subdivisionLevel(parent->subdivisionLevel + 1),
    parent(parent),
    m_octree(octree)

{
    double hs = size * 0.5;
    for (int i=0; i<3; i++)
    {
        if (subdivision.s[i] == 0)
            center.x[i] = parent->center.x[i] - hs;
        else
            center.x[i] = parent->center.x[i] + hs;
    }
    calculateCorners();
    updateDiameter();
}

Node::Node(Octree* octree, Position center, double size) :
        center(center), size(size), subdivisionLevel(0), m_octree(octree)
{
    calculateCorners();
    updateDiameter();
}

void Node::addElement(std::shared_ptr<Element> e)
{
    if (!hasSubnodes)
    {
        // If this node is empty, adding rlement directly here
        if (element == nullptr)
        {
            element = e;
            element->parent = this;
            if (m_octree->centerMassUpdatingEnabled())
                updateMassCenterReqursiveUp();
            updateDiameter();
            return;
        }

        // So this node is not empty, but it has no subnodes
        // and hods element by itself. Giving holded element to subnodes
        // and giving element e to subnodes too. This node became
        // subnodes-holding
        if (element->pos == e->pos)
        {
            updateDiameter();
            throw std::runtime_error("Cannot work with 2 elements at one place");
        }
        element->parent = nullptr;
        giveElementToSubnodes(element);
        element.reset();
    }
    giveElementToSubnodes(e);
    updateDiameter();
}


size_t Node::elementsCount() const
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

DistToNode Node::getDistsToNode(Position pos) const
{
    DistToNode result;
    if (element != nullptr)
    {
        result.nearest = result.farest = element->pos.distTo(pos);
        return result;
    }
    if (isInside(pos))
    {
        result.nearest = 0.0;
        result.farest = m_corners[0].distTo(pos);
        for (int i=1; i<8; i++)
        {
            double dist = m_corners[i].distTo(pos);
            if (result.farest < dist)
                result.farest = dist;
        }
        return result;
    }

    result.nearest = result.farest = m_corners[0].distTo(pos);
    for (int i=1; i<8; i++)
    {
        double dist = m_corners[i].distTo(pos);
        if (result.farest < dist)
            result.farest = dist;
        if (result.nearest > dist)
            result.nearest = dist;
    }

    /*
    result.nearest = -1;
    result.farest = -1;
    Position corner;

    double hs = size*0.5;

    corner.x[0] = center.x[0] + hs;
    corner.x[1] = center.x[1] + hs;
    corner.x[2] = center.x[2] + hs;
    result.nearest = result.farest = (corner - pos).len();

    for (int x = -1; x <=1; x += 2)
        for (int y = -1; y <=1; y += 2)
            for (int z = -1; z <=1; z += 2)
            {
                corner.x[0] = center.x[0] + x*hs;
                corner.x[1] = center.x[1] + y*hs;
                corner.x[2] = center.x[2] + z*hs;
                double dist = (corner - pos).len();
                if (result.farest < dist)
                    result.farest = dist;
                if (result.nearest > dist)
                    result.nearest = dist;
            }*/
    return result;
}

double Node::getMinDist(const Position& pos) const
{
    if (element != nullptr)
    {
        return element->pos.distTo(pos);
    }
    double minDist = m_corners[0].distTo(pos);
    for (int i=1; i<8; i++)
    {
        double dist = m_corners[i].distTo(pos);
        if (minDist > dist)
            minDist = dist;
    }
    return minDist;
}

double Node::getDistToCenter(const Position& pos) const
{
    return pos.distTo(center);
}

bool Node::isInside(const Position& pos) const
{
    const double *p = pos.x;
    const double *c = center.x;
    double hs = size*0.5;
    return (p[0] >= c[0] - hs) & (p[0] < c[0] + hs)
            & (p[1] >= c[1] - hs) & (p[1] < c[1] + hs)
            & (p[2] >= c[2] - hs) & (p[2] < c[2] + hs);
}

void Node::dbgOutCoords(std::ostream& s) const
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

void Node::updateMassCenter()
{
    if (element != nullptr)
    {
        massCenter = element->pos;
        mass = element->value;
        return;
    }
    massCenter = {0.0, 0.0, 0.0};
    mass = 0.0;
    for (int i=0; i<8; i++)
    {
        if (subnodes[i] != nullptr)
        {
            double nodeMass = subnodes[i]->mass;
            massCenter += subnodes[i]->massCenter * nodeMass;
            mass += nodeMass;
        }
    }
    if (mass != 0.0)
        massCenter /= mass;
    else
    {
        massCenter = center;
    }
}

void Node::updateMassCenterReqursiveUp()
{
    updateMassCenter();
    if (parent != nullptr)
        parent->updateMassCenterReqursiveUp();
}

void Node::updateMassCenterReqursiveDown()
{
    /// @todo May be optimized by using one cycle for calls and sum calculations, but this will duplicate code
    for (int i=0; i<8; i++)
        if (subnodes[i] != nullptr)
            subnodes[i]->updateMassCenterReqursiveDown();
    updateMassCenter();
}

void Node::giveElementToSubnodes(std::shared_ptr<Element> e)
{
    SubdivisionPos targerSubdivision(center, e->pos);

    int index = targerSubdivision.index();
    if (subnodes[index] == nullptr)
    {
        subnodes[index].reset(new Node(m_octree, targerSubdivision, this));
        hasSubnodes = true;
    }
    subnodes[index]->addElement(e);
}

void Node::calculateCorners()
{
    double hs = size*0.5;
    int i = 0;
    for (int x = -1; x <=1; x += 2)
        for (int y = -1; y <=1; y += 2)
            for (int z = -1; z <=1; z += 2)
            {
                m_corners[i].x[0] = center.x[0] + x*hs;
                m_corners[i].x[1] = center.x[1] + y*hs;
                m_corners[i].x[2] = center.x[2] + z*hs;
                i++;
            }
}

void Node::updateDiameter()
{
    if (element == nullptr)
        dia = size * sqrt(3.0);
    else
        dia = 0.0;
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

void Octree::clear()
{
    m_root.reset();
    m_centerIsSet = false;
}

bool Octree::empty() const
{
    return m_root == nullptr;
}

void Octree::add(std::shared_ptr<Element> e)
{
    // Creating root if no
    if (m_root == nullptr)
    {
        if (!m_centerIsSet)
        {
            m_center = e->pos;

            /**
             * We should not put grid center directly into the point due to double
             * computetion errors: it may be concerned as a point from mode than one subnodes,
             * because subnodes centers are not inaccurate.
             *
             * If you know better way to get rid of floating point errors, do it.
             */
            m_center[0] -= m_initialSize * 0.13;
            m_center[1] -= m_initialSize * 0.13;
            m_center[2] -= m_initialSize * 0.13;
            m_centerIsSet = true;
        }

        m_root.reset(
            new Node(this, m_center, m_initialSize)
        );
    }
    // Enlarging root cell
    while (!m_root->isInside(e->pos))
    {
        enlargeSpaceIteration(e->pos);
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

const Element& Octree::getNearest(Position pos)
{
    using namespace std;
    if (m_root == nullptr)
        throw(std::runtime_error("Octree is empty"));
    // ndp = node-distance pair
    using ndp = pair<const Node*, DistToNode>;
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
            const Node& n = *(it->first);
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
    return *(nodes.front().first->element);
}

void Octree::getClose(std::vector<Element*>& target, const Position& pos, double dist) const
{
    std::vector<const Node*> nodesVector;
    nodesVector.reserve(200);
    if (empty())
        return;

    nodesVector.push_back(&root());
    for (size_t i=0; i != nodesVector.size(); i++)
    {
        const Node *n = nodesVector[i];
        DistToNode nodeDist = n->getDistsToNode(pos);
        // All node is too far
        if (nodeDist.nearest > dist)
            continue;
        // All node is enough close
        if (nodeDist.farest <= dist)
        {
            n->pushBackAllElements(target);
            continue;
        }

        // Some parts are close and some are far. Need division
        n->pushBackSubnodes(nodesVector);
    }
}

const Node& Octree::root() const
{
    return *m_root;
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
    if (empty())
        return;
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
    std::unique_ptr<Node> n(new Node(this, newRootCenter, m_root->size * 2));
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

LinearScales::LinearScales(double k) :
    m_k(k)
{
}

double LinearScales::findScale(double distance) const
{
    if (distance < 0.0)
        return 0.0;
    return distance*m_k;
}

///////////////////////////
/// ScalesConfig
DiscreteScales::DiscreteScales()
{
    addScale(0.0, 0.0);
}

void DiscreteScales::addScale(double minDistance, double averagingScale)
{
    m_distsScales.push_back(std::pair<double, double>(minDistance, averagingScale));
    sortDistsScales(); // Not very quick solutions, but we should not rewrite ScalesConfig often
}

void DiscreteScales::sortDistsScales()
{
    std::sort(m_distsScales.begin(), m_distsScales.end(),
        [](const std::pair<double, double> p1, std::pair<double, double> p2)
        { return p1.first < p2.first; }
    );
}

double DiscreteScales::findScale(double distance) const
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
