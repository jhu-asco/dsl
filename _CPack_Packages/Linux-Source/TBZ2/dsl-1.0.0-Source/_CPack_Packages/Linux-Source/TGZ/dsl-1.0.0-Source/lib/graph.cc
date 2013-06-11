#include "graph.h"
#include "search.h"
#include <iostream>

using namespace dsl;
using namespace std;

Graph::Graph() : search(0)
{
}

Graph::~Graph()
{
}

void Graph::AddVertex(Vertex &v)
{
  vertices[v.id] = &v;
  
}

void Graph::RemoveVertex(Vertex &v, bool re, bool del)
{
  //  cout << "RemoveVretex: <" << endl;
  if (re) {
    std::map<int, Edge*>::iterator i;
    // remove all incoming edges
    for (i = v.in.begin(); i != v.in.end(); ++i)
      RemoveEdge(*i->second, true, del);
    // remove all outgoing edges
    for (i = v.out.begin(); i != v.out.end(); ++i)
      RemoveEdge(*i->second, true, del);
  }

  // remove from list of vertices
  vertices.erase(v.id);

  //  cout << "RemoveVretex: ." << endl;
  //  cout << v << endl;

  if (search && search->last)
    search->Remove(v);

  if (del)
    delete &v;
}


void Graph::AddEdge(Edge &e)
{
  // attach to start node
  if (e.from)
    e.from->out[e.id] = &e;
  // attach to end node
  if (e.to)
    e.to->in[e.id] = &e;
  // add to list of graph edges
  edges[e.id] = &e;
 
  // if there's an active search
  if (search && search->last) {
    double cost = e.cost;
    e.cost = INF;
    search->ChangeCost(e, cost);
  }
}


void Graph::RemoveEdge(Edge &e, bool update, bool del)
{
  // if there's an active search
  if (update && search && search->last) {
    Vertex *u = e.from;
    Vertex *v = e.to;
    
    double cost = e.cost;
    e.costChange = INF - cost;
    e.cost = INF;

    //    cout << e.costChange << " " << e.cost << " " << e.cost - e.costChange << endl;
    //    cout << "u->rhs=" << u->rhs << " v->g=" << v->g << " cost=" << cost << endl;

    // update the start vertex

    
    if (u && u->openListNode && v && v->openListNode) {
      if (search->Eq(u->rhs, cost + v->g) && u != search->goal) {
        search->MinSucc(&u->rhs, *u);        
      }
      search->UpdateVertex(*u);
    }
  }

  // remove from list of edges
  edges.erase(e.id);
  // remove from start node
  if (e.from)
    e.from->out.erase(e.id);
  // remove from end node
  if (e.to)
    e.to->in.erase(e.id);

  if (del)
    delete &e;
}

bool Graph::Exists(const Vertex &v) const
{
  return (vertices.find(v.id) != vertices.end());
}
