// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include "fibheap.h"
#include "search.h"

using namespace std;
using namespace dsl;

#define DSL_MIN(a,b) ((a<b)?(a):(b))

// these are needed by fibheap
int FIBHEAPKEY_SIZE = 2*sizeof(double);
fibheapkey_t FIBHEAPKEY_MIN = (void*)(double[2]){-INF, -INF};

extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b)
{
  assert(a); assert(b);
  double* af = (double *)a;
  double* bf = (double *)b;
  if ((af[0] < bf[0]) || (af[0] == bf[0] && af[1] < bf[1]))
    return -1;
  if (af[0] == bf[0] && af[1] == bf[1])
    return 0;
  return 1;
}


Search::Search(Graph &graph, const Cost &cost) : 
  graph(graph),
  cost(cost),
  start(0), 
  goal(0), 
  last(0), 
  km(0), 
  eps(1e-10)
{
  openList = fibheap_new();
}


Search::~Search()
{
  fibheap_delete(openList);
  changedEdges.clear();
}


void Search::Reset()
{
  std::map<int,Vertex*>::iterator vi;
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->Reset();
  }
  fibheap_clear(openList);
  km = 0;
  last = 0;
}


void Search::SetStart(const Vertex &s)
{
  start = (Vertex*)&s;
}


void Search::SetGoal(const Vertex &s)
{
  if (!start) {
    cout << "[W] Search::SetGoal: start should be set first!" << endl;
    return;
  }
  
  // reset planner
  Reset();
  // set goal
  goal = (Vertex*)&s;
  goal->rhs = 0;
  goal->key[0] = cost.Heur(*start, *goal);
  goal->key[1] = 0;
  InsertExt(*goal, goal->key);
}


void Search::ChangeCost(Edge &edge, double cost)
{
  if (Eq(cost, edge.cost))
    return;
  edge.costChange = cost - edge.cost;
  edge.cost = cost;
  changedEdges.push_back(&edge);
}


//#define STDOUT_DEBUG

void Search::UpdateVertex(Vertex &u)
{
#ifdef STDOUT_DEBUG
  printf("UpdateVertex: begin\n");
#endif
  if (u.g != u.rhs && u.t == Vertex::OPEN) {
#ifdef STDOUT_DEBUG
    printf("UpdateVertex: 1\n");
#endif
    Update(u);
  } else if (u.g != u.rhs && u.t != Vertex::OPEN) {
#ifdef STDOUT_DEBUG
    printf("UpdateVertex: 2\n");
#endif
    Insert(u);
  } else if (u.g == u.rhs && u.t == Vertex::OPEN) {
#ifdef STDOUT_DEBUG
    printf("UpdateVertex: 3\n");
#endif
    Remove(u);
  }
}


void Search::ComputeShortestPath()
{
  Vertex *u;
  Vertex *s;
  double kold[2];
  double gOld;
  map<int, Edge*>::iterator ei;
  Edge* edge;

  graph.search = this;

#ifdef STDOUT_DEBUG
  printf("ComputeShortestPath: begin\n");
#endif  
 
  //  assert(!fibheap_empty(openList));
  if (fibheap_empty(openList)) {
    cerr << "[W] Search::ComputeShortestPath: openList is empty -- most likely this means that a there is no path between start and goal!" << endl;
    return;
  }
  
  while(TopKey() && (fibkey_compare(TopKey(), CalculateKey(*start)) < 0 || start->rhs != start->g)) {
    u = Top();
#ifdef STDOUT_DEBUG
    printf("ComputeShortestPath:Top() -> ");// u->Print(stdout);
#endif

    kold[0] = u->key[0];
    kold[1] = u->key[1];
    
    if (fibkey_compare(kold, CalculateKey(*u)) < 0) {
#ifdef STDOUT_DEBUG  
      printf("ComputeShortestPath: 1 -> ");
      printf("old:[%.2f %.2f]\nnew:[%.2f %.2f]\n", kold[0],kold[1], u->key[0], u->key[1]);
#endif

      Update(*u);
    } else
    if (u->g > u->rhs) {
#ifdef STDOUT_DEBUG
      printf("ComputeShortestPath: 2\n");
#endif

      u->g = u->rhs;
      Remove(*u);

      for (ei = u->in.begin(); ei != u->in.end(); ++ei) {
        edge = (Edge*)ei->second;
        s = edge->from;
	if (s != goal)
	  s->rhs = DSL_MIN(s->rhs, edge->cost + u->g);
	UpdateVertex(*s);
      }
    } else
    {
#ifdef STDOUT_DEBUG
      printf("ComputeShortestPath: 3\n");
#endif
      
      gOld = u->g;	
      u->g = INF;
      
      for (ei = u->in.begin(); ei != u->in.end(); ++ei) {
        edge = (Edge*)ei->second;
        s = edge->from;
	if (Eq(s->rhs, edge->cost + gOld))
	  if (s != goal)
	    MinSucc(&s->rhs, *s);
	UpdateVertex(*s);
      }

      if (u != goal)
	MinSucc(&u->rhs, *u);
      UpdateVertex(*u);
    }
  }
}



int Search::Plan()
{
  Vertex* cur = start;
  Vertex *u, *v;
  int count = 1;
  
  assert(start);
  
  if (!start->out.size()) {
    cerr << "[W] Search::Plan: start vertex has no outgoing edges!" << endl;
    return 0;
  }

  if (!goal->in.size()) {
    cerr << "[W] Search::Plan: goal vertex has no incoming edges!" << endl;
    return 0;
  }

  if (!last)
    last = start;

  if (changedEdges.size()) {
    km += (cost.Heur(*last, *start));
    last = start;
    vector<Edge*>::iterator ei;
    for (ei = changedEdges.begin(); ei != changedEdges.end(); ++ei) {
      Edge* edge = *ei;
      u = edge->from;
      v = edge->to;
      
      if (edge->costChange < 0) {
        if (u != goal) 
          // new cost
          u->rhs = DSL_MIN(u->rhs, edge->cost + v->g);
      } else {
        // old cost
        if (Eq(u->rhs, edge->cost - edge->costChange + v->g)) {
          if (u != goal) {
            MinSucc(&u->rhs, *u);
          }
        }
      }
      UpdateVertex(*u); 
    }
    changedEdges.clear();
  }

  ComputeShortestPath();
  
  do {
    Vertex *next = MinSucc(0, *cur);
    cur->next = next;
    if (!next) {
      break;
    }
    next->prev = cur;
    cur = next;
    count++;
  } while(cur != goal);

  /*
  do {
    cur->next = MinSucc(0, *cur);
    cur = cur->next;
    if (!cur) {
      break;
    }
    count++;
  } while(cur != goal);
  */
  return count;
}


Vertex* Search::MinSucc(double *minRhs, const Vertex &s)
{
  double minVal = INF;
  Vertex* minSucc = NULL;
  Vertex* s_;
  map<int,Edge*>::const_iterator ei;
  Edge *edge;

  for (ei = s.out.begin(); ei != s.out.end(); ++ei) {
    edge = (Edge*)ei->second;
    s_ = edge->to;
    
    double val = edge->cost + s_->g;
#ifdef DSL_GOAL_BIAS
    if (val < minVal || 
	(val == minVal && minSucc && cost.Real(*s_, *goal) < cost.Real(*minSucc, *goal))) {
#else
   if (val < minVal) {
#endif
      minVal = val;
      minSucc = s_;
    }
  }
  if (minRhs)
    *minRhs = minVal;
  return minSucc;
}


double* Search::CalculateExtKey(double *key, Vertex &s)
{
  double m = DSL_MIN(s.g, s.rhs);
  if (m == INF) {
    key[0] = INF;
    key[1] = INF;
  } else {
    assert(start);
    key[0] = m + cost.Heur(*start, s) + km;
    key[1] = m;
  }

  return key;
}

double* Search::CalculateKey(Vertex &s)
{
  return CalculateExtKey(s.key, s);
}


void Search::Insert(Vertex &s)
{
  InsertExt(s, CalculateExtKey(s.key, s));
}

void Search::InsertExt(Vertex &s, double *key)
{

#ifdef DSL_DSTAR_MIN
  s.r = start;
#endif
  s.openListNode = fibheap_insert(openList, (void*)key, &s);
  s.t = Vertex::OPEN;
}

void Search::Update(Vertex &s)
{
  double key[2];
  assert(s.t == Vertex::OPEN);
  CalculateExtKey(key, s);

#ifdef DSL_DSTAR_MIN
  s.r = start;
#endif

  if (fibkey_compare(key, s.key) > 0) {
    fibheap_delete_node(openList, s.openListNode);
    s.openListNode = fibheap_insert(openList, key, &s); 
  } else 
  if (!fibkey_compare(key, s.key)) {
    return;
  } else {
    fibheap_replace_key(openList, s.openListNode, key);
  }
  // copy back to s's own key
  // s->openListNode->key = s->key;
  s.key[0] = key[0]; s.key[1] = key[1];
}


void Search::Remove(Vertex &s)
{
  if (s.t == Vertex::CLOSED)
    return;

  s.t = Vertex::CLOSED;
  if (s.openListNode)
    fibheap_delete_node(openList, s.openListNode);
}

Vertex* Search::Pop()
{
  printf("Pop:deprecated.\n");
  assert(0);
  /*
#ifdef INSPECT
  inspector.heapPercs++;
#endif

  if (fibheap_empty(openList))
    return NULL;
  State* s = (State*)fibheap_extract_min(openList);
  s->t = Vertex::CLOSED;
  return s;
  */
}


Vertex* Search::Top() 
{
#ifdef DSL_DSTAR_MIN
  Vertex* s;
  while((s = (Vertex*)fibheap_min(openList))) {
    if (s->r != start)
      Update(*s);
    else
      return s;
  }
  return NULL;
#else  
  return (Vertex*)fibheap_min(openList);
#endif
}


double* Search::TopKey()
{  
  Vertex* s = Top();
  if (!s)
    return NULL;
  return s->key;
}
