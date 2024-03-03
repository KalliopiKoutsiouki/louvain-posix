#pragma once
#ifndef DSPL_HPP
#define DSPL_HPP

#include <chrono>
#include <pthread.h>
#include "graph.hpp"
#include "utils.hpp"

struct Comm {
  GraphElem size;
  GraphWeight degree;

  Comm() : size(0), degree(0.0) {};
};

struct CommInfo {
    GraphElem community;
    GraphElem size;
    GraphWeight degree;
};

struct clmap_t {
  GraphElem f;
  GraphElem s;
};


struct ThreadInfo {
    GraphElem start;
    GraphElem end;
    const GraphElem *edge_indices;
    const Edge *edge_list;
    const GraphElem *currComm;
    GraphElem *targetComm;
    const GraphWeight *vDegree;
    Comm *localCinfo;
    Comm *localCupdate;
    GraphWeight constantForSecondTerm;
    GraphWeight *clusterWeight;
};


pthread_mutex_t localCinfoMutex;
pthread_mutex_t localCupdateMutex;
pthread_mutex_t clusterWeightMutex;

#define CLMAP_MAX_NUM 32
#define COUNT_MAX_NUM 32

void sumVertexDegree(const Graph &g, std::vector<GraphWeight> &vDegree, std::vector<Comm> &localCinfo)
{
  const GraphElem nv = g.get_nv();

  for (GraphElem i = 0; i < nv; i++) {
    GraphElem e0, e1;

    g.edge_range(i, e0, e1);

    for (GraphElem k = e0; k < e1; k++) {
      const Edge &edge = g.get_edge(k);
      vDegree[i] += edge.weight_;
    }

    localCinfo[i].degree = vDegree[i];
    localCinfo[i].size = 1L;
  }
} // sumVertexDegree

GraphWeight calcConstantForSecondTerm(const std::vector<GraphWeight> &vDegree)
{
  GraphWeight localWeight = 0.0;

  const size_t vsz = vDegree.size();

  for (GraphElem i = 0; i < vsz; i++) {
    localWeight += vDegree[i]; // Local reduction
  }

  return (1.0 / static_cast<GraphWeight>(localWeight));
} // calcConstantForSecondTerm

void initComm(std::vector<GraphElem> &pastComm, std::vector<GraphElem> &currComm)
{
  const size_t csz = currComm.size();

#ifdef DEBUG_PRINTF  
  assert(csz == pastComm.size());
#endif

  for (GraphElem i = 0L; i < csz; i++) {
    pastComm[i] = i;
    currComm[i] = i;
  }
} // initComm

void initLouvain(const Graph &g, std::vector<GraphElem> &pastComm, 
        std::vector<GraphElem> &currComm, std::vector<GraphWeight> &vDegree, 
        std::vector<GraphWeight> &clusterWeight, std::vector<Comm> &localCinfo, 
        std::vector<Comm> &localCupdate, GraphWeight &constantForSecondTerm)
{
  const GraphElem nv = g.get_nv();

  vDegree.resize(nv);
  pastComm.resize(nv);
  currComm.resize(nv);
  clusterWeight.resize(nv);
  localCinfo.resize(nv);
  localCupdate.resize(nv);
 
  sumVertexDegree(g, vDegree, localCinfo);
  constantForSecondTerm = calcConstantForSecondTerm(vDegree);

  initComm(pastComm, currComm);
} // initLouvain

GraphElem getMaxIndex(clmap_t *clmap, int &clmap_size, GraphWeight *counter, int &counter_size,
		      const GraphWeight selfLoop, const Comm *localCinfo, const GraphWeight vDegree, 
		      const GraphElem currSize, const GraphWeight currDegree, const GraphElem currComm,
		      const GraphWeight constant)
{
  clmap_t *storedAlready;
  GraphElem maxIndex = currComm;
  GraphWeight curGain = 0.0, maxGain = 0.0;
  GraphWeight eix = counter[0] - selfLoop;

  GraphWeight ax = currDegree - vDegree;
  GraphWeight eiy = 0.0, ay = 0.0;

  GraphElem maxSize = currSize; 
  GraphElem size = 0;

  for (storedAlready = clmap; storedAlready != clmap + clmap_size; storedAlready++) {
      if (currComm != storedAlready->f) {
          ay = localCinfo[storedAlready->f].degree;
          size = localCinfo[storedAlready->f].size;   

          if (storedAlready->s < counter_size) {
            eiy = counter[storedAlready->s];
	  }

          curGain = 2.0 * (eiy - eix) - 2.0 * vDegree * (ay - ax) * constant;

          if ((curGain > maxGain) || ((curGain == maxGain) && (curGain != 0.0) && (storedAlready->f < maxIndex))) {
              maxGain = curGain;
              maxIndex = storedAlready->f;
              maxSize = size;
          }
      }
  }

  if ((maxSize == 1) && (currSize == 1) && (maxIndex > currComm)) {
    maxIndex = currComm;
  }

  return maxIndex;
} // getMaxIndex

GraphWeight buildLocalMapCounter(const GraphElem e0, const GraphElem e1, clmap_t *clmap, int &clmap_size, 
				 GraphWeight *counter, int &counter_size, const Edge *edge_list, const GraphElem *currComm,
				 const GraphElem vertex)
{
  GraphElem numUniqueClusters = 1L;
  GraphWeight selfLoop = 0;
  clmap_t *storedAlready;
  for (GraphElem j = e0; j < e1; j++) {
        
    const Edge &edge = edge_list[j];
    const GraphElem &tail_ = edge.tail_;
    const GraphWeight &weight = edge.weight_;
    GraphElem tcomm;

    if (tail_ == vertex) {
      selfLoop += weight;
    }

    tcomm = currComm[tail_];

    storedAlready = clmap;
    for (int i = 0; i < clmap_size; i++, storedAlready++) {
      if (clmap[i].f == tcomm) {
        break;
      }
    }
    
    if (storedAlready != clmap + clmap_size && storedAlready->s < counter_size) {
      counter[storedAlready->s] += weight;
    } else {
        if (clmap_size < CLMAP_MAX_NUM) {
          clmap[clmap_size].f = tcomm;
          clmap[clmap_size].s = numUniqueClusters;
          clmap_size++;
        }
        if (counter_size < COUNT_MAX_NUM) {
          counter[counter_size] = weight;
          counter_size++;
        }
        numUniqueClusters++;
    }
  }

  return selfLoop;
} // buildLocalMapCounter



void execLouvainIteration(const GraphElem i, const GraphElem *edge_indices, const Edge *edge_list,
			  const GraphElem *currComm, GraphElem *targetComm, const GraphWeight *vDegree, Comm *localCinfo, Comm *localCupdate,
			  const GraphWeight constantForSecondTerm, GraphWeight *clusterWeight)
{
    GraphElem localTarget = -1;
    GraphElem e0, e1, selfLoop = 0;
    clmap_t clmap[CLMAP_MAX_NUM];
    int clmap_size = 0;
    GraphWeight counter[COUNT_MAX_NUM];
    int counter_size = 0;

    const GraphElem cc = currComm[i];
    GraphWeight ccDegree;
    GraphElem ccSize;

    ccDegree = localCinfo[cc].degree;
    ccSize = localCinfo[cc].size;

    e0 = edge_indices[i];
    e1 = edge_indices[i+1];

    if (e0 != e1) {
        clmap[0].f = cc;
        clmap[0].s = 0;
        clmap_size++;
        counter[0] = 0.0;
        counter_size++;

        selfLoop =  buildLocalMapCounter(e0, e1, clmap, clmap_size, counter, counter_size, edge_list, currComm, i);

        // Lock mutex to protect clusterWeight
        pthread_mutex_lock(&clusterWeightMutex);
        clusterWeight[i] += counter[0];
        // Unclock mutex
        pthread_mutex_unlock(&clusterWeightMutex);

        localTarget = getMaxIndex(clmap, clmap_size, counter, counter_size, selfLoop, localCinfo,
                        vDegree[i], ccSize, ccDegree, cc, constantForSecondTerm);
    } else {
        localTarget = cc;
    }

    if ((localTarget != cc) && (localTarget != -1)) {
        // Lock mutexes to protect localCinfo and localCupdate 
        pthread_mutex_lock(&localCinfoMutex);
        pthread_mutex_lock(&localCupdateMutex);

        localCupdate[localTarget].degree += vDegree[i];
        localCupdate[localTarget].size++;
        localCupdate[cc].degree -= vDegree[i];
        localCupdate[cc].size--;

        // Unlock mutexes
        pthread_mutex_unlock(&localCinfoMutex);
        pthread_mutex_unlock(&localCupdateMutex);
    }

#ifdef DEBUG_PRINTF  
    assert(localTarget != -1);
#endif
    targetComm[i] = localTarget;
} // execLouvainIteration

GraphWeight computeModularity(const Graph &g, Comm *localCinfo,
			      const GraphWeight *clusterWeight,
			      const GraphWeight constantForSecondTerm)
{
  const GraphElem nv = g.get_nv();
  GraphWeight le_xx = 0.0, la2_x = 0.0;

  for (GraphElem i = 0L; i < nv; i++) {
    le_xx += clusterWeight[i];
    la2_x += localCinfo[i].degree * localCinfo[i].degree; 
  } 

  GraphWeight currMod = (le_xx * constantForSecondTerm) - (la2_x * constantForSecondTerm * constantForSecondTerm);
#ifdef DEBUG_PRINTF  
  std::cout << "le_xx: " << le_xx << ", la2_x: " << la2_x << std::endl;
#endif

  return currMod;
} // computeModularity

void updateLocalCinfo(const GraphElem nv, Comm *localCinfo, const Comm *localCupdate)
{
    for (GraphElem i = 0L; i < nv; i++) {
        localCinfo[i].size += localCupdate[i].size;
        localCinfo[i].degree += localCupdate[i].degree;
    }
}

void cleanCWandCU(const GraphElem nv, GraphWeight *clusterWeight,
        Comm *localCupdate)
{
    for (GraphElem i = 0L; i < nv; i++) {
        clusterWeight[i] = 0;
        localCupdate[i].degree = 0;
        localCupdate[i].size = 0;
    }
} // distCleanCWandCU


void *execThreadMethod(void *arg) {
    ThreadInfo *info = (ThreadInfo*)(arg);
    for (GraphElem i = info->start; i < info->end; i++) {
        execLouvainIteration(i, info->edge_indices, info->edge_list, info->currComm,
                             info->targetComm, info->vDegree, info->localCinfo,
                             info->localCupdate, info->constantForSecondTerm,
                             info->clusterWeight);
    }
    pthread_exit(nullptr);
}


void execLouvainParallel(const GraphElem n, const GraphElem *edge_indices, const Edge *edge_list, const GraphElem *currComm,
                                   GraphElem *targetComm, const GraphWeight *vDegree, Comm *localCinfo, Comm *localCupdate,
                                   const GraphWeight constantForSecondTerm, GraphWeight *clusterWeight, const int numThreads) {
    pthread_t threads[numThreads];
    std::vector<ThreadInfo> threadInfo(numThreads);
    GraphElem chunkSize = n / numThreads;
    GraphElem start = 0, end;


    // Create threads
    for (int i = 0; i < numThreads; i++) {
        if (i==numThreads-1) {
          end = n;
        }
        else {
          end = start + chunkSize;
        }
        threadInfo[i].start = start;
        threadInfo[i].end = end;
        threadInfo[i].edge_indices = edge_indices;
        threadInfo[i].edge_list = edge_list;
        threadInfo[i].currComm = currComm;
        threadInfo[i].targetComm = targetComm;
        threadInfo[i].vDegree = vDegree;
        threadInfo[i].localCinfo = localCinfo;
        threadInfo[i].localCupdate = localCupdate;
        threadInfo[i].constantForSecondTerm = constantForSecondTerm;
        threadInfo[i].clusterWeight = clusterWeight;


        pthread_create(&threads[i], nullptr, execThreadMethod, (void *)&threadInfo[i]);
        
        start = end;
    }


    // Join threads
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], nullptr);
    }
}

GraphWeight louvainMethod(const Graph &g, const GraphWeight lower, const GraphWeight thresh, int &iters,int numThreads)
{
  // Times
  std::chrono::time_point<std::chrono::high_resolution_clock> t_start[6], t_end[6];
  std::chrono::duration<double> dur[6];
  std::vector<GraphElem> pastComm, currComm, targetComm;
  std::vector<GraphWeight> vDegree;
  std::vector<GraphWeight> clusterWeight;
  std::vector<Comm> localCinfo, localCupdate;
 
  const GraphElem nv = g.get_nv();

  GraphWeight constantForSecondTerm;
  GraphWeight prevMod = lower;
  GraphWeight currMod = -1.0;
  int numIters = 0;
  
  t_start[0] = std::chrono::high_resolution_clock::now();
  initLouvain(g, pastComm, currComm, vDegree, clusterWeight, localCinfo, localCupdate, constantForSecondTerm);
  t_end[0] = std::chrono::high_resolution_clock::now();
  dur[0] = std::chrono::duration_cast<std::chrono::microseconds>(t_end[0] - t_start[0]);
  targetComm.resize(nv);

#ifdef DEBUG_PRINTF  
  std::cout << "constantForSecondTerm: " << constantForSecondTerm << std::endl;
  std::cout << "Threshold: " << thresh << std::endl;
#endif

  const GraphElem *d_edge_indices = &g.edge_indices_[0];
  const Edge *d_edge_list = &g.edge_list_[0];
  GraphElem *d_currComm = &currComm[0];
  const GraphWeight *d_vDegree = &vDegree[0];
  GraphElem *d_targetComm = &targetComm[0];
  Comm *d_localCinfo = &localCinfo[0];
  Comm *d_localCupdate = &localCupdate[0];
  GraphWeight *d_clusterWeight = &clusterWeight[0];

  std::chrono::time_point<std::chrono::high_resolution_clock> t_start_all = std::chrono::high_resolution_clock::now();

  // start Louvain iteration
  while(true) {
#ifdef DEBUG_PRINTF  
    std::cout << "Starting Louvain iteration: " << numIters << std::endl;
#endif
    numIters++;
    t_start[1] = std::chrono::high_resolution_clock::now();    
    cleanCWandCU(nv, d_clusterWeight, d_localCupdate);
    t_end[1] = std::chrono::high_resolution_clock::now();
    dur[1]  += std::chrono::duration_cast<std::chrono::microseconds>(t_end[1] - t_start[1]);

    t_start[2] = std::chrono::high_resolution_clock::now();
    execLouvainParallel(nv, d_edge_indices, d_edge_list, d_currComm, d_targetComm, d_vDegree, d_localCinfo,d_localCupdate, constantForSecondTerm, d_clusterWeight,numThreads);
    t_end[2] = std::chrono::high_resolution_clock::now();
    dur[2]  += std::chrono::duration_cast<std::chrono::microseconds>(t_end[2] - t_start[2]);

    t_start[3] = std::chrono::high_resolution_clock::now();
    updateLocalCinfo(nv, d_localCinfo, d_localCupdate);
    t_end[3] = std::chrono::high_resolution_clock::now();
    dur[3]  += std::chrono::duration_cast<std::chrono::microseconds>(t_end[3] - t_start[3]);

    t_start[4] = std::chrono::high_resolution_clock::now();
    // compute modularity
    currMod = computeModularity(g, d_localCinfo, d_clusterWeight, constantForSecondTerm);
    t_end[4] = std::chrono::high_resolution_clock::now();
    dur[4]  += std::chrono::duration_cast<std::chrono::microseconds>(t_end[4] - t_start[4]);

    // exit criteria
    if (currMod - prevMod < thresh) {
        break;
    }

    prevMod = currMod;
    if (prevMod < lower) {
        prevMod = lower;
    }

    t_start[5] = std::chrono::high_resolution_clock::now();
    for (GraphElem i = 0; i < nv; i++) {
        GraphElem tmp = pastComm[i];
        pastComm[i] = d_currComm[i];
        d_currComm[i] = d_targetComm[i];
        d_targetComm[i] = tmp;
    }
    t_end[5] = std::chrono::high_resolution_clock::now();
    dur[5]  += std::chrono::duration_cast<std::chrono::microseconds>(t_end[5] - t_start[5]);
  } // end of Louvain iteration
  
  std::chrono::time_point<std::chrono::high_resolution_clock> t_end_all = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dur_all = std::chrono::duration_cast<std::chrono::microseconds>(t_end_all - t_start_all);

  std::cout << "Louvain initLouvain time: " << dur[0].count() << std::endl;
  std::cout << "Louvain cleanCWandCU time: " << dur[1].count() << std::endl;
  std::cout << "Louvain execLouvainIteration time: " << dur[2].count() << std::endl;
  std::cout << "Louvain updateLocalCinfo time: " << dur[3].count() << std::endl;
  std::cout << "Louvain computeModularity time: " << dur[4].count() << std::endl;
  std::cout << "Louvain update time (host): " << dur[5].count() << std::endl;
  std::cout << "Louvain execution time: " << dur_all.count() << std::endl;


  iters = numIters;

  vDegree.clear();
  pastComm.clear();
  currComm.clear();
  targetComm.clear();
  clusterWeight.clear();
  localCinfo.clear();
  localCupdate.clear();
  
  return prevMod;
} // louvainMethod

#endif // __DSPL
