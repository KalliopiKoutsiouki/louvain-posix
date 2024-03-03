#pragma once
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <climits>

#include <unistd.h>

#include "utils.hpp"

class Graph
{
    public:
        Graph(): edge_indices_(nullptr), edge_list_(nullptr), nv_(-1), ne_(-1)
        {}
                
        Graph(GraphElem nv): edge_list_(nullptr), nv_(-1), ne_(-1)
        {
            edge_indices_   = new GraphElem[nv_+1];
        }

        Graph(GraphElem nv, GraphElem ne): 
            nv_(nv), ne_(ne) 
        {
            edge_indices_   = new GraphElem[nv_+1];
            edge_list_      = new Edge[ne_];
        }

        ~Graph() 
        {

            delete [] edge_indices_;
            delete [] edge_list_;
        }
       
        Graph(const Graph &other) = delete;
        Graph& operator=(const Graph& d) = delete;
 
        void set_edge_index(GraphElem const vertex, GraphElem const e0)
        {
            edge_indices_[vertex] = e0;
        } 
        
        void edge_range(GraphElem const vertex, GraphElem& e0, GraphElem& e1) const
        {
            e0 = edge_indices_[vertex];
            e1 = edge_indices_[vertex+1];
        } 

        void set_nedges(GraphElem ne) 
        { 
            ne_ = ne;
            edge_list_      = new Edge[ne_];
        }

        GraphElem get_nv() const { return nv_; }
        GraphElem get_ne() const { return ne_; }
       
        // return edge and active info
        // ----------------------------
       
        Edge const& get_edge(GraphElem const index) const
        { return edge_list_[index]; }
         
        Edge& set_edge(GraphElem const index)
        { return edge_list_[index]; }       
                
        // print edge list (with weights)
        void print(bool print_weight = true) const
        {
            if (ne_ < MAX_PRINT_NEDGE)
            {
                for (GraphElem i = 0; i < nv_; i++)
                {
                    GraphElem e0, e1;
                    edge_range(i, e0, e1);
                    if (print_weight) { // print weights (default)
                        for (GraphElem e = e0; e < e1; e++)
                        {
                            Edge const& edge = get_edge(e);
                            std::cout << i << " " << edge.tail_ << " " << edge.weight_ << std::endl;
                        }
                    }
                    else { // don't print weights
                        for (GraphElem e = e0; e < e1; e++)
                        {
                            Edge const& edge = get_edge(e);
                            std::cout << i << " " << edge.tail_ << std::endl;
                        }
                    }
                }
            }
            else
            {
                std::cout << "Graph size is {" << nv_ << ", " << ne_ << "}, which will overwhelm STDOUT." << std::endl;
            }
        }

        void print_preview() const
        {
            std::cout << "Printing vertex#0 and all associated edges." << std::endl;
            GraphElem e0, e1;
            for (GraphElem i = 0; i < nv_; i++)
            {
                edge_range(i, e0, e1);
                if ((e1 - e0) > 0)
                {
                    for (GraphElem e = e0; e < e1; e++)
                    {
                        Edge const& edge = get_edge(e);
                        std::cout << 0 << " " << edge.tail_ << " " << edge.weight_ << std::endl;
                    }
                    break;
                }
            }
        }
       
        // print statistics about edge distribution
        void print_stats()
        {
          std::vector<GraphElem> pdeg(nv_, 0);
          for (GraphElem v = 0; v < nv_; v++)
          {
            GraphElem e0, e1;
            edge_range(v, e0, e1);
            for (GraphElem e = e0; e < e1; e++)
              pdeg[v] += 1;
          }

          std::sort(pdeg.begin(), pdeg.end());
          GraphWeight loc = (GraphWeight)(nv_ + 1)/2.0;
          GraphElem median;
          if (fmod(loc, 1) != 0)
            median = pdeg[(GraphElem)loc]; 
          else
            median = (pdeg[(GraphElem)floor(loc)] + pdeg[((GraphElem)floor(loc)+1)]) / 2;
          GraphElem spdeg = std::accumulate(pdeg.begin(), pdeg.end(), 0);
          GraphElem mpdeg = *(std::max_element(pdeg.begin(), pdeg.end()));
          std::transform(pdeg.cbegin(), pdeg.cend(), pdeg.cbegin(),
              pdeg.begin(), std::multiplies<GraphElem>{});

          GraphElem psum_sq = std::accumulate(pdeg.begin(), pdeg.end(), 0);

          GraphWeight paverage = (GraphWeight) spdeg / nv_;
          GraphWeight pavg_sq  = (GraphWeight) psum_sq / nv_;
          GraphWeight pvar     = std::abs(pavg_sq - (paverage*paverage));
          GraphWeight pstddev  = sqrt(pvar);

          std::cout << std::endl;
          std::cout << "--------------------------------------" << std::endl;
          std::cout << "Graph characteristics" << std::endl;
          std::cout << "--------------------------------------" << std::endl;
          std::cout << "Number of vertices: " << nv_ << std::endl;
          std::cout << "Number of edges: " << ne_ << std::endl;
          std::cout << "Maximum number of edges: " << mpdeg << std::endl;
          std::cout << "Median number of edges: " << median << std::endl;
          std::cout << "Expected value of X^2: " << pavg_sq << std::endl;
          std::cout << "Variance: " << pvar << std::endl;
          std::cout << "Standard deviation: " << pstddev << std::endl;
          std::cout << "--------------------------------------" << std::endl;
        }
                
        GraphElem *edge_indices_;
        Edge *edge_list_;
        
        GraphElem get_num_vertices() { return nv_;};
        GraphElem get_num_edges() {return ne_;};
        GraphElem* get_index_ranges() {return edge_indices_;};
        void* get_edge_list() {return edge_list_;};

    private:
        GraphElem nv_, ne_;
};

// read in binary edge list files using POSIX I/O
class BinaryEdgeList
{
    public:
        BinaryEdgeList() : 
            M_(-1), N_(-1)
        {}
        
        // read a file and return a graph
        Graph* read(std::string binfile, bool isUnitEdgeWeight)
        {
            std::ifstream file;

            file.open(binfile.c_str(), std::ios::in | std::ios::binary); 

            if (!file.is_open()) 
            {
                std::cout << " Error opening file! " << std::endl;
                std::abort();
            }

            // read the dimensions 
            file.read(reinterpret_cast<char*>(&M_), sizeof(GraphElem));
            file.read(reinterpret_cast<char*>(&N_), sizeof(GraphElem));

            // create local graph
            Graph *g = new Graph(M_, N_);

            uint64_t tot_bytes=(M_+1)*sizeof(GraphElem);
            ptrdiff_t offset = 2*sizeof(GraphElem);

            if (tot_bytes < INT_MAX)
                file.read(reinterpret_cast<char*>(&g->edge_indices_[0]), tot_bytes);
            else 
            {
                int chunk_bytes=INT_MAX;
                uint8_t *curr_pointer = (uint8_t*) &g->edge_indices_[0];
                uint64_t transf_bytes = 0;

                while (transf_bytes < tot_bytes)
                {
                    file.read(reinterpret_cast<char*>(&curr_pointer[offset]), chunk_bytes);
                    transf_bytes += chunk_bytes;
                    offset += chunk_bytes;
                    curr_pointer += chunk_bytes;

                    if ((tot_bytes - transf_bytes) < INT_MAX)
                        chunk_bytes = tot_bytes - transf_bytes;
                } 
            }    

            N_ = g->edge_indices_[M_] - g->edge_indices_[0];
            g->set_nedges(N_);
            tot_bytes = N_*(sizeof(Edge));
            offset = 2*sizeof(GraphElem) + (M_+1)*sizeof(GraphElem) + g->edge_indices_[0]*(sizeof(Edge));

            if (tot_bytes < INT_MAX)
                file.read(reinterpret_cast<char*>(&g->edge_list_[0]), tot_bytes);
            else 
            {
                int chunk_bytes=INT_MAX;
                uint8_t *curr_pointer = (uint8_t*)&g->edge_list_[0];
                uint64_t transf_bytes = 0;

                while (transf_bytes < tot_bytes)
                {
                    file.read(reinterpret_cast<char*>(&curr_pointer[offset]), tot_bytes);
                    transf_bytes += chunk_bytes;
                    offset += chunk_bytes;
                    curr_pointer += chunk_bytes;

                    if ((tot_bytes - transf_bytes) < INT_MAX)
                        chunk_bytes = (tot_bytes - transf_bytes);
                } 
            }   

            file.close();

            for(GraphElem i=1;  i < M_+1; i++)
                g->edge_indices_[i] -= g->edge_indices_[0];   
            g->edge_indices_[0] = 0;
            
            return g;
        }
    private:
        GraphElem M_, N_;
};
#endif
