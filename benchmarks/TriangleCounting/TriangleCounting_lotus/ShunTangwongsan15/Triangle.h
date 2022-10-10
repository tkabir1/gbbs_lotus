// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#include <algorithm>

#include "gbbs/gbbs.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
using namespace std;
#define HUB_PERCENTAGE 10
namespace gbbs {
template <class Graph>
struct Lotus{
  Graph &hub;//csr
  Graph &non_hub;
  bool **hub_array;///CSR format
  //sequence<sequence<uintE>> hubs;
  //sequence<sequence<uintE>> non_hubs;
  //vector <parlay::sequence<typename int, typename Allocator> non_hubs[N];
  uintE hub_count;
};
template <class Graph>
symmetric_graph<symmetric_vertex, gbbs::empty> construct_lotus(Graph& G) {

  return symmetric_graph<symmetric_vertex, gbbs::empty>(
      G);
}
template <class Graph>
struct countF {
  Graph& G;
  size_t* counts;
  countF(Graph& G, size_t* _counts) : G(G), counts(_counts) {}

  inline bool update(uintE s, uintE d) {
    auto d_neighbors = G.get_vertex(d).out_neighbors();
    gbbs::write_add(&counts[s], G.get_vertex(s).out_neighbors().intersect(
                                    &d_neighbors, s, d));
    return 1;
  }

  inline bool updateAtomic(uintE s, uintE d) {
    auto d_neighbors = G.get_vertex(d).out_neighbors();
    gbbs::write_add(&counts[s], G.get_vertex(s).out_neighbors().intersect(
                                    &d_neighbors, s, d));
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = gbbs::new_array_no_init<uintE>(n);
  sequence<uintE> o = sequence<uintE>::uninitialized(n);

  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { o[i] = i; });
  parlay::sample_sort_inplace(make_slice(o), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).out_degree() < G.get_vertex(v).out_degree();
  });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { r[o[i]] = i; });
  for (int i=0;i<n;i++)
  {
    //cout<<"i: "<<r[i]<<" ";
  }
  return r;
}

// Directly call edgemap dense-forward.
template <class Graph, class VS, class F>
inline vertexSubset emdf(Graph& G, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<gbbs::empty>(G, vs, f, fl);
}

template <class Graph>
inline size_t CountDirected(Graph& DG, size_t* counts, vertexSubset& Frontier) {
  using W = typename Graph::weight_type;
  emdf(DG, Frontier, wrap_em_f<W>(countF<Graph>(DG, counts)), no_output);
  auto count_seq =
      parlay::delayed_seq<size_t>(DG.n, [&](size_t i) { return counts[i]; });
  size_t count = parlay::reduce(count_seq);
  return count;
}

// Returns the number of directed triangles in the input graph of the following
// orientation:
//        w
//       ^ ^
//      /   \.
//     u --> v
//
// Arguments:
//   DG
//     Graph on which we'll count triangles.
//   f: (uintE, uintE, uintE) -> void
//     Function that's run each triangle. On a directed triangle like the one
//     pictured above, we run `f(u, v, w)`.
template <class Graph, class F>
inline size_t CountDirectedBalanced_lotus(Graph& DG1, Graph& DG2, size_t* counts, const F& f) {
  using W = typename Graph::weight_type;
  debug(std::cout << "Starting counting"
                  << "\n";);
  size_t n = DG1.n;
  //cout<<"in lotus value of n"<<n<<"\n";

  auto parallel_work = sequence<size_t>::uninitialized(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      //cout<<"inside function v"<<v<<" and "<<DG2.get_vertex(v).out_degree()<<"\n";
      return DG1.get_vertex(v).out_degree();//first time
    };
    parallel_for(0, n, [&](size_t i) {
      auto monoid = parlay::addm<size_t>();
      parallel_work[i] = DG1.get_vertex(i).out_neighbors().reduce(map_f, monoid);//second time
      //cout<<i<<"in seconf inside funstion: "<<parallel_work[i]<<"\n";
    });
  }
  size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

  size_t block_size = 50000;
  size_t n_blocks = total_work / block_size + 1;
  size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    //cout<<"code here\n";
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      //cout<<"value of "<<i<<" ";
      auto our_neighbors = DG2.get_vertex(i).out_neighbors();//third time
      //cout<<"neighbo\n";
      //cout<<"neighboors\n"<<our_neighbors<<"\n";
      //cout<<"neighbos"<<DG.get_vertex(i).out_neighbors()<<"\n";
      size_t total_ct = 0;
      auto map_f = [&](uintE u, uintE v, W wgh) {
        auto neighbors2=DG1.get_vertex(u).out_neighbors();
        auto their_neighbors = DG1.get_vertex(v).out_neighbors();//forth time
        //cout<<"total before"<<total_ct<<"\n";
        total_ct += neighbors2.intersect_f_par(&their_neighbors, f);
        //cout<<"total after"<<total_ct<<"\n";
      };
      our_neighbors.map(map_f, false);  // run map sequentially
      counts[i] = total_ct;
      //cout<<"counts values"<<counts[i]<<"\n";
    }
  };

  parallel_for(0, n_blocks, 1, [&](size_t i) {
    //cout<<"entered in parallel for i\n"<<i<<" and work"<<work_per_block<<"\n";
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
    size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
    //cout<<"start"<<start_ind<<" and end "<<end_ind<<"\n";
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = gbbs::make_slice<size_t>(counts, DG1.n);//fifth time
  size_t count = parlay::reduce(count_seq);
  //cout<<"count"<<count<<"\n";
  return count;
}


template <class Graph, class F>
inline size_t hub_modified_CountDirectedBalanced_lotus(Graph& DG1, Graph& DG2, size_t* counts, const F& f) {
  using W = typename Graph::weight_type;
  debug(std::cout << "Starting counting"
                  << "\n";);
  size_t n = DG1.n;
  //cout<<"in lotus value of n"<<n<<"\n";

  auto parallel_work = sequence<size_t>::uninitialized(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      //cout<<"inside function v"<<v<<" and "<<DG1.get_vertex(v).out_degree()<<"\n";
      return DG2.get_vertex(v).out_degree();//first time
    };
    parallel_for(0, n, [&](size_t i) {
      auto monoid = parlay::addm<size_t>();
      parallel_work[i] = DG2.get_vertex(i).out_neighbors().reduce(map_f, monoid);//second time
      //cout<<i<<"in seconf inside funstion: "<<parallel_work[i]<<"\n";
    });
  }
  size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

  size_t block_size = 50000;
  size_t n_blocks = total_work / block_size + 1;
  size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    //cout<<"code here\n";
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      //cout<<"value of "<<i<<" ";
      auto our_neighbors = DG2.get_vertex(i).out_neighbors();//third time
      //auto our_neighbors2=DG1.get_vertex(i).out_neighbors();//third time
      //cout<<"neighbo\n";
      //cout<<"neighboors\n"<<((uintE*)(our_neighbors->neighbors))<<"\n";
      //cout<<"neighbos"<<DG.get_vertex(i).out_neighbors()<<"\n";
      size_t total_ct = 0;

      //cout<<"original vertex"<<DG1.get_vertex(i)<<"\n";
      auto map_f = [&](uintE u, uintE v, W wgh) {
        //cout<<""
        auto their_neighbors = DG1.get_vertex(v).out_neighbors();//forth time
        //cout<<"total before"<<total_ct<<"\n";
        //cout<<"neighbor"<<&their_neighbors<<"n";
        //uintT nA = our_neighbors;
        //cout<<"\nnA"<<nA<<"Done\n";

        total_ct += our_neighbors.intersect_f_par(&their_neighbors, f);
        //cout<<"total after"<<total_ct<<"\n";
      };
      our_neighbors.map(map_f, false);  // run map sequentially
      counts[i] = total_ct;
      //cout<<"counts values"<<counts[i]<<"\n";
    }
  };

  parallel_for(0, n_blocks, 1, [&](size_t i) {
    //cout<<"entered in parallel for i\n"<<i<<" and work"<<work_per_block<<"\n";
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
    size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
    //cout<<"start"<<start_ind<<" and end "<<end_ind<<"\n";
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = gbbs::make_slice<size_t>(counts, DG1.n);//fifth time
  size_t count = parlay::reduce(count_seq);
  //cout<<"count"<<count<<"\n";
  return count;
}


template <class Graph, class F>
inline size_t hub_CountDirectedBalanced_lotus(Graph& DG1, Graph& DG2,size_t* counts, const F& f) {
  using W = typename Graph::weight_type;
  cout<<"edges number: "<<DG1.m;
  auto hub2hub_array=to_edge_array<gbbs::empty>(DG1);
  cout<<"vertex: "<<hub2hub_array.size()<<"\n";
  debug(std::cout << "Starting counting"
                  << "\n";);
  size_t n = DG1.n;
  //cout<<"in lotus value of n"<<n<<"\n";

  auto parallel_work = sequence<size_t>::uninitialized(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      //cout<<"inside function v"<<v<<" and "<<DG1.get_vertex(v).out_degree()<<"\n";
      return DG2.get_vertex(v).out_degree();//first time
    };
    parallel_for(0, n, [&](size_t i) {
      auto monoid = parlay::addm<size_t>();
      parallel_work[i] = DG2.get_vertex(i).out_neighbors().reduce(map_f, monoid);//second time
      //cout<<i<<"in seconf inside funstion: "<<parallel_work[i]<<"\n";
    });
  }
  size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

  size_t block_size = 50000;
  size_t n_blocks = total_work / block_size + 1;
  size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  //hub2hub_array.print();
  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    //cout<<"code here\n";
    //cout<<"start: "<<start_ind<<" and end: "<<end_ind<<"\n";
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      //cout<<"value of "<<i<<" ";
      auto our_neighbors = DG2.get_vertex(i).out_neighbors();//third time
      size_t total_ct = 0;
      //auto our_neighbors2=DG1.get_vertex(i).out_neighbors();//third time
      //cout<<"neighbo\n";
      //cout<<"neighboors\n"<<((uintE*)(our_neighbors->neighbors))<<"\n";
      //cout<<"neighbos"<<DG.get_vertex(i).out_neighbors()<<"\n";
      //cout<<"\n";
      uintT nA = our_neighbors.degree;
      //cout<<"original vertex"<<DG1.get_vertex(i)<<"\n";
      uintE* nghA = (uintE*)(our_neighbors.neighbors);
      //cout<<"\nneighbors of: "<<i<<" : ";
      /*for (int i=0;i<nA;i++)
      {
        std::cout<<nghA[i]<<" , ";
        //cout<<"check: "<<hub2hub_array.find_edges(nghA[i],nghA[i+1])<<" next";
        //total_ct=total_ct+hub2hub_array.find_edges(nghA[i],nghA[i+1]);
        //nghA[i]->neighbors
      }*/
      if (nA>=2){
        //cout<<" Done ";
        for (int i=0;i<nA;i++)
        {
          //std::cout<<nghA[i]<<" , "<<nghA[i+1];
          uintE u=nghA[i];
          for (int j=(i+1);j<nA;j++){
            //std::cout<<nghA[i]<<" , "<<nghA[j];
            uintE v=nghA[j];
            //cout<<"check: "<<hub2hub_array.find_edges(nghA[i],nghA[j]);
            total_ct=total_ct+hub2hub_array.find_edges(nghA[i],nghA[j]);
            total_ct=total_ct+hub2hub_array.find_edges(nghA[j],nghA[i]);
          }
          //nghA[i]->neighbors
        }
      }
      //cout<<"\n";
      auto map_f = [&](uintE u, uintE v, W wgh) {
        //cout<<""
        //cout<<"vertex: "<<hub2hub_array.find_edges(u,v)<<"\n";
        //size_t check=hub2hub_array.find_edges(u,v);
        //cout<<"u: "<<u<<"v: "<<v<<" ";
        //total_ct += our_neighbors.check_edge_lotus(&their_neighbors,f);
        //cout<<"total after"<<total_ct<<"\n";
      };
      our_neighbors.map(map_f, false);  // run map sequentially
      counts[i] = total_ct;
      //cout<<"counts values"<<counts[i]<<"\n";
    }
    //cout<<"total: "<<total_ct<<"\n";
  };

  parallel_for(0, n_blocks, 1, [&](size_t i) {
    //cout<<"entered in parallel for i\n"<<i<<" and work"<<work_per_block<<"\n";
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
    size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
    //cout<<"start"<<start_ind<<" and end "<<end_ind<<"\n";
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = gbbs::make_slice<size_t>(counts, DG1.n);//fifth time
  size_t count = parlay::reduce(count_seq);
  //cout<<"count"<<count<<"\n";
  return count;
}

// Counts the number of triangles in the input graph.
//
// Implementation note: this converts the input graph to a directed graph in
// which we point edges from lower-degree vertices to higher-degree vertices,
// hence the function name.
//
// Arguments:
//   G
//     Graph on which we'll count triangles.
//   f: (uintE, uintE, uintE) -> void
//     Function that's run each triangle. On a triangle with vertices {u, v, w},
//     we run `f(u, v, w)`.
//
// Returns:
//   The number of triangles in `G`.
template <class Graph, class F>
inline size_t Triangle_degree_ordering(Graph& G, const F& f) {
  //cout<<"in degree function\n";
  //cout<<"Graph"<<G<<"\n";
  using W = typename Graph::weight_type;
  //cout<<"W"<<W<<"\n";
  timer gt;
  gt.start();
  uintT n = G.n;
  //cout<<"vertex number: "<<n<<"\n";
  auto counts = sequence<size_t>::uninitialized(n);
  /*for (int i=0;i<counts.size();i++)
    cout<<"couts"<<counts.at(i)<<" ";
  cout<<"\n";*/
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });

  // 1. Rank vertices based on degree
  timer rt;
  rt.start();
  uintE* rank = rankNodes(G, G.n);
  rt.stop();
  rt.next("rank time");
  //cout<<"rank"<<rank<<"\n";

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto hub_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    int num=(int)((HUB_PERCENTAGE*1.0/100.0)*G.n);
    //cout<<"num: "<<num<<"\n";
    //cout<<" u "<<rank[u]<<" v "<<rank[v]<<" ";
    bool a=((rank[u]>=(G.n-num))||(rank[v]>=(G.n-num)))&&(rank[u] < rank[v]);
    //cout<<"return"<<a<<" ";
    return a;
  };
  auto hub2hub_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    int num=(int)((HUB_PERCENTAGE*1.0/100.0)*G.n);
    //cout<<"num: "<<num<<"\n";
    //cout<<" u "<<rank[u]<<" v "<<rank[v]<<" ";
    bool a=((rank[u]>=(G.n-num))&&(rank[v]>=(G.n-num)))&&(rank[u] < rank[v]);
    //cout<<"return"<<a<<" ";
    return a;
  };

  auto nonhub_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    int num=(int)((HUB_PERCENTAGE*1.0/100.0)*G.n);
    ////cout<<"num: "<<num<<"\n";
    //cout<<" u "<<rank[u]<<" v "<<rank[v]<<" ";
    bool a=((rank[u]<(G.n-num))&&(rank[v]<(G.n-num)))&&(rank[u] < rank[v]);
    //cout<<"return"<<a<<" ";
    return a;
  };
  //auto DG=filterGraph(G, pack_predicate);
  auto hub2hub=filterGraph(G, hub2hub_predicate);
  /*for (int i=0;i<hub2hub.n;i++)
  {
    auto a=hub2hub.get_vertex(i);
    cout<<"hell"<<a<<" , ";
  }*/
  auto hub = filterGraph(G, hub_predicate);
  auto non_hub = filterGraph(G, nonhub_predicate);
  //auto hub_array=create_bit_array(G,hub2hub_predicate);

  //cout<<"type of array: "<<hub2hub_array.size()<<"\n";
  //auto non_hub=filterGraph(G, pack_predicate);
  //cout<<"new n in hub graph:"<<hub.m<<" and in non hub"<<non_hub.m<<" and prev"<<DG.m<<"\n";
  //auto DG = Graph::filterGraph(G, pack_predicate);
  gt.stop();
  gt.next("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();
  //size_t count=0;
  size_t count_hub = hub_CountDirectedBalanced_lotus(hub2hub, hub,counts.begin(), f);
  size_t count_hub_mod = hub_modified_CountDirectedBalanced_lotus(hub2hub, hub,counts.begin(), f);
  //size_t count_hub = CountDirectedBalanced_lotus(hub, counts.begin(), f);
  size_t count_nonhub1 = CountDirectedBalanced_lotus(hub,non_hub, counts.begin(), f);
  size_t count_nonhub2 = CountDirectedBalanced_lotus(non_hub,non_hub, counts.begin(), f);
  cout<<"hub"<<count_hub<<"\n";
  cout<<"hub modified"<<count_hub_mod<<"\n";
  cout<<"nonhub1: "<<count_nonhub1<<"\n";
  cout<<"nonhub2: "<<count_nonhub2<<"\n";
  size_t count=count_hub+count_nonhub1+count_nonhub2;
  //size_t count=count_hub;
  std::cout << "### Num triangles = " << count << "\n";
  ct.stop();
  ct.next("count time");
  gbbs::free_array(rank, G.n);
  return count;
}


template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, const std::string& ordering,
                       commandLine& P) {
  if (ordering == "degree") {
    //cout<<"In ordering\n";
    return Triangle_degree_ordering<Graph, F>(G, f);
  } else {
    std::cerr << "Unexpected ordering: " << ordering << '\n';
    exit(1);
  }
}

}  // namespace gbbs
