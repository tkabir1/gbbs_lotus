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

#include "bridge.h"
#include "macros.h"
#include <bitset>

namespace gbbs {


// Edge Array Representation
template <class W>
struct edge_array {
  using weight_type = W;
  using edge = std::tuple<uintE, uintE, W>;

  // A sequence of edge tuples.
  sequence<edge> E;

  size_t n;  // num vertices.

  edge_array(sequence<edge>&& _E, size_t _n) : E(_E), n(_n) {}

  edge_array() {}

  size_t size() { return E.size(); }

  // Clears the edge array.
  sequence<edge>&& to_seq() {
    n = 0;
    return std::move(E);
  }
  size_t find_edges(uintE u, uintE v) {
    size_t m = size();
    size_t flag=0;
    parallel_for(0, m, [&](size_t i) {
      if ((std::get<0>(E[i])==u) && (std::get<1>(E[i])==v)){
        flag=1;
      }
      else if ((std::get<0>(E[i])==v) && (std::get<1>(E[i])==u)){
        flag=1;
      }
    });
    return flag;
  }
  void print() {
    size_t m = size();
    for (int i=0;i<m;i++){
      std::cout<<"edges number: "<<i<<" vertex: "<<get<0>(E[i])<<" , "<<get<1>(E[i])<<" ";
    }
  }
  bool** convert_to_array(int num_v) {
    bool **hub_array;
    hub_array = (bool**)malloc(num_v * sizeof(bool*));
    /*for (int i = 0; i < hub_count; i++){
        hub_array[i] = (bool*)malloc(hub_count* sizeof(bool));
    }*/
    for (int i = 0; i < num_v; i++){
        hub_array[i] = (bool*)malloc(i* sizeof(bool));
    }
    size_t m = size();
    parallel_for(0, m, [&](size_t i) {
      uintE first=get<0>(E[i]);
      uintE second=get<1>(E[i]);
      if (first>second){
        first=get<1>(E[i]);
        second=get<0>(E[i]);
      }
      hub_array[first][second]=1;
    
      //std::cout<<"edges number: "<<i<<" vertex: "<<get<0>(E[i])<<" , "<<get<1>(E[i])<<" ";
    });
    return hub_array;
  }

 uint32_t* convert_to_bit(uint32_t num_v,uintE* rank,int num) {
    uint32_t *hub_array;
    //std::cout<<"hub percentage: "<<num<<"\n";
    //size_t space=(num_v/8)+1;
    size_t space=((num*num)/32)+1;
    hub_array = (uint32_t*)calloc(space,sizeof(uint32_t));
    /*for (int i = 0; i < hub_count; i++){
        hub_array[i] = (bool*)malloc(hub_count* sizeof(bool));
    }*/
    size_t m = size();
    size_t a=10;
    int count_edge=0;
    //std::cout<<"The size of integer is" <<sizeof(a)<<"\n";
    for (int i=0;i<m;i++) {
      uint32_t first=get<0>(E[i]);
      uint32_t second=get<1>(E[i]);
      //std::cout<<"first: "<<first<<" and rank: "<<rank[first];
      //std::cout<<"second: "<<second<<" and rank: "<<rank[second];
      bool a=((rank[first]>=(num_v-num))&&(rank[second]>=(num_v-num)));
      //std::cout<<"a: "<<a;
      if (first>second){
        first=get<1>(E[i]);
        second=get<0>(E[i]);
      }
      if (a){
        count_edge++;
        first=num_v-1-rank[first];
        second=num_v-1-rank[second];
        //std::cout<<"after: "<<first<<" and second "<<second;
        uint32_t index1=((first*num)+second)/32;
        //std::cout<<"index: "<<unsigned(index1);
        uint32_t or1=((first*num)+second)%32;
        uint32_t next=1;
        //std::cout<<" index "<<unsigned(index1)<<" ";
        //std::cout<<"before: "<<std::bitset<sizeof(uint32_t)*8>(hub_array[index1])<<" and sec"<<(second%32)<<" with bit rep: "<<std::bitset<sizeof(uint32_t)*8>(next<<or1)<<"\n";
        hub_array[index1]=((next << or1) | hub_array[index1]);
        //std::cout<<"after: "<<std::bitset<sizeof(uint32_t)*8>(hub_array[index1])<<"\n";
        //hub_array[first][second]=1;
      
        //std::cout<<"edges number: "<<i<<" vertex: "<<get<0>(E[i])<<" , "<<get<1>(E[i])<<" ";
      }
      }
      //std::cout<<"first: "<<first<<" second: "<<second;
      //size_t index1=((first/64)*space)+(second/64);
    //std::cout<<" \nedge count "<<count_edge<<"\n";
    return hub_array;
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    size_t m = size();
    parallel_for(0, m, [&](size_t i) {
      uintE u, v;
      W w;
      std::tie(u, v, w) = E[i];
      f(u, v, w);
    });
  }
};

template <class W, class Graph>
inline edge_array<W> to_edge_array(Graph& G) {
  using edge = std::tuple<uintE, uintE, W>;

  size_t n = G.n;
  auto sizes = sequence<uintT>::uninitialized(n);
  parallel_for(0, n,
               [&](size_t i) { sizes[i] = G.get_vertex(i).out_degree(); });
  size_t m = parlay::scan_inplace(make_slice(sizes));
  assert(m == G.m);

  auto arr = sequence<edge>::uninitialized(m);
  parallel_for(0, n, [&](size_t i) {
    size_t idx = 0;
    uintT offset = sizes[i];
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      arr[offset + idx] = std::make_tuple(u, v, wgh);
      idx++;
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  });
  return edge_array<W>(std::move(arr), n);
}

}  // namespace gbbs
