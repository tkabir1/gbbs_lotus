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

#include "SpanningForest.h"
#include "ligra/ligra.h"
#include "benchmarks/SpanningForest/common.h"
#include "benchmarks/SpanningForest/check.h"

template <class Graph>
double SF_runner(Graph& G, commandLine P) {
  std::cout << "### Application: SpanningForest (LabelPropagation-based)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;

  timer t;
  t.start();
  pbbs::sequence<edge> edges;
  if (P.getOptionValue("-permute")) {
    edges = labelprop_sf::SpanningForest</*use_permutation=*/true>(G);
  } else {
    edges = labelprop_sf::SpanningForest</*use_permutation=*/false>(G);
  }
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  if (P.getOptionValue("-check")) {
    spanning_forest::sf_compute_and_check(G, edges);
  }

  return tt;
}

generate_main(SF_runner, false);
