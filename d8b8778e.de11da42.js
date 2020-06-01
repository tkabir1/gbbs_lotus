(window.webpackJsonp=window.webpackJsonp||[]).push([[32],{134:function(e,a,t){"use strict";t.r(a),t.d(a,"frontMatter",(function(){return m})),t.d(a,"metadata",(function(){return p})),t.d(a,"rightToc",(function(){return b})),t.d(a,"default",(function(){return i}));var n=t(1),s=t(6),c=(t(0),t(145)),m={id:"positive_weight_sssp",title:"Positive-Weight SSSP (Delta Stepping)"},p={id:"benchmarks/sssp/positive_weight_sssp",title:"Positive-Weight SSSP (Delta Stepping)",description:"## Problem Specification",source:"@site/docs/benchmarks/sssp/positive_weight_sssp.md",permalink:"/gbbs/docs/benchmarks/sssp/positive_weight_sssp",editUrl:"https://github.com/facebook/docusaurus/edit/master/website/docs/benchmarks/sssp/positive_weight_sssp.md",sidebar:"docs",previous:{title:"Integral-Weight SSSP (weighted BFS)",permalink:"/gbbs/docs/benchmarks/sssp/integral_weight_sssp"},next:{title:"General-Weight SSSP (Bellman-Ford)",permalink:"/gbbs/docs/benchmarks/sssp/general_weight_sssp"}},b=[{value:"Problem Specification",id:"problem-specification",children:[]},{value:"Algorithm Implementations",id:"algorithm-implementations",children:[]},{value:"Cost Bounds",id:"cost-bounds",children:[]},{value:"Compiling and Running",id:"compiling-and-running",children:[]},{value:"References",id:"references",children:[]}],r={rightToc:b};function i(e){var a=e.components,t=Object(s.a)(e,["components"]);return Object(c.b)("wrapper",Object(n.a)({},r,t,{components:a,mdxType:"MDXLayout"}),Object(c.b)("h2",{id:"problem-specification"},"Problem Specification"),Object(c.b)("h4",{id:"input"},"Input"),Object(c.b)("p",null,Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"G"),Object(c.b)("mo",{parentName:"mrow"},"="),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{stretchy:"false"}),"("),Object(c.b)("mi",{parentName:"mrow"},"V"),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{separator:"true"}),","),Object(c.b)("mi",{parentName:"mrow"},"E"),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{separator:"true"}),","),Object(c.b)("mi",{parentName:"mrow"},"w"),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{stretchy:"false"}),")")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"G=(V, E, w)")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault"}),"G"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mrel"}),"="),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.2777777777777778em"}}))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mopen"}),"("),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.22222em"}}),"V"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mpunct"}),","),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.16666666666666666em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.05764em"}}),"E"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mpunct"}),","),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.16666666666666666em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.02691em"}}),"w"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mclose"}),")"))))),", a weighted graph with positive edge weights, and a\nsource, ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"s"),Object(c.b)("mo",{parentName:"mrow"},"\u2208"),Object(c.b)("mi",{parentName:"mrow"},"V")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"s \\in V")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.5782em",verticalAlign:"-0.0391em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault"}),"s"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mrel"}),"\u2208"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mspace",style:{marginRight:"0.2777777777777778em"}}))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.22222em"}}),"V"))))),". The input graph can either be undirected, or\ndirected."),Object(c.b)("h4",{id:"output"},"Output"),Object(c.b)("p",null,"Output: ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"D")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"D")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.02778em"}}),"D"))))),", a mapping where ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"D"),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{stretchy:"false"}),"["),Object(c.b)("mi",{parentName:"mrow"},"v"),Object(c.b)("mo",Object(n.a)({parentName:"mrow"},{stretchy:"false"}),"]")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"D[v]")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.02778em"}}),"D"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mopen"}),"["),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.03588em"}}),"v"),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mclose"}),"]")))))," is the shortest path distance from\n",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"s")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"s")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault"}),"s")))))," to ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"v")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"v")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.03588em"}}),"v")))))," in ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"G")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"G")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault"}),"G")))))," and ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",Object(n.a)({parentName:"mrow"},{mathvariant:"normal"}),"\u221e")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"\\infty")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord"}),"\u221e")))))," if ",Object(c.b)("span",Object(n.a)({parentName:"p"},{className:"math math-inline"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-mathml"}),Object(c.b)("math",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/1998/Math/MathML"}),Object(c.b)("semantics",{parentName:"math"},Object(c.b)("mrow",{parentName:"semantics"},Object(c.b)("mi",{parentName:"mrow"},"v")),Object(c.b)("annotation",Object(n.a)({parentName:"semantics"},{encoding:"application/x-tex"}),"v")))),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"katex-html","aria-hidden":"true"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"base"}),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}})),Object(c.b)("span",Object(n.a)({parentName:"span"},{className:"mord mathdefault",style:{marginRight:"0.03588em"}}),"v")))))," is unreachable."),Object(c.b)("h2",{id:"algorithm-implementations"},"Algorithm Implementations"),Object(c.b)("p",null,"The code for our implemenation is available\n",Object(c.b)("a",Object(n.a)({parentName:"p"},{href:"https://github.com/ldhulipala/gbbs/tree/master/benchmarks/PositiveWeightSSSP/DeltaStepping"}),"here"),".\nThe implementation is from the Julienne paper ","[1]","."),Object(c.b)("h2",{id:"cost-bounds"},"Cost Bounds"),Object(c.b)("p",null,"Please ","[1]"," for details."),Object(c.b)("h2",{id:"compiling-and-running"},"Compiling and Running"),Object(c.b)("p",null,"The benchmark can be compiled by running:"),Object(c.b)("pre",null,Object(c.b)("code",Object(n.a)({parentName:"pre"},{}),"bazel build -c opt //benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping\n")),Object(c.b)("p",null,"It can then be run on a test input graph in the ",Object(c.b)("em",{parentName:"p"},"uncompressed format")," as follows:"),Object(c.b)("pre",null,Object(c.b)("code",Object(n.a)({parentName:"pre"},{}),"numactl -i all ./bazel-bin/benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping_main -s -m -src 1 inputs/rMatGraph_J_5_100\n")),Object(c.b)("p",null,"It can then be run on a test input graph in the ",Object(c.b)("em",{parentName:"p"},"compressed format")," as follows:"),Object(c.b)("pre",null,Object(c.b)("code",Object(n.a)({parentName:"pre"},{}),"numactl -i all ./bazel-bin/benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda\n")),Object(c.b)("h2",{id:"references"},"References"),Object(c.b)("p",null,"[1]"," Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. ",Object(c.b)("a",Object(n.a)({parentName:"p"},{href:"https://ldhulipala.github.io/papers/Bucketing.pdf"}),"Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing"),". Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017."))}i.isMDXComponent=!0},145:function(e,a,t){"use strict";t.d(a,"a",(function(){return l})),t.d(a,"b",(function(){return o}));var n=t(0),s=t.n(n);function c(e,a,t){return a in e?Object.defineProperty(e,a,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[a]=t,e}function m(e,a){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);a&&(n=n.filter((function(a){return Object.getOwnPropertyDescriptor(e,a).enumerable}))),t.push.apply(t,n)}return t}function p(e){for(var a=1;a<arguments.length;a++){var t=null!=arguments[a]?arguments[a]:{};a%2?m(Object(t),!0).forEach((function(a){c(e,a,t[a])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):m(Object(t)).forEach((function(a){Object.defineProperty(e,a,Object.getOwnPropertyDescriptor(t,a))}))}return e}function b(e,a){if(null==e)return{};var t,n,s=function(e,a){if(null==e)return{};var t,n,s={},c=Object.keys(e);for(n=0;n<c.length;n++)t=c[n],a.indexOf(t)>=0||(s[t]=e[t]);return s}(e,a);if(Object.getOwnPropertySymbols){var c=Object.getOwnPropertySymbols(e);for(n=0;n<c.length;n++)t=c[n],a.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(s[t]=e[t])}return s}var r=s.a.createContext({}),i=function(e){var a=s.a.useContext(r),t=a;return e&&(t="function"==typeof e?e(a):p({},a,{},e)),t},l=function(e){var a=i(e.components);return s.a.createElement(r.Provider,{value:a},e.children)},O={inlineCode:"code",wrapper:function(e){var a=e.children;return s.a.createElement(s.a.Fragment,{},a)}},j=Object(n.forwardRef)((function(e,a){var t=e.components,n=e.mdxType,c=e.originalType,m=e.parentName,r=b(e,["components","mdxType","originalType","parentName"]),l=i(t),j=n,o=l["".concat(m,".").concat(j)]||l[j]||O[j]||c;return t?s.a.createElement(o,p({ref:a},r,{components:t})):s.a.createElement(o,p({ref:a},r))}));function o(e,a){var t=arguments,n=a&&a.mdxType;if("string"==typeof e||n){var c=t.length,m=new Array(c);m[0]=j;var p={};for(var b in a)hasOwnProperty.call(a,b)&&(p[b]=a[b]);p.originalType=e,p.mdxType="string"==typeof e?e:n,m[1]=p;for(var r=2;r<c;r++)m[r]=t[r];return s.a.createElement.apply(null,m)}return s.a.createElement.apply(null,t)}j.displayName="MDXCreateElement"}}]);