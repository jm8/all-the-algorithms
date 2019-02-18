# all-the-algorithms

This project is a work in progress. It will probably never be finished, but if it is, it will be a single python library with on [**The List**](https://en.wikipedia.org/w/index.php?title=List_of_algorithms&oldid=881652734). (It does not include algorithms that are on a seperate list which is linked to on **The List**.)

A checklist directly from Wikipedia, showing our progress:

Combinatorial algorithms
------------------------

### General combinatorial algorithms

- [ ]   [Brent's algorithm](https://en.wikipedia.org/wiki/Brent's_algorithm): finds a cycle in
    function value iterations using only two iterators
- [ ]   [Floyd's cycle-finding
    algorithm](https://en.wikipedia.org/wiki/Floyd's_cycle-finding_algorithm): finds a
    cycle in function value iterations
- [ ]   [Gale–Shapley algorithm](https://en.wikipedia.org/wiki/stable_marriage_problem): solves
    the stable marriage problem
- [ ]   [Pseudorandom number
    generators](https://en.wikipedia.org/wiki/Pseudorandom_number_generator) (uniformly
    distributed):
    - [x]   [Blum Blum Shub](https://en.wikipedia.org/wiki/Blum_Blum_Shub)
    - [x]   [Lagged Fibonacci
        generator](https://en.wikipedia.org/wiki/Lagged_Fibonacci_generator)
    - [x]   [Linear congruential
        generator](https://en.wikipedia.org/wiki/Linear_congruential_generator)
    - [ ]   [Mersenne Twister](https://en.wikipedia.org/wiki/Mersenne_Twister)

### Graph algorithms

- [ ]   [Coloring algorithm](https://en.wikipedia.org/wiki/Coloring_algorithm): Graph coloring
    algorithm.
- [ ]   [Hopcroft–Karp algorithm](https://en.wikipedia.org/wiki/Hopcroft–Karp_algorithm):
    convert a bipartite graph to a [maximum cardinality
    matching](https://en.wikipedia.org/wiki/maximum_cardinality_matching)
- [ ]   [Hungarian algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm): algorithm for
    finding a [perfect matching](https://en.wikipedia.org/wiki/perfect_matching)
- [ ]   [Prüfer coding](https://en.wikipedia.org/wiki/Prüfer_sequence): conversion between a
    labeled tree and its [Prüfer sequence](https://en.wikipedia.org/wiki/Prüfer_sequence)
- [ ]   [Tarjan's off-line lowest common ancestors
    algorithm](https://en.wikipedia.org/wiki/Tarjan's_off-line_lowest_common_ancestors_algorithm):
    compute [lowest common ancestors](https://en.wikipedia.org/wiki/lowest_common_ancestor)
    for pairs of nodes in a tree
- [ ]   [Topological sort](https://en.wikipedia.org/wiki/Topological_sorting): finds linear
    order of nodes (e.g. jobs) based on their dependencies.

#### Graph drawing

- [ ]   [Force-based
    algorithms](Force-based_algorithms_https://en.wikipedia.org/wiki/(graph_drawing)) (also
    known as force-directed algorithms or spring-based algorithm)
- [ ]   [Spectral layout](https://en.wikipedia.org/wiki/Spectral_layout)

#### Network theory

- [ ]   Network analysis
    - [ ]   Link analysis
        - [ ]   [Girvan–Newman
            algorithm](https://en.wikipedia.org/wiki/Girvan–Newman_algorithm): detect
            communities in complex systems
        - [ ]   Web link analysis
            - [ ]   [Hyperlink-Induced Topic
                Search](https://en.wikipedia.org/wiki/Hyperlink-Induced_Topic_Search)
                (HITS) (also known as [Hubs and
                authorities](https://en.wikipedia.org/wiki/Hubs_and_authorities))
            - [ ]   [PageRank](https://en.wikipedia.org/wiki/PageRank)
            - [ ]   [TrustRank](https://en.wikipedia.org/wiki/TrustRank)
- [ ]   [Flow networks](https://en.wikipedia.org/wiki/Flow_network)
    - [ ]   [Dinic's algorithm](https://en.wikipedia.org/wiki/Dinic's_algorithm): is a
        [strongly polynomial](https://en.wikipedia.org/wiki/strongly_polynomial) algorithm
        for computing the [maximum flow](https://en.wikipedia.org/wiki/maximum_flow) in a
        [flow network](https://en.wikipedia.org/wiki/flow_network).
    - [ ]   [Edmonds–Karp algorithm](https://en.wikipedia.org/wiki/Edmonds–Karp_algorithm):
        implementation of Ford–Fulkerson
    - [ ]   [Ford–Fulkerson algorithm](https://en.wikipedia.org/wiki/Ford–Fulkerson_algorithm):
        computes the [maximum flow](https://en.wikipedia.org/wiki/maximum_flow_problem) in
        a graph
    - [ ]   [Karger's algorithm](https://en.wikipedia.org/wiki/Karger's_algorithm): a Monte
        Carlo method to compute the [minimum
        cut](https://en.wikipedia.org/wiki/minimum_cut) of a connected graph
    - [ ]   [Push–relabel algorithm](https://en.wikipedia.org/wiki/Push–relabel_algorithm):
        computes a [maximum flow](https://en.wikipedia.org/wiki/maximum_flow_problem) in a
        graph

#### Routing for graphs

- [ ]   [Edmonds' algorithm](https://en.wikipedia.org/wiki/Edmonds'_algorithm) (also known as
    Chu–Liu/Edmonds' algorithm): find maximum or minimum branchings
- [ ]   [Euclidean minimum spanning
    tree](https://en.wikipedia.org/wiki/Euclidean_minimum_spanning_tree): algorithms for
    computing the minimum spanning tree of a set of points in the plane
- [ ]   [Euclidean shortest path
    problem](https://en.wikipedia.org/wiki/Euclidean_shortest_path_problem): find the
    shortest path between two points that does not intersect any
    obstacle
- [ ]   [Longest path problem](https://en.wikipedia.org/wiki/Longest_path_problem): find a
    simple path of maximum length in a given graph
- [ ]   [Minimum spanning tree](https://en.wikipedia.org/wiki/Minimum_spanning_tree)
    - [ ]   [Borůvka's algorithm](https://en.wikipedia.org/wiki/Borůvka's_algorithm)
    - [ ]   [Kruskal's algorithm](https://en.wikipedia.org/wiki/Kruskal's_algorithm)
    - [ ]   [Prim's algorithm](https://en.wikipedia.org/wiki/Prim's_algorithm)
    - [ ]   [Reverse-delete algorithm](https://en.wikipedia.org/wiki/Reverse-delete_algorithm)
- [ ]   [Nonblocking minimal spanning
    switch](https://en.wikipedia.org/wiki/Nonblocking_minimal_spanning_switch) say, for a
    [telephone exchange](https://en.wikipedia.org/wiki/telephone_exchange)
- [ ]   [Shortest path problem](https://en.wikipedia.org/wiki/Shortest_path_problem)
    - [ ]   [Bellman–Ford algorithm](https://en.wikipedia.org/wiki/Bellman–Ford_algorithm):
        computes [shortest paths](https://en.wikipedia.org/wiki/shortest_path_problem) in a
        weighted graph (where some of the edge weights may be negative)
    - [ ]   [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra's_algorithm):
        computes [shortest paths](https://en.wikipedia.org/wiki/shortest_path_problem) in a
        graph with non-negative edge weights
    - [ ]   [Floyd–Warshall algorithm](https://en.wikipedia.org/wiki/Floyd–Warshall_algorithm):
        solves the [all pairs shortest
        path](https://en.wikipedia.org/wiki/all_pairs_shortest_path) problem in a weighted,
        directed graph
    - [ ]   [Johnson's algorithm](https://en.wikipedia.org/wiki/Johnson's_algorithm): All pairs
        shortest path algorithm in sparse weighted directed graph

<!-- -->

- [ ]   [Transitive closure](https://en.wikipedia.org/wiki/Transitive_closure) problem: find
    the [transitive closure](https://en.wikipedia.org/wiki/transitive_closure) of a given
    binary relation
- [ ]   [Traveling salesman problem](https://en.wikipedia.org/wiki/Traveling_salesman_problem)
    - [ ]   [Christofides algorithm](https://en.wikipedia.org/wiki/Christofides_algorithm)
    - [ ]   [Nearest neighbour
        algorithm](https://en.wikipedia.org/wiki/Nearest_neighbour_algorithm)
- [ ]   [Warnsdorff's rule](https://en.wikipedia.org/wiki/Warnsdorff's_rule): A heuristic
    method for solving the [Knight's tour](https://en.wikipedia.org/wiki/Knight's_tour)
    problem.

#### Graph search

- [ ]   [A\*](https://en.wikipedia.org/wiki/A*): special case of best-first search that uses
    heuristics to improve speed
- [ ]   [B\*](https://en.wikipedia.org/wiki/B*): a best-first graph search algorithm that finds
    the least-cost path from a given initial node to any goal node (out
    of one or more possible goals)
- [ ]   [Backtracking](https://en.wikipedia.org/wiki/Backtracking): abandons partial solutions
    when they are found not to satisfy a complete solution
- [ ]   [Beam search](https://en.wikipedia.org/wiki/Beam_search): is a heuristic search
    algorithm that is an optimization of [best-first
    search](https://en.wikipedia.org/wiki/best-first_search) that reduces its memory
    requirement
- [ ]   [Beam stack search](https://en.wikipedia.org/wiki/Beam_stack_search): integrates
    backtracking with [beam search](https://en.wikipedia.org/wiki/beam_search)
- [ ]   [Best-first search](https://en.wikipedia.org/wiki/Best-first_search): traverses a graph
    in the order of likely importance using a [priority
    queue](https://en.wikipedia.org/wiki/priority_queue)
- [ ]   [Bidirectional search](https://en.wikipedia.org/wiki/Bidirectional_search): find the
    shortest path from an initial vertex to a goal vertex in a directed
    graph
- [ ]   [Breadth-first search](https://en.wikipedia.org/wiki/Breadth-first_search): traverses a
    graph level by level
- [ ]   [Brute-force search](https://en.wikipedia.org/wiki/Brute-force_search): An exhaustive
    and reliable search method, but computationally inefficient in many
    applications.
- [ ]   [D\*](https://en.wikipedia.org/wiki/D*): an [incremental heuristic
    search](https://en.wikipedia.org/wiki/incremental_heuristic_search) algorithm
- [ ]   [Depth-first search](https://en.wikipedia.org/wiki/Depth-first_search): traverses a
    graph branch by branch
- [ ]   [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra's_algorithm): A special
    case of A\* for which no heuristic function is used
- [ ]   [General Problem Solver](https://en.wikipedia.org/wiki/General_Problem_Solver): a
    seminal theorem-proving algorithm intended to work as a universal
    problem solver machine.
- [ ]   [Iterative deepening depth-first
    search](https://en.wikipedia.org/wiki/Iterative_deepening_depth-first_search) (IDDFS):
    a state space search strategy
- [ ]   [Jump point search](https://en.wikipedia.org/wiki/Jump_point_search): An optimization
    to A\* which may reduce computation time by an order of magnitude
    using further heuristics.
- [ ]   [Lexicographic breadth-first
    search](https://en.wikipedia.org/wiki/Lexicographic_breadth-first_search) (also known
    as Lex-BFS): a linear time algorithm for ordering the vertices of a
    graph
- [ ]   [Uniform-cost search](https://en.wikipedia.org/wiki/Uniform-cost_search): a [tree
    search](https://en.wikipedia.org/wiki/Tree_traversal) that finds the lowest-cost route
    where costs vary
- [ ]   [SSS\*](https://en.wikipedia.org/wiki/SSS*): state space search traversing a game tree
    in a best-first fashion similar to that of the A\* search algorithm

#### Subgraphs

- [ ]   [Cliques](Clique_https://en.wikipedia.org/wiki/(graph_theory))
    - [ ]   [Bron–Kerbosch algorithm](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm): a
        technique for finding [maximal
        cliques](https://en.wikipedia.org/wiki/maximal_clique) in an undirected graph
    - [ ]   [MaxCliqueDyn maximum clique
        algorithm](https://en.wikipedia.org/wiki/MaxCliqueDyn_maximum_clique_algorithm):
        find a [maximum clique](https://en.wikipedia.org/wiki/maximum_clique) in an
        undirected graph
- [ ]   [Strongly connected
    components](https://en.wikipedia.org/wiki/Strongly_connected_components)
    - [ ]   [Path-based strong component
        algorithm](https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm)
    - [ ]   [Kosaraju's algorithm](https://en.wikipedia.org/wiki/Kosaraju's_algorithm)
    - [ ]   [Tarjan's strongly connected components
        algorithm](https://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm)
- [ ]   [Subgraph isomorphism
    problem](https://en.wikipedia.org/wiki/Subgraph_isomorphism_problem)


#### Approximate sequence matching

- [ ]   [Bitap algorithm](https://en.wikipedia.org/wiki/Bitap_algorithm): fuzzy algorithm that
    determines if strings are approximately equal.
- [ ]   [Phonetic algorithms](https://en.wikipedia.org/wiki/Phonetic_algorithm)
    - [ ]   [Daitch–Mokotoff Soundex](https://en.wikipedia.org/wiki/Daitch–Mokotoff_Soundex): a
        [Soundex](https://en.wikipedia.org/wiki/Soundex) refinement which allows matching
        of Slavic and Germanic surnames
    - [ ]   [Double Metaphone](https://en.wikipedia.org/wiki/Double_Metaphone): an improvement
        on Metaphone
    - [ ]   [Match rating approach](https://en.wikipedia.org/wiki/Match_rating_approach): a
        phonetic algorithm developed by Western Airlines
    - [ ]   [Metaphone](https://en.wikipedia.org/wiki/Metaphone): an algorithm for indexing
        words by their sound, when pronounced in English
    - [ ]   [NYSIIS](https://en.wikipedia.org/wiki/New_York_State_Identification_and_Intelligence_System):
        [phonetic algorithm](https://en.wikipedia.org/wiki/phonetic_algorithm), improves on
        [Soundex](https://en.wikipedia.org/wiki/Soundex)
    - [ ]   [Soundex](https://en.wikipedia.org/wiki/Soundex): a phonetic algorithm for indexing
        names by sound, as pronounced in English
- [ ]   [String metrics](https://en.wikipedia.org/wiki/String_metric): compute a similarity or
    dissimilarity (distance) score between two pairs of text strings
    - [ ]   [Damerau–Levenshtein
        distance](https://en.wikipedia.org/wiki/Damerau–Levenshtein_distance) compute a
        distance measure between two strings, improves on [Levenshtein
        distance](https://en.wikipedia.org/wiki/Levenshtein_distance)
    - [ ]   [Dice's coefficient](https://en.wikipedia.org/wiki/Dice's_coefficient) (also known
        as the Dice coefficient): a similarity measure related to the
        [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)
    - [ ]   [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance): sum number of
        positions which are different
    - [ ]   [Jaro–Winkler distance](https://en.wikipedia.org/wiki/Jaro–Winkler_distance): is a
        measure of similarity between two strings
    - [ ]   [Levenshtein edit distance](https://en.wikipedia.org/wiki/Levenshtein_distance):
        compute a metric for the amount of difference between two
        sequences
- [ ]   [Trigram search](https://en.wikipedia.org/wiki/Trigram_search): search for text when
    the exact syntax or spelling of the target object is not precisely
    known

#### Selection algorithms

- [ ]   [Quickselect](https://en.wikipedia.org/wiki/Quickselect)
- [ ]   [Introselect](https://en.wikipedia.org/wiki/Introselect)

#### Sequence search

- [ ]   [Linear search](https://en.wikipedia.org/wiki/Linear_search): finds an item in an
    unsorted sequence
- [ ]   [Selection algorithm](https://en.wikipedia.org/wiki/Selection_algorithm): finds the
    *k*th largest item in a sequence
- [ ]   [Ternary search](https://en.wikipedia.org/wiki/Ternary_search): a technique for finding
    the minimum or maximum of a function that is either strictly
    increasing and then strictly decreasing or vice versa
- [ ]   [Sorted lists](https://en.wikipedia.org/wiki/Sorted_list)
    - [ ]   [Binary search algorithm](https://en.wikipedia.org/wiki/Binary_search_algorithm):
        locates an item in a sorted sequence
    - [ ]   [Fibonacci search
        technique](https://en.wikipedia.org/wiki/Fibonacci_search_technique): search a
        sorted sequence using a divide and conquer algorithm that
        narrows down possible locations with the aid of [Fibonacci
        numbers](https://en.wikipedia.org/wiki/Fibonacci_numbers)
    - [ ]   [Jump search](https://en.wikipedia.org/wiki/Jump_search) (or block search): linear
        search on a smaller subset of the sequence
    - [ ]   [Predictive search](https://en.wikipedia.org/wiki/Interpolation_search):
        binary-like search which factors in
        [magnitude](magnitude_https://en.wikipedia.org/wiki/(mathematics)) of search term
        versus the high and low values in the search. Sometimes called
        dictionary search or interpolated search.
    - [ ]   [Uniform binary search](https://en.wikipedia.org/wiki/Uniform_binary_search): an
        optimization of the classic binary search algorithm

#### Sequence merging

- [ ]   Simple merge algorithm
- [ ]   k-way merge algorithm
- [ ]   Union (merge, with elements on the output not repeated)

#### Sequence permutations

- [ ]   [Fisher–Yates shuffle](https://en.wikipedia.org/wiki/Fisher–Yates_shuffle) (also known
    as the Knuth shuffle): randomly shuffle a finite set
- [ ]   [Schensted algorithm](https://en.wikipedia.org/wiki/Schensted_algorithm): constructs a
    pair of [Young tableaux](https://en.wikipedia.org/wiki/Young_tableau) from a
    permutation
- [ ]   [Steinhaus–Johnson–Trotter
    algorithm](https://en.wikipedia.org/wiki/Steinhaus–Johnson–Trotter_algorithm) (also
    known as the Johnson–Trotter algorithm): generate permutations by
    transposing elements
- [ ]   [Heap's permutation generation
    algorithm](https://en.wikipedia.org/wiki/Heap's_algorithm): interchange elements to
    generate next permutation

#### Sequence alignment

- [ ]   [Dynamic time warping](https://en.wikipedia.org/wiki/Dynamic_time_warping): measure
    similarity between two sequences which may vary in time or speed
- [ ]   [Hirschberg's algorithm](https://en.wikipedia.org/wiki/Hirschberg's_algorithm): finds
    the least cost [sequence alignment](https://en.wikipedia.org/wiki/sequence_alignment)
    between two sequences, as measured by their [Levenshtein
    distance](https://en.wikipedia.org/wiki/Levenshtein_distance)
- [ ]   [Needleman–Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm):
    find global alignment between two sequences
- [ ]   [Smith–Waterman algorithm](https://en.wikipedia.org/wiki/Smith–Waterman_algorithm):
    find local sequence alignment

#### Sequence sorting

- [x]   Exchange sorts
    - [x]   [Bubble sort](https://en.wikipedia.org/wiki/Bubble_sort): for each pair of indices,
        swap the items if out of order
    - [x]   [Cocktail shaker sort](https://en.wikipedia.org/wiki/Cocktail_shaker_sort) or
        bidirectional bubble sort, a bubble sort traversing the list
        alternately from front to back and back to front
    - [x]   [Comb sort](https://en.wikipedia.org/wiki/Comb_sort)
    - [x]   [Gnome sort](https://en.wikipedia.org/wiki/Gnome_sort)
    - [x]   [Odd–even sort](https://en.wikipedia.org/wiki/Odd–even_sort)
    - [x]   [Quicksort](https://en.wikipedia.org/wiki/Quicksort): divide list into two, with
        all items on the first list coming before all items on the
        second list.; then sort the two lists. Often the method of
        choice
- [x]   Humorous or ineffective
    - [x]   [Bogosort](https://en.wikipedia.org/wiki/Bogosort)
    - [x]   [Stooge sort](https://en.wikipedia.org/wiki/Stooge_sort)
- [ ]   Hybrid
    - [ ]   [Flashsort](https://en.wikipedia.org/wiki/Flashsort)
    - [ ]   [Introsort](https://en.wikipedia.org/wiki/Introsort): begin with quicksort and
        switch to heapsort when the recursion depth exceeds a certain
        level
    - [ ]   [Timsort](https://en.wikipedia.org/wiki/Timsort): adaptative algorithm derived from
        merge sort and insertion sort. Used in Python 2.3 and up, and
        Java SE 7.
- [ ]   Insertion sorts
    - [ ]   [Insertion sort](https://en.wikipedia.org/wiki/Insertion_sort): determine where the
        current item belongs in the list of sorted ones, and insert it
        there
    - [ ]   [Library sort](https://en.wikipedia.org/wiki/Library_sort)
    - [ ]   [Patience sorting](https://en.wikipedia.org/wiki/Patience_sorting)
    - [ ]   [Shell sort](https://en.wikipedia.org/wiki/Shellsort): an attempt to improve
        insertion sort
    - [ ]   [Tree sort](https://en.wikipedia.org/wiki/Tree_sort) (binary tree sort): build
        binary tree, then traverse it to create sorted list
    - [ ]   [Cycle sort](https://en.wikipedia.org/wiki/Cycle_sort): in-place with theoretically
        optimal number of writes
- [ ]   Merge sorts
    - [ ]   [Merge sort](https://en.wikipedia.org/wiki/Merge_sort): sort the first and second
        half of the list separately, then merge the sorted lists
    - [ ]   [Strand sort](https://en.wikipedia.org/wiki/Strand_sort)
- [ ]   Non-comparison sorts
    - [ ]   [Bead sort](https://en.wikipedia.org/wiki/Bead_sort)
    - [ ]   [Bucket sort](https://en.wikipedia.org/wiki/Bucket_sort)
    - [ ]   [Burstsort](https://en.wikipedia.org/wiki/Burstsort): build a compact, cache
        efficient [burst trie](https://en.wikipedia.org/wiki/burst_trie) and then traverse
        it to create sorted output
    - [ ]   [Counting sort](https://en.wikipedia.org/wiki/Counting_sort)
    - [ ]   [Pigeonhole sort](https://en.wikipedia.org/wiki/Pigeonhole_sort)
    - [ ]   [Postman sort](https://en.wikipedia.org/wiki/Postman_sort): variant of Bucket sort
        which takes advantage of hierarchical structure
    - [ ]   [Radix sort](https://en.wikipedia.org/wiki/Radix_sort): sorts strings letter by
        letter
- [ ]   Selection sorts
    - [ ]   [Heapsort](https://en.wikipedia.org/wiki/Heapsort): convert the list into a heap,
        keep removing the largest element from the heap and adding it to
        the end of the list
    - [ ]   [Selection sort](https://en.wikipedia.org/wiki/Selection_sort): pick the smallest
        of the remaining elements, add it to the end of the sorted list
    - [ ]   [Smoothsort](https://en.wikipedia.org/wiki/Smoothsort)
- [ ]   Other
    - [ ]   [Bitonic sorter](https://en.wikipedia.org/wiki/Bitonic_sorter)
    - [ ]   [Pancake sorting](https://en.wikipedia.org/wiki/Pancake_sorting)
    - [ ]   [Spaghetti sort](https://en.wikipedia.org/wiki/Spaghetti_sort)
    - [ ]   [Topological sort](https://en.wikipedia.org/wiki/Topological_sorting)
- [ ]   Unknown class
    - [ ]   [Samplesort](https://en.wikipedia.org/wiki/Samplesort)

#### Subsequences

- [ ]   [Kadane's algorithm](https://en.wikipedia.org/wiki/Kadane's_algorithm): finds maximum
    sub-array of any size
- [ ]   [Longest common subsequence
    problem](https://en.wikipedia.org/wiki/Longest_common_subsequence_problem): Find the
    longest subsequence common to all sequences in a set of sequences
- [ ]   [Longest increasing subsequence
    problem](https://en.wikipedia.org/wiki/Longest_increasing_subsequence_problem): Find
    the longest increasing subsequence of a given sequence
- [ ]   [Shortest common supersequence
    problem](https://en.wikipedia.org/wiki/Shortest_common_supersequence_problem): Find the
    shortest supersequence that contains two or more sequences as
    subsequences

#### Substrings

- [ ]   [Longest common substring
    problem](https://en.wikipedia.org/wiki/Longest_common_substring_problem): find the
    longest string (or strings) that is a substring (or are substrings)
    of two or more strings
- [ ]   [Substring search](https://en.wikipedia.org/wiki/Substring_search)
    - [ ]   [Aho–Corasick string matching
        algorithm](https://en.wikipedia.org/wiki/Aho–Corasick_string_matching_algorithm):
        [trie](https://en.wikipedia.org/wiki/trie) based algorithm for finding all
        substring matches to any of a finite set of strings
    - [ ]   [Boyer–Moore string-search
        algorithm](https://en.wikipedia.org/wiki/Boyer–Moore_string-search_algorithm):
        amortized linear ([sublinear](https://en.wikipedia.org/wiki/sublinear) in most
        times) algorithm for substring search
    - [ ]   [Boyer–Moore–Horspool
        algorithm](https://en.wikipedia.org/wiki/Boyer–Moore–Horspool_algorithm):
        Simplification of Boyer–Moore
    - [ ]   [Knuth–Morris–Pratt
        algorithm](https://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm): substring
        search which bypasses reexamination of matched characters
    - [ ]   [Rabin–Karp string search
        algorithm](https://en.wikipedia.org/wiki/Rabin–Karp_string_search_algorithm):
        searches multiple patterns efficiently
    - [ ]   [Zhu–Takaoka string matching
        algorithm](https://en.wikipedia.org/wiki/Zhu–Takaoka_string_matching_algorithm): a
        variant of Boyer–Moore
- [ ]   [Ukkonen's algorithm](https://en.wikipedia.org/wiki/Ukkonen's_algorithm): a
    [linear-time](https://en.wikipedia.org/wiki/linear-time), [online
    algorithm](https://en.wikipedia.org/wiki/online_algorithm) for constructing [suffix
    trees](https://en.wikipedia.org/wiki/suffix_tree)
- [ ]   [Matching wildcards](https://en.wikipedia.org/wiki/Matching_wildcards)
    - [ ]   [Rich Salz](https://en.wikipedia.org/wiki/InterNetNews)'
        [wildmat](https://en.wikipedia.org/wiki/wildmat): a widely used
        [open-source](https://en.wikipedia.org/wiki/Open-source_software)
        [recursive](https://en.wikipedia.org/wiki/recursion) algorithm
    - [ ]   [Krauss matching wildcards
        algorithm](https://en.wikipedia.org/wiki/Krauss_matching_wildcards_algorithm): an
        open-source non-recursive algorithm

Computational mathematics
-------------------------

### Abstract algebra

- [ ]   [Chien search](https://en.wikipedia.org/wiki/Chien_search): a recursive algorithm for
    determining roots of polynomials defined over a finite field
- [ ]   [Schreier–Sims algorithm](https://en.wikipedia.org/wiki/Schreier–Sims_algorithm):
    computing a base and [strong generating
    set](https://en.wikipedia.org/wiki/strong_generating_set) (BSGS) of a [permutation
    group](https://en.wikipedia.org/wiki/permutation_group)
- [ ]   [Todd–Coxeter algorithm](https://en.wikipedia.org/wiki/Todd–Coxeter_algorithm):
    Procedure for generating [cosets](https://en.wikipedia.org/wiki/coset).

### Computer algebra

- [ ]   [Buchberger's algorithm](https://en.wikipedia.org/wiki/Buchberger's_algorithm): finds a
    [Gröbner basis](https://en.wikipedia.org/wiki/Gröbner_basis)
- [ ]   [Cantor–Zassenhaus
    algorithm](https://en.wikipedia.org/wiki/Cantor–Zassenhaus_algorithm): factor
    polynomials over finite fields
- [ ]   [Faugère F4 algorithm](https://en.wikipedia.org/wiki/Faugère_F4_algorithm): finds a
    Gröbner basis (also mentions the F5 algorithm)
- [ ]   [Gosper's algorithm](https://en.wikipedia.org/wiki/Gosper's_algorithm): find sums of
    hypergeometric terms that are themselves hypergeometric terms
- [ ]   [Knuth–Bendix completion
    algorithm](https://en.wikipedia.org/wiki/Knuth–Bendix_completion_algorithm): for
    [rewriting](https://en.wikipedia.org/wiki/rewriting) rule systems
- [ ]   [Multivariate division
    algorithm](https://en.wikipedia.org/wiki/Multivariate_division_algorithm): for
    [polynomials](https://en.wikipedia.org/wiki/polynomial) in several indeterminates
- [ ]   [Pollard's kangaroo
    algorithm](https://en.wikipedia.org/wiki/Pollard's_kangaroo_algorithm) (also known as
    Pollard's lambda algorithm ): an algorithm for solving the discrete
    logarithm problem
- [ ]   [Polynomial long division](https://en.wikipedia.org/wiki/Polynomial_long_division): an
    algorithm for dividing a polynomial by another polynomial of the
    same or lower degree
- [ ]   [Risch algorithm](https://en.wikipedia.org/wiki/Risch_algorithm): an algorithm for the
    calculus operation of indefinite integration (i.e. finding
    [antiderivatives](https://en.wikipedia.org/wiki/antiderivatives))

### Geometry

- [ ]   [Closest pair problem](https://en.wikipedia.org/wiki/Closest_pair_problem): find the
    pair of points (from a set of points) with the smallest distance
    between them
- [ ]   [Collision detection](https://en.wikipedia.org/wiki/Collision_detection) algorithms:
    check for the collision or intersection of two given solids
- [ ]   [Cone algorithm](https://en.wikipedia.org/wiki/Cone_algorithm): identify surface points
- [ ]   [Convex hull algorithms](https://en.wikipedia.org/wiki/Convex_hull_algorithms):
    determining the [convex hull](https://en.wikipedia.org/wiki/convex_hull) of a
    [set](Set_https://en.wikipedia.org/wiki/(mathematics)) of points
    - [ ]   [Graham scan](https://en.wikipedia.org/wiki/Graham_scan)
    - [ ]   [Quickhull](https://en.wikipedia.org/wiki/Quickhull)
    - [ ]   [Gift wrapping algorithm](https://en.wikipedia.org/wiki/Gift_wrapping_algorithm) or
        Jarvis march
    - [ ]   [Chan's algorithm](https://en.wikipedia.org/wiki/Chan's_algorithm)
    - [ ]   [Kirkpatrick–Seidel
        algorithm](https://en.wikipedia.org/wiki/Kirkpatrick–Seidel_algorithm)
- [ ]   [Euclidean Distance Transform](https://en.wikipedia.org/wiki/Euclidean_distance_map) -
    Computes the distance between every point in a grid and a discrete
    collection of points.
- [ ]   [Geometric hashing](https://en.wikipedia.org/wiki/Geometric_hashing): a method for
    efficiently finding two-dimensional objects represented by discrete
    points that have undergone an [affine
    transformation](https://en.wikipedia.org/wiki/affine_transformation)
- [ ]   [Gilbert–Johnson–Keerthi distance
    algorithm](https://en.wikipedia.org/wiki/Gilbert–Johnson–Keerthi_distance_algorithm):
    determining the smallest distance between two
    [convex](https://en.wikipedia.org/wiki/convex_set) shapes.
- [ ]   [Jump-and-Walk algorithm](https://en.wikipedia.org/wiki/Jump-and-Walk_algorithm): an
    algorithm for point location in triangulations
- [ ]   [Laplacian smoothing](https://en.wikipedia.org/wiki/Laplacian_smoothing): an algorithm
    to smooth a polygonal mesh
- [ ]   [Line segment intersection](https://en.wikipedia.org/wiki/Line_segment_intersection):
    finding whether lines intersect, usually with a [sweep line
    algorithm](https://en.wikipedia.org/wiki/sweep_line_algorithm)
    - [ ]   [Bentley–Ottmann
        algorithm](https://en.wikipedia.org/wiki/Bentley–Ottmann_algorithm)
    - [ ]   [Shamos–Hoey algorithm](https://en.wikipedia.org/wiki/Shamos–Hoey_algorithm)
- [ ]   [Minimum bounding box
    algorithms](https://en.wikipedia.org/wiki/Minimum_bounding_box_algorithms): find the
    [oriented minimum bounding
    box](https://en.wikipedia.org/wiki/Minimum_bounding_box#Arbitrarily_oriented_minimum_bounding_box)
    enclosing a set of points
- [ ]   [Nearest neighbor search](https://en.wikipedia.org/wiki/Nearest_neighbor_search): find
    the nearest point or points to a query point
- [ ]   [Point in polygon](https://en.wikipedia.org/wiki/Point_in_polygon) algorithms: tests
    whether a given point lies within a given polygon
- [ ]   [Point set registration](https://en.wikipedia.org/wiki/Point_set_registration)
    algorithms: finds the transformation between two [point
    sets](https://en.wikipedia.org/wiki/point_cloud) to optimally align them.
- [ ]   [Rotating calipers](https://en.wikipedia.org/wiki/Rotating_calipers): determine all
    [antipodal](https://en.wikipedia.org/wiki/antipodal_point) pairs of points and vertices
    on a [convex polygon](https://en.wikipedia.org/wiki/convex_polygon) or [convex
    hull](https://en.wikipedia.org/wiki/convex_hull).
- [ ]   [Shoelace algorithm](https://en.wikipedia.org/wiki/Shoelace_algorithm): determine the
    area of a polygon whose vertices are described by ordered pairs in
    the plane
- [ ]   [Triangulation](Triangulation_https://en.wikipedia.org/wiki/(geometry))
    - [ ]   [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation)
        - [ ]   [Ruppert's algorithm](https://en.wikipedia.org/wiki/Ruppert's_algorithm) (also
            known as Delaunay refinement): create quality Delaunay
            triangulations
        - [ ]   [Chew's second
            algorithm](https://en.wikipedia.org/wiki/Chew's_second_algorithm): create
            quality [constrained Delaunay
            triangulations](https://en.wikipedia.org/wiki/constrained_Delaunay_triangulation)
    - [ ]   [Marching triangles](https://en.wikipedia.org/wiki/Marching_triangles): reconstruct
        two-dimensional surface geometry from an unstructured [point
        cloud](https://en.wikipedia.org/wiki/point_cloud)
    - [ ]   [Polygon triangulation](https://en.wikipedia.org/wiki/Polygon_triangulation)
        algorithms: decompose a polygon into a set of triangles
    - [ ]   [Voronoi diagrams](https://en.wikipedia.org/wiki/Voronoi_diagram), geometric
        [dual](duality_https://en.wikipedia.org/wiki/(mathematics)) of [Delaunay
        triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation)
        - [ ]   [Bowyer–Watson
            algorithm](https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm): create
            voronoi diagram in any number of dimensions
        - [ ]   [Fortune's Algorithm](https://en.wikipedia.org/wiki/Fortune's_Algorithm):
            create voronoi diagram
    - [ ]   [Quasitriangulation](https://en.wikipedia.org/wiki/Quasitriangulation)

### Number theoretic algorithms

- [x]   [Binary GCD algorithm](https://en.wikipedia.org/wiki/Binary_GCD_algorithm): Efficient
    way of calculating GCD.
- [ ]   [Booth's multiplication
    algorithm](https://en.wikipedia.org/wiki/Booth's_multiplication_algorithm)
- [ ]   [Chakravala method](https://en.wikipedia.org/wiki/Chakravala_method): a cyclic
    algorithm to solve indeterminate quadratic equations, including
    [Pell's equation](https://en.wikipedia.org/wiki/Pell's_equation)
- [ ]   [Discrete logarithm](https://en.wikipedia.org/wiki/Discrete_logarithm):
    - [ ]   [Baby-step giant-step](https://en.wikipedia.org/wiki/Baby-step_giant-step)
    - [ ]   [Index calculus algorithm](https://en.wikipedia.org/wiki/Index_calculus_algorithm)
    - [ ]   [Pollard's rho algorithm for
        logarithms](https://en.wikipedia.org/wiki/Pollard's_rho_algorithm_for_logarithms)
    - [ ]   [Pohlig&ndash;Hellman
        algorithm](https://en.wikipedia.org/wiki/Pohlig&ndash;Hellman_algorithm)
- [ ]   [Euclidean algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm): computes the
    [greatest common divisor](https://en.wikipedia.org/wiki/greatest_common_divisor)
- [ ]   [Extended Euclidean
    algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm): Also solves the
    equation *ax* + *by* = *c*.
- [ ]   [Integer factorization](https://en.wikipedia.org/wiki/Integer_factorization): breaking
    an integer into its [prime](https://en.wikipedia.org/wiki/prime_number) factors
    - [ ]   [Congruence of squares](https://en.wikipedia.org/wiki/Congruence_of_squares)
    - [ ]   [Dixon's algorithm](https://en.wikipedia.org/wiki/Dixon's_algorithm)
    - [ ]   [Fermat's factorization
        method](https://en.wikipedia.org/wiki/Fermat's_factorization_method)
    - [ ]   [General number field
        sieve](https://en.wikipedia.org/wiki/General_number_field_sieve)
    - [ ]   [Lenstra elliptic curve
        factorization](https://en.wikipedia.org/wiki/Lenstra_elliptic_curve_factorization)
    - [ ]   [Pollard's *p* − 1
        algorithm](https://en.wikipedia.org/wiki/Pollard's_p_&minus;_1_algorithm)
    - [ ]   [Pollard's rho algorithm](https://en.wikipedia.org/wiki/Pollard's_rho_algorithm)
    - [ ]   [prime factorization
        algorithm](https://en.wikipedia.org/wiki/prime_factorization_algorithm)
    - [ ]   [Quadratic sieve](https://en.wikipedia.org/wiki/Quadratic_sieve)
    - [ ]   [Shor's algorithm](https://en.wikipedia.org/wiki/Shor's_algorithm)
    - [ ]   [Special number field
        sieve](https://en.wikipedia.org/wiki/Special_number_field_sieve)
    - [ ]   [Trial division](https://en.wikipedia.org/wiki/Trial_division)
- [ ]   [Multiplication algorithms](https://en.wikipedia.org/wiki/Multiplication_algorithm):
    fast multiplication of two numbers
    - [ ]   [Karatsuba algorithm](https://en.wikipedia.org/wiki/Karatsuba_algorithm)
    - [ ]   [Schönhage–Strassen
        algorithm](https://en.wikipedia.org/wiki/Schönhage–Strassen_algorithm)
    - [ ]   [Toom–Cook multiplication](https://en.wikipedia.org/wiki/Toom–Cook_multiplication)
- [ ]   [Modular square root](https://en.wikipedia.org/wiki/Modular_square_root): computing
    square roots modulo a prime number
    - [ ]   [Tonelli–Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli–Shanks_algorithm)
    - [ ]   [Cipolla's algorithm](https://en.wikipedia.org/wiki/Cipolla's_algorithm)
- [ ]   [Odlyzko&ndash;Schönhage
    algorithm](https://en.wikipedia.org/wiki/Odlyzko&ndash;Schönhage_algorithm): calculates
    nontrivial zeroes of the [Riemann zeta
    function](https://en.wikipedia.org/wiki/Riemann_zeta_function)
- [ ]   [Lenstra–Lenstra–Lovász
    algorithm](https://en.wikipedia.org/wiki/Lenstra–Lenstra–Lovász_lattice_basis_reduction_algorithm)
    (also known as LLL algorithm): find a short, nearly orthogonal
    [lattice](Lattice_https://en.wikipedia.org/wiki/(group))
    [basis](Basis_https://en.wikipedia.org/wiki/(linear_algebra)) in polynomial time
- [ ]   [Primality tests](https://en.wikipedia.org/wiki/Primality_test): determining whether a
    given number is [prime](https://en.wikipedia.org/wiki/prime_number)
    - [ ]   [AKS primality test](https://en.wikipedia.org/wiki/AKS_primality_test)
    - [ ]   [Baillie-PSW primality
        test](https://en.wikipedia.org/wiki/Baillie-PSW_primality_test)
    - [ ]   [Fermat primality test](https://en.wikipedia.org/wiki/Fermat_primality_test)
    - [ ]   [Lucas primality test](https://en.wikipedia.org/wiki/Lucas_primality_test)
    - [ ]   [Miller&ndash;Rabin primality
        test](https://en.wikipedia.org/wiki/Miller&ndash;Rabin_primality_test)
    - [ ]   [Sieve of Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin)
    - [ ]   [Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes)
    - [ ]   [Sieve of Sundaram](https://en.wikipedia.org/wiki/Sieve_of_Sundaram)

### Numerical algorithms

#### Differential equation solving

- [ ]   [Euler method](https://en.wikipedia.org/wiki/Euler_method)
- [ ]   [Backward Euler method](https://en.wikipedia.org/wiki/Backward_Euler_method)
- [ ]   [Trapezoidal rule (differential
    equations)](Trapezoidal_rule_https://en.wikipedia.org/wiki/(differential_equations))
- [ ]   [Linear multistep methods](https://en.wikipedia.org/wiki/Linear_multistep_method)
- [ ]   [Runge–Kutta methods](https://en.wikipedia.org/wiki/Runge–Kutta_methods)
    - [ ]   [Euler integration](https://en.wikipedia.org/wiki/Euler_integration)
- [ ]   [Multigrid methods](https://en.wikipedia.org/wiki/Multigrid_method) (MG methods), a
    group of algorithms for solving differential equations using a
    hierarchy of discretizations
- [ ]   [Partial differential
    equation](https://en.wikipedia.org/wiki/Partial_differential_equation):
    - [ ]   [Finite difference method](https://en.wikipedia.org/wiki/Finite_difference_method)
    - [ ]   [Crank–Nicolson method](https://en.wikipedia.org/wiki/Crank–Nicolson_method) for
        diffusion equations
    - [ ]   [Lax-Wendroff](https://en.wikipedia.org/wiki/Lax–Wendroff_method) for wave
        equations
- [ ]   [Verlet integration](https://en.wikipedia.org/wiki/Verlet_integration) (): integrate
    Newton's equations of motion

#### Elementary and special functions

- [ ]   [Computation of π](https://en.wikipedia.org/wiki/Computing_π):
    - [ ]   [Borwein's algorithm](https://en.wikipedia.org/wiki/Borwein's_algorithm): an
        algorithm to calculate the value of 1/π
    - [ ]   [Gauss–Legendre algorithm](https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm):
        computes the digits of [pi](https://en.wikipedia.org/wiki/pi)
    - [ ]   [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm): A fast
        method for calculating the digits of π
    - [ ]   [Bailey–Borwein–Plouffe
        formula](https://en.wikipedia.org/wiki/Bailey–Borwein–Plouffe_formula): (BBP
        formula) a spigot algorithm for the computation of the nth
        binary digit of π
- [ ]   [Division algorithms](https://en.wikipedia.org/wiki/Division_algorithm): for computing
    quotient and/or remainder of two numbers
    - [ ]   [Long division](https://en.wikipedia.org/wiki/Long_division)
    - [ ]   [Restoring division](https://en.wikipedia.org/wiki/Restoring_division)
    - [ ]   [Non-restoring division](https://en.wikipedia.org/wiki/Non-restoring_division)
    - [ ]   [SRT division](https://en.wikipedia.org/wiki/SRT_division)
    - [ ]   [Newton–Raphson division](https://en.wikipedia.org/wiki/Newton–Raphson_division):
        uses [Newton's method](https://en.wikipedia.org/wiki/Newton's_method) to find the
        [reciprocal](https://en.wikipedia.org/wiki/Multiplicative_inverse) of D, and
        multiply that reciprocal by N to find the final quotient Q.
    - [ ]   [Goldschmidt division](https://en.wikipedia.org/wiki/Goldschmidt_division)
- [ ]   Hyperbolic and Trigonometric Functions:
    - [ ]   [BKM algorithm](https://en.wikipedia.org/wiki/BKM_algorithm): compute [elementary
        functions](Elementary_function_https://en.wikipedia.org/wiki/(differential_algebra))
        using a table of logarithms
    - [ ]   [CORDIC](https://en.wikipedia.org/wiki/CORDIC): compute hyperbolic and
        trigonometric functions using a table of arctangents
- [ ]   Exponentiation:
    - [ ]   [Addition-chain
        exponentiation](https://en.wikipedia.org/wiki/Addition-chain_exponentiation)
        exponentiation by positive integer powers that requires a
        minimal number of multiplications
    - [ ]   [Exponentiating by
        squaring](https://en.wikipedia.org/wiki/Exponentiating_by_squaring): an algorithm
        used for the fast computation of [large
        integer](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic) powers of a
        number
- [ ]   [Montgomery reduction](https://en.wikipedia.org/wiki/Montgomery_reduction): an
    algorithm that allows [modular
    arithmetic](https://en.wikipedia.org/wiki/modular_arithmetic) to be performed
    efficiently when the modulus is large
- [ ]   [Multiplication algorithms](https://en.wikipedia.org/wiki/Multiplication_algorithm):
    fast multiplication of two numbers
    - [ ]   [Booth's multiplication
        algorithm](https://en.wikipedia.org/wiki/Booth's_multiplication_algorithm): a
        multiplication algorithm that multiplies two signed binary
        numbers in two's complement notation
    - [ ]   [Fürer's algorithm](https://en.wikipedia.org/wiki/Fürer's_algorithm): an integer
        multiplication algorithm for very large numbers possessing a
        very low [asymptotic
        complexity](https://en.wikipedia.org/wiki/Computational_complexity_theory)
    - [ ]   [Karatsuba algorithm](https://en.wikipedia.org/wiki/Karatsuba_algorithm): an
        efficient procedure for multiplying large numbers
    - [ ]   [Schönhage–Strassen
        algorithm](https://en.wikipedia.org/wiki/Schönhage–Strassen_algorithm): an
        asymptotically fast multiplication algorithm for large integers
    - [ ]   [Toom–Cook multiplication](https://en.wikipedia.org/wiki/Toom–Cook_multiplication):
        (Toom3) a multiplication algorithm for large integers
- [ ]   [Multiplicative inverse
    Algorithms](https://en.wikipedia.org/wiki/Multiplicative_inverse#Algorithms): for
    computing a number's multiplicative inverse (reciprocal).
    - [ ]   [Newton's
        method](https://en.wikipedia.org/wiki/Newton's_method#Multiplicative_inverses_of_numbers_and_power_series)
- [ ]   [Rounding functions](https://en.wikipedia.org/wiki/Rounding_functions): the classic
    ways to round numbers
- [ ]   [Spigot algorithm](https://en.wikipedia.org/wiki/Spigot_algorithm): A way to compute
    the value of a [mathematical
    constant](https://en.wikipedia.org/wiki/mathematical_constant) without knowing
    preceding digits
- [ ]   Square and Nth root of a number:
    - [ ]   [Alpha max plus beta min
        algorithm](https://en.wikipedia.org/wiki/Alpha_max_plus_beta_min_algorithm): an
        approximation of the square-root of the sum of two squares
    - [ ]   [Methods of computing square
        roots](https://en.wikipedia.org/wiki/Methods_of_computing_square_roots)
    - [ ]   [*n*th root algorithm](https://en.wikipedia.org/wiki/Nth_root_algorithm)
    - [ ]   [Shifting nth-root
        algorithm](https://en.wikipedia.org/wiki/Shifting_nth-root_algorithm): digit by
        digit root extraction
- [ ]   Summation:
    - [ ]   [Binary splitting](https://en.wikipedia.org/wiki/Binary_splitting): a [divide and
        conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithm) technique
        which speeds up the numerical evaluation of many types of series
        with rational terms
    - [ ]   [Kahan summation
        algorithm](https://en.wikipedia.org/wiki/Kahan_summation_algorithm): a more
        accurate method of summing floating-point numbers
- [ ]   [Unrestricted algorithm](https://en.wikipedia.org/wiki/Unrestricted_algorithm)

#### Geometric

- [ ]   [Filtered
    back-projection](https://en.wikipedia.org/wiki/Radon_transform#Filtered_back-projection):
    efficiently compute the inverse 2-dimensional [Radon
    transform](https://en.wikipedia.org/wiki/Radon_transform).
- [ ]   [Level set method](https://en.wikipedia.org/wiki/Level_set_method) (LSM): a numerical
    technique for tracking interfaces and shapes

#### Interpolation and extrapolation

- [ ]   [Birkhoff interpolation](https://en.wikipedia.org/wiki/Birkhoff_interpolation): an
    extension of polynomial interpolation
- [ ]   [Cubic interpolation](https://en.wikipedia.org/wiki/Cubic_interpolation)
- [ ]   [Hermite interpolation](https://en.wikipedia.org/wiki/Hermite_interpolation)
- [ ]   [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_interpolation):
    interpolation using [Lagrange
    polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial)
- [ ]   [Linear interpolation](https://en.wikipedia.org/wiki/Linear_interpolation): a method of
    curve fitting using linear polynomials
- [ ]   [Monotone cubic
    interpolation](https://en.wikipedia.org/wiki/Monotone_cubic_interpolation): a variant
    of cubic interpolation that preserves monotonicity of the data set
    being interpolated.
- [ ]   [Multivariate interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation)
    - [ ]   [Bicubic interpolation](https://en.wikipedia.org/wiki/Bicubic_interpolation), a
        generalization of [cubic
        interpolation](https://en.wikipedia.org/wiki/cubic_interpolation) to two dimensions
    - [ ]   [Bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation): an
        extension of [linear
        interpolation](https://en.wikipedia.org/wiki/linear_interpolation) for
        interpolating functions of two variables on a regular grid
    - [ ]   [Lanczos resampling](https://en.wikipedia.org/wiki/Lanczos_resampling) (“Lanzosh”):
        a multivariate interpolation method used to compute new values
        for any digitally sampled data
    - [ ]   [Nearest-neighbor
        interpolation](https://en.wikipedia.org/wiki/Nearest-neighbor_interpolation)
    - [ ]   [Tricubic interpolation](https://en.wikipedia.org/wiki/Tricubic_interpolation), a
        generalization of [cubic
        interpolation](https://en.wikipedia.org/wiki/cubic_interpolation) to three
        dimensions
- [ ]   [Pareto interpolation](https://en.wikipedia.org/wiki/Pareto_interpolation): a method of
    estimating the median and other properties of a population that
    follows a [Pareto distribution](https://en.wikipedia.org/wiki/Pareto_distribution).
- [ ]   [Polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation)
    - [ ]   [Neville's algorithm](https://en.wikipedia.org/wiki/Neville's_algorithm)
- [ ]   [Spline interpolation](https://en.wikipedia.org/wiki/Spline_interpolation): Reduces
    error with [Runge's phenomenon](https://en.wikipedia.org/wiki/Runge's_phenomenon).
    - [ ]   [De Boor algorithm](https://en.wikipedia.org/wiki/De_Boor_algorithm):
        [B-splines](https://en.wikipedia.org/wiki/B-spline)
    - [ ]   [De Casteljau's algorithm](https://en.wikipedia.org/wiki/De_Casteljau's_algorithm):
        [Bézier curves](https://en.wikipedia.org/wiki/Bézier_curve)
- [ ]   [Trigonometric
    interpolation](https://en.wikipedia.org/wiki/Trigonometric_interpolation)

#### Linear algebra

- [ ]   [Eigenvalue algorithms](https://en.wikipedia.org/wiki/Eigenvalue_algorithm)
    - [ ]   [Arnoldi iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration)
    - [ ]   [Inverse iteration](https://en.wikipedia.org/wiki/Inverse_iteration)
    - [ ]   [Jacobi method](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm)
    - [ ]   [Lanczos iteration](https://en.wikipedia.org/wiki/Lanczos_iteration)
    - [ ]   [Power iteration](https://en.wikipedia.org/wiki/Power_iteration)
    - [ ]   [QR algorithm](https://en.wikipedia.org/wiki/QR_algorithm)
    - [ ]   [Rayleigh quotient
        iteration](https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration)
- [ ]   [Gram–Schmidt process](https://en.wikipedia.org/wiki/Gram–Schmidt_process):
    orthogonalizes a set of vectors
- [ ]   [Matrix multiplication
    algorithms](https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm)
    - [ ]   [Cannon's algorithm](https://en.wikipedia.org/wiki/Cannon's_algorithm): a
        [distributed algorithm](https://en.wikipedia.org/wiki/distributed_algorithm) for
        [matrix multiplication](https://en.wikipedia.org/wiki/matrix_multiplication)
        especially suitable for computers laid out in an N × N mesh
    - [ ]   [Coppersmith–Winograd
        algorithm](https://en.wikipedia.org/wiki/Coppersmith–Winograd_algorithm): square
        [matrix multiplication](https://en.wikipedia.org/wiki/matrix_multiplication)
    - [ ]   [Freivalds' algorithm](https://en.wikipedia.org/wiki/Freivalds'_algorithm): a
        randomized algorithm used to verify matrix multiplication
    - [ ]   [Strassen algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm): faster
        [matrix multiplication](https://en.wikipedia.org/wiki/matrix_multiplication)

- [ ]   Solving [systems of linear
    equations](https://en.wikipedia.org/wiki/system_of_linear_equations)
    - [ ]   [Biconjugate gradient
        method](https://en.wikipedia.org/wiki/Biconjugate_gradient_method): solves systems
        of linear equations
    - [ ]   [Conjugate gradient](https://en.wikipedia.org/wiki/Conjugate_gradient): an
        algorithm for the numerical solution of particular systems of
        linear equations
    - [ ]   [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination)
    - [ ]   [Gauss–Jordan elimination](https://en.wikipedia.org/wiki/Gauss–Jordan_elimination):
        solves systems of linear equations
    - [ ]   [Gauss–Seidel method](https://en.wikipedia.org/wiki/Gauss–Seidel_method): solves
        systems of linear equations iteratively
    - [ ]   [Levinson recursion](https://en.wikipedia.org/wiki/Levinson_recursion): solves
        equation involving a [Toeplitz
        matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix)
    - [ ]   [Stone's method](https://en.wikipedia.org/wiki/Stone's_method): also known as the
        strongly implicit procedure or SIP, is an algorithm for solving
        a sparse linear system of equations
    - [ ]   [Successive
        over-relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation) (SOR):
        method used to speed up convergence of the [Gauss–Seidel
        method](https://en.wikipedia.org/wiki/Gauss–Seidel_method)
    - [ ]   [Tridiagonal matrix
        algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) (Thomas
        algorithm): solves systems of tridiagonal equations
- [ ]   [Sparse matrix](https://en.wikipedia.org/wiki/Sparse_matrix) algorithms
    - [ ]   [Cuthill–McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm):
        reduce the [bandwidth](bandwidth_https://en.wikipedia.org/wiki/(matrix_theory)) of
        a [symmetric sparse matrix](https://en.wikipedia.org/wiki/symmetric_sparse_matrix)
    - [ ]   [Minimum degree algorithm](https://en.wikipedia.org/wiki/Minimum_degree_algorithm):
        permute the rows and columns of a [symmetric sparse
        matrix](https://en.wikipedia.org/wiki/symmetric_sparse_matrix) before applying the
        [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition)
    - [ ]   [Symbolic Cholesky
        decomposition](https://en.wikipedia.org/wiki/Symbolic_Cholesky_decomposition):
        Efficient way of storing [sparse
        matrix](https://en.wikipedia.org/wiki/sparse_matrix)

#### Monte Carlo

- [ ]   [Gibbs sampling](https://en.wikipedia.org/wiki/Gibbs_sampling): generate a sequence of
    samples from the joint probability distribution of two or more
    random variables
- [ ]   [Hybrid Monte Carlo](https://en.wikipedia.org/wiki/Hybrid_Monte_Carlo): generate a
    sequence of samples using
    [Hamiltonian](https://en.wikipedia.org/wiki/Hamiltonian_dynamics) weighted [Markov
    chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo), from a
    probability distribution which is difficult to sample directly.
- [ ]   [Metropolis–Hastings
    algorithm](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm): used to
    generate a sequence of samples from the [probability
    distribution](https://en.wikipedia.org/wiki/probability_distribution) of one or more
    variables
- [ ]   [Wang and Landau algorithm](https://en.wikipedia.org/wiki/Wang_and_Landau_algorithm):
    an extension of [Metropolis–Hastings
    algorithm](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm) sampling

#### Numerical integration

- [ ]   [MISER algorithm](https://en.wikipedia.org/wiki/MISER_algorithm): Monte Carlo
    simulation, [numerical
    integration](https://en.wikipedia.org/wiki/numerical_integration)

#### Root finding

- [ ]   [Bisection method](https://en.wikipedia.org/wiki/Bisection_method)
- [ ]   [False position method](https://en.wikipedia.org/wiki/False_position_method):
    approximates roots of a function
- [ ]   [Newton's method](https://en.wikipedia.org/wiki/Newton's_method): finds zeros of
    functions with [calculus](https://en.wikipedia.org/wiki/calculus)
- [ ]   [Halley's method](https://en.wikipedia.org/wiki/Halley's_method): uses first and second
    derivatives
- [ ]   [Secant method](https://en.wikipedia.org/wiki/Secant_method): 2-point, 1-sided
- [ ]   [False position method](https://en.wikipedia.org/wiki/False_position_method) and
    Illinois method: 2-point, bracketing
- [ ]   [Ridder's method](https://en.wikipedia.org/wiki/Ridder's_method): 3-point, exponential
    scaling
- [ ]   [Muller's method](https://en.wikipedia.org/wiki/Muller's_method): 3-point, quadratic
    interpolation

### Optimization algorithms

- [ ]   [Alpha-beta pruning](https://en.wikipedia.org/wiki/Alpha-beta_pruning): search to
    reduce number of nodes in minimax algorithm
- [ ]   [Branch and bound](https://en.wikipedia.org/wiki/Branch_and_bound)
- [ ]   [Bruss algorithm](https://en.wikipedia.org/wiki/Bruss_algorithm): see [odds
    algorithm](https://en.wikipedia.org/wiki/odds_algorithm)
- [ ]   [Chain matrix
    multiplication](https://en.wikipedia.org/wiki/Chain_matrix_multiplication)
- [ ]   [Combinatorial optimization](https://en.wikipedia.org/wiki/Combinatorial_optimization):
    optimization problems where the set of feasible solutions is
    discrete
    - [ ]   [Greedy randomized adaptive search
        procedure](https://en.wikipedia.org/wiki/Greedy_randomized_adaptive_search_procedure)
        (GRASP): successive constructions of a greedy randomized
        solution and subsequent iterative improvements of it through a
        local search
    - [ ]   [Hungarian method](https://en.wikipedia.org/wiki/Hungarian_method): a combinatorial
        optimization algorithm which solves the [assignment
        problem](https://en.wikipedia.org/wiki/assignment_problem) in polynomial time
- [ ]   [Constraint satisfaction](https://en.wikipedia.org/wiki/Constraint_satisfaction)
    - [ ]   General algorithms for the constraint satisfaction
        - [ ]   [AC-3 algorithm](https://en.wikipedia.org/wiki/AC-3_algorithm)
        - [ ]   [Difference map
            algorithm](https://en.wikipedia.org/wiki/Difference_map_algorithm)
        - [ ]   [Min conflicts
            algorithm](https://en.wikipedia.org/wiki/Min_conflicts_algorithm)
    - [ ]   [Chaff algorithm](https://en.wikipedia.org/wiki/Chaff_algorithm): an algorithm for
        solving instances of the boolean satisfiability problem
    - [ ]   [Davis–Putnam algorithm](https://en.wikipedia.org/wiki/Davis–Putnam_algorithm):
        check the validity of a first-order logic formula
    - [ ]   [Davis–Putnam–Logemann–Loveland
        algorithm](https://en.wikipedia.org/wiki/DPLL_algorithm) (DPLL): an algorithm for
        deciding the satisfiability of propositional logic formula in
        [conjunctive normal form](https://en.wikipedia.org/wiki/conjunctive_normal_form),
        i.e. for solving the [CNF-SAT](https://en.wikipedia.org/wiki/CNF-SAT) problem
    - [ ]   [Exact cover](https://en.wikipedia.org/wiki/Exact_cover) problem
        - [ ]   [Algorithm X](https://en.wikipedia.org/wiki/Algorithm_X): a [nondeterministic
            algorithm](https://en.wikipedia.org/wiki/nondeterministic_algorithm)
        - [ ]   [Dancing Links](https://en.wikipedia.org/wiki/Dancing_Links): an efficient
            implementation of Algorithm X
- [ ]   [Cross-entropy method](https://en.wikipedia.org/wiki/Cross-entropy_method): a general
    Monte Carlo approach to combinatorial and continuous multi-extremal
    optimization and [importance
    sampling](https://en.wikipedia.org/wiki/importance_sampling)
- [ ]   [Differential evolution](https://en.wikipedia.org/wiki/Differential_evolution)
- [ ]   [Dynamic Programming](https://en.wikipedia.org/wiki/Dynamic_Programming): problems
    exhibiting the properties of [overlapping
    subproblems](https://en.wikipedia.org/wiki/overlapping_subproblem) and [optimal
    substructure](https://en.wikipedia.org/wiki/optimal_substructure)
- [ ]   [Ellipsoid method](https://en.wikipedia.org/wiki/Ellipsoid_method): is an algorithm for
    solving convex optimization problems
- [ ]   [Evolutionary computation](https://en.wikipedia.org/wiki/Evolutionary_computation):
    optimization inspired by biological mechanisms of evolution
    - [ ]   [Evolution strategy](https://en.wikipedia.org/wiki/Evolution_strategy)
    - [ ]   [Gene expression
        programming](https://en.wikipedia.org/wiki/Gene_expression_programming)
    - [ ]   [Genetic algorithms](https://en.wikipedia.org/wiki/Genetic_algorithms)
        - [ ]   [Fitness proportionate
            selection](https://en.wikipedia.org/wiki/Fitness_proportionate_selection) -
            also known as roulette-wheel selection
        - [ ]   [Stochastic universal
            sampling](https://en.wikipedia.org/wiki/Stochastic_universal_sampling)
        - [ ]   [Truncation selection](https://en.wikipedia.org/wiki/Truncation_selection)
        - [ ]   [Tournament selection](https://en.wikipedia.org/wiki/Tournament_selection)
    - [ ]   [Memetic algorithm](https://en.wikipedia.org/wiki/Memetic_algorithm)
    - [ ]   [Swarm intelligence](https://en.wikipedia.org/wiki/Swarm_intelligence)
        - [ ]   [Ant colony
            optimization](https://en.wikipedia.org/wiki/Ant_colony_optimization)
        - [ ]   [Bees algorithm](https://en.wikipedia.org/wiki/Bees_algorithm): a search
            algorithm which mimics the food foraging behavior of swarms
            of honey bees
        - [ ]   [Particle swarm](https://en.wikipedia.org/wiki/Particle_swarm_optimization)
- [ ]   [golden section search](https://en.wikipedia.org/wiki/golden_section_search): an
    algorithm for finding the maximum of a real function
- [ ]   [Gradient descent](https://en.wikipedia.org/wiki/Gradient_descent)
- [ ]   [Harmony search](https://en.wikipedia.org/wiki/Harmony_search) (HS): a
    [metaheuristic](https://en.wikipedia.org/wiki/metaheuristic) algorithm mimicking the
    improvisation process of musicians
- [ ]   [Interior point method](https://en.wikipedia.org/wiki/Interior_point_method)
- [ ]   [Linear programming](https://en.wikipedia.org/wiki/Linear_programming)

    - [ ]   [Benson's algorithm](https://en.wikipedia.org/wiki/Benson's_algorithm): an
        algorithm for solving linear [vector
        optimization](https://en.wikipedia.org/wiki/vector_optimization) problems
    - [ ]   [Dantzig–Wolfe
        decomposition](https://en.wikipedia.org/wiki/Dantzig–Wolfe_decomposition): an
        algorithm for solving linear programming problems with special
        structure
    - [ ]   [Delayed column
        generation](https://en.wikipedia.org/wiki/Delayed_column_generation)
    - [ ]   [Integer linear
        programming](https://en.wikipedia.org/wiki/Integer_linear_programming): solve
        linear programming problems where some or all the unknowns are
        restricted to integer values
        - [ ]   [Branch and cut](https://en.wikipedia.org/wiki/Branch_and_cut)
        - [ ]   [Cutting-plane method](https://en.wikipedia.org/wiki/Cutting-plane_method)
    - [ ]   [Karmarkar's algorithm](https://en.wikipedia.org/wiki/Karmarkar's_algorithm): The
        first reasonably efficient algorithm that solves the [linear
        programming](https://en.wikipedia.org/wiki/linear_programming) problem in
        [polynomial time](https://en.wikipedia.org/wiki/polynomial_time).
    - [ ]   [Simplex algorithm](https://en.wikipedia.org/wiki/Simplex_algorithm): An algorithm
        for solving [linear programming](https://en.wikipedia.org/wiki/linear_programming)
        problems

- [ ]   [Line search](https://en.wikipedia.org/wiki/Line_search)
- [ ]   [Local search](Local_search_https://en.wikipedia.org/wiki/(optimization)): a
    metaheuristic for solving computationally hard optimization problems
    - [ ]   [Random-restart hill
        climbing](https://en.wikipedia.org/wiki/Random-restart_hill_climbing)
    - [ ]   [Tabu search](https://en.wikipedia.org/wiki/Tabu_search)
- [ ]   [Minimax](https://en.wikipedia.org/wiki/Minimax#Minimax_algorithm_with_alternate_moves)
    used in game programming
- [ ]   [Nearest neighbor search](https://en.wikipedia.org/wiki/Nearest_neighbor_search) (NNS):
    find closest points in a [metric space](https://en.wikipedia.org/wiki/metric_space)
    - [ ]   [Best Bin First](https://en.wikipedia.org/wiki/Best_Bin_First): find an approximate
        solution to the [Nearest neighbor
        search](https://en.wikipedia.org/wiki/Nearest_neighbor_search) problem in
        very-high-dimensional spaces
- [ ]   [Newton's method in
    optimization](https://en.wikipedia.org/wiki/Newton's_method_in_optimization)
- [ ]   [Nonlinear optimization](https://en.wikipedia.org/wiki/Nonlinear_optimization)
    - [ ]   [BFGS method](https://en.wikipedia.org/wiki/BFGS_method): A [nonlinear
        optimization](https://en.wikipedia.org/wiki/nonlinear_optimization) algorithm
    - [ ]   [Gauss–Newton algorithm](https://en.wikipedia.org/wiki/Gauss–Newton_algorithm): An
        algorithm for solving nonlinear [least
        squares](https://en.wikipedia.org/wiki/least_squares) problems.
    - [ ]   [Levenberg–Marquardt
        algorithm](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm): An
        algorithm for solving nonlinear [least
        squares](https://en.wikipedia.org/wiki/least_squares) problems.
    - [ ]   [Nelder–Mead method](https://en.wikipedia.org/wiki/Nelder–Mead_method) (downhill
        simplex method): A [nonlinear
        optimization](https://en.wikipedia.org/wiki/nonlinear_optimization) algorithm
- [ ]   [Odds algorithm](https://en.wikipedia.org/wiki/Odds_algorithm) (Bruss algorithm) :
    Finds the optimal strategy to predict a last specific event in a
    random sequence event
- [ ]   [Simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing)
- [ ]   [Stochastic tunneling](https://en.wikipedia.org/wiki/Stochastic_tunneling)
- [ ]   [Subset sum](https://en.wikipedia.org/wiki/Subset_sum_problem) algorithm

Computational science
---------------------

### Astronomy

- [ ]   [Doomsday algorithm](https://en.wikipedia.org/wiki/Doomsday_algorithm): day of the week
- [ ]   [Zeller's congruence](https://en.wikipedia.org/wiki/Zeller's_congruence) is an
    algorithm to calculate the day of the week for any Julian or
    Gregorian calendar date
- [ ]   various [Easter algorithms](https://en.wikipedia.org/wiki/Computus) are used to
    calculate the day of Easter

### Bioinformatics

- [ ]   [Basic Local Alignment Search
    Tool](https://en.wikipedia.org/wiki/Basic_Local_Alignment_Search_Tool) also known as
    BLAST: an algorithm for comparing primary biological sequence
    information
- [ ]   [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm): calculate the
    optimal alignment of two sets of points in order to compute the
    [root mean squared deviation](https://en.wikipedia.org/wiki/RMSD) between two protein
    structures.
- [ ]   [Velvet](Velvet_https://en.wikipedia.org/wiki/(algorithm)): a set of algorithms
    manipulating [de Bruijn graphs](https://en.wikipedia.org/wiki/de_Bruijn_graph) for
    genomic [sequence assembly](https://en.wikipedia.org/wiki/sequence_assembly)
- [ ]   [Sorting by signed
    reversals](https://en.wikipedia.org/wiki/Sorting_by_signed_reversals): an algorithm for
    understanding genomic evolution.
- [ ]   [Maximum parsimony
    (phylogenetics)](Maximum_parsimony_https://en.wikipedia.org/wiki/(phylogenetics)): an
    algorithm for finding the simplest phylogenetic tree to explain a
    given character matrix.
- [ ]   [UPGMA](https://en.wikipedia.org/wiki/UPGMA): a distance-based phylogenetic tree
    construction algorithm.

### Geoscience

- [ ]   [Vincenty's formulae](https://en.wikipedia.org/wiki/Vincenty's_formulae): a fast
    algorithm to calculate the distance between two latitude/longitude
    points on an ellipsoid
- [ ]   [Geohash](https://en.wikipedia.org/wiki/Geohash): a public domain algorithm that
    encodes a decimal latitude/longitude pair as a hash string

### Linguistics

- [ ]   [Lesk algorithm](https://en.wikipedia.org/wiki/Lesk_algorithm): word sense
    disambiguation
- [ ]   [Stemming algorithm](https://en.wikipedia.org/wiki/Stemming): a method of reducing
    words to their stem, base, or root form
- [ ]   [Sukhotin's algorithm](https://en.wikipedia.org/wiki/Sukhotin's_algorithm): a
    statistical classification algorithm for classifying characters in a
    text as vowels or consonants

### Medicine

- [ ]   [ESC algorithm](https://en.wikipedia.org/wiki/ESC_algorithm) for the diagnosis of heart
    failure
- [ ]   [Manning Criteria](https://en.wikipedia.org/wiki/Manning_Criteria) for irritable bowel
    syndrome
- [ ]   [Pulmonary embolism](https://en.wikipedia.org/wiki/Pulmonary_embolism#Algorithms)
    diagnostic algorithms
- [ ]   [Texas Medication Algorithm
    Project](https://en.wikipedia.org/wiki/Texas_Medication_Algorithm_Project)

### Physics

- [ ]   [Constraint algorithm](https://en.wikipedia.org/wiki/Constraint_algorithm): a class of
    algorithms for satisfying constraints for bodies that obey Newton's
    equations of motion
- [ ]   [Demon algorithm](https://en.wikipedia.org/wiki/Demon_algorithm): a [Monte Carlo
    method](https://en.wikipedia.org/wiki/Monte_Carlo_method) for efficiently sampling
    members of a [microcanonical
    ensemble](https://en.wikipedia.org/wiki/microcanonical_ensemble) with a given energy
- [ ]   [Featherstone's algorithm](https://en.wikipedia.org/wiki/Featherstone's_algorithm):
    compute the effects of forces applied to a structure of joints and
    links
- [ ]   [Ground state](https://en.wikipedia.org/wiki/Ground_state) approximation
    - [ ]   [Variational method](https://en.wikipedia.org/wiki/Variational_method)
        - [ ]   [Ritz method](https://en.wikipedia.org/wiki/Ritz_method)
- [ ]   [N-body problems](https://en.wikipedia.org/wiki/N-body_problem)
    - [ ]   [Barnes–Hut simulation](https://en.wikipedia.org/wiki/Barnes–Hut_simulation):
        Solves the n-body problem in an approximate way that has the
        order instead of as in a direct-sum simulation.
    - [ ]   [Fast multipole method](https://en.wikipedia.org/wiki/Fast_multipole_method) (FMM):
        speeds up the calculation of long-ranged forces
- [ ]   [Rainflow-counting
    algorithm](https://en.wikipedia.org/wiki/Rainflow-counting_algorithm): Reduces a
    complex [stress](stress_https://en.wikipedia.org/wiki/(physics)) history to a count of
    elementary stress-reversals for use in
    [fatigue](fatigue_https://en.wikipedia.org/wiki/(material)) analysis
- [ ]   [Sweep and prune](https://en.wikipedia.org/wiki/Sweep_and_prune): a broad phase
    algorithm used during [collision
    detection](https://en.wikipedia.org/wiki/collision_detection) to limit the number of
    pairs of solids that need to be checked for collision
- [ ]   [VEGAS algorithm](https://en.wikipedia.org/wiki/VEGAS_algorithm): a method for reducing
    error in [Monte Carlo
    simulations](https://en.wikipedia.org/wiki/Monte_Carlo_simulation)

### Statistics

- [ ]   [Algorithms for calculating
    variance](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance): avoiding
    instability and numerical overflow
- [ ]   [Approximate counting
    algorithm](https://en.wikipedia.org/wiki/Approximate_counting_algorithm): Allows
    counting large number of events in a small register
- [ ]   [Bayesian statistics](https://en.wikipedia.org/wiki/Bayesian_statistics)
    - [ ]   [Nested sampling
        algorithm](https://en.wikipedia.org/wiki/Nested_sampling_algorithm): a
        computational approach to the problem of comparing models in
        Bayesian statistics
- [ ]   [Clustering Algorithms](https://en.wikipedia.org/wiki/Data_clustering)
    - [ ]   [Average-linkage clustering](https://en.wikipedia.org/wiki/UPGMA): a simple
        agglomerative clustering algorithm
    - [ ]   [Canopy clustering
        algorithm](https://en.wikipedia.org/wiki/Canopy_clustering_algorithm): an
        unsupervised pre-clustering algorithm related to the K-means
        algorithm
    - [ ]   [Complete-linkage
        clustering](https://en.wikipedia.org/wiki/Complete-linkage_clustering): a simple
        agglomerative clustering algorithm
    - [ ]   [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN): a density based clustering
        algorithm
    - [ ]   [Expectation-maximization
        algorithm](https://en.wikipedia.org/wiki/Expectation-maximization_algorithm)
    - [ ]   [Fuzzy clustering](https://en.wikipedia.org/wiki/Fuzzy_clustering): a class of
        clustering algorithms where each point has a degree of belonging
        to clusters
        - [ ]   [Fuzzy
            c-means](https://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_c-means_clustering)
        - [ ]   [FLAME clustering](https://en.wikipedia.org/wiki/FLAME_clustering) (Fuzzy
            clustering by Local Approximation of MEmberships): define
            clusters in the dense parts of a dataset and perform cluster
            assignment solely based on the neighborhood relationships
            among objects
    - [ ]   [KHOPCA clustering
        algorithm](https://en.wikipedia.org/wiki/KHOPCA_clustering_algorithm): a local
        clustering algorithm, which produces hierarchical multi-hop
        clusters in static and mobile environments.
    - [ ]   [k-means clustering](https://en.wikipedia.org/wiki/k-means_clustering): cluster
        objects based on attributes into partitions
    - [ ]   [k-means++](https://en.wikipedia.org/wiki/k-means++): a variation of this, using
        modified random seeds
    - [ ]   [k-medoids](https://en.wikipedia.org/wiki/k-medoids): similar to k-means, but
        chooses datapoints or [medoids](https://en.wikipedia.org/wiki/medoid) as centers
    - [ ]   [Linde–Buzo–Gray
        algorithm](https://en.wikipedia.org/wiki/Linde–Buzo–Gray_algorithm): a vector
        quantization algorithm to derive a good codebook
    - [ ]   [Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd's_algorithm) (Voronoi
        iteration or relaxation): group data points into a given number
        of categories, a popular algorithm for [k-means
        clustering](https://en.wikipedia.org/wiki/k-means_clustering)
    - [ ]   [OPTICS](https://en.wikipedia.org/wiki/OPTICS_algorithm): a density based
        clustering algorithm with a visual evaluation method
    - [ ]   [Single-linkage
        clustering](https://en.wikipedia.org/wiki/Single-linkage_clustering): a simple
        agglomerative clustering algorithm
    - [ ]   [SUBCLU](https://en.wikipedia.org/wiki/SUBCLU): a subspace clustering algorithm
    - [ ]   [Ward's method](https://en.wikipedia.org/wiki/Ward's_method) : an agglomerative
        clustering algorithm, extended to more general Lance–Williams
        algorithms
    - [ ]   [WACA clustering
        algorithm](https://en.wikipedia.org/wiki/WACA_clustering_algorithm): a local
        clustering algorithm with potentially multi-hop structures; for
        dynamic networks
- [ ]   [Estimation Theory](https://en.wikipedia.org/wiki/Estimation_theory)
    - [ ]   [Expectation-maximization
        algorithm](https://en.wikipedia.org/wiki/Expectation-maximization_algorithm) A
        class of related algorithms for finding maximum likelihood
        estimates of parameters in probabilistic models
        - [ ]   [Ordered subset expectation
            maximization](https://en.wikipedia.org/wiki/Ordered_subset_expectation_maximization)
            (OSEM): used in [medical
            imaging](https://en.wikipedia.org/wiki/medical_imaging) for [positron emission
            tomography](https://en.wikipedia.org/wiki/positron_emission_tomography),
            [single photon emission computed
            tomography](https://en.wikipedia.org/wiki/single_photon_emission_computed_tomography)
            and [X-ray](https://en.wikipedia.org/wiki/X-ray) computed tomography.
    - [ ]   [Odds algorithm](https://en.wikipedia.org/wiki/Odds_algorithm) (Bruss algorithm)
        Optimal online search for distinguished value in sequential
        random input
    - [ ]   [Kalman filter](https://en.wikipedia.org/wiki/Kalman_filter): estimate the state of
        a linear [dynamic system](https://en.wikipedia.org/wiki/Dynamical_system) from a
        series of noisy measurements
- [ ]   [False nearest neighbor
    algorithm](https://en.wikipedia.org/wiki/False_nearest_neighbor_algorithm) (FNN)
    estimates [fractal dimension](https://en.wikipedia.org/wiki/fractal_dimension)
- [ ]   [Hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model)
    - [ ]   [Baum–Welch algorithm](https://en.wikipedia.org/wiki/Baum–Welch_algorithm): compute
        maximum likelihood estimates and [posterior
        mode](https://en.wikipedia.org/wiki/Maximum_a_posteriori) estimates for the
        parameters of a [hidden markov
        model](https://en.wikipedia.org/wiki/hidden_markov_model)
    - [ ]   [Forward-backward
        algorithm](https://en.wikipedia.org/wiki/Forward-backward_algorithm) a dynamic
        programming algorithm for computing the probability of a
        particular observation sequence
    - [ ]   [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm): find the most
        likely sequence of hidden states in a [hidden markov
        model](https://en.wikipedia.org/wiki/hidden_markov_model)
- [ ]   [Partial least squares
    regression](https://en.wikipedia.org/wiki/Partial_least_squares_regression): finds a
    linear model describing some predicted variables in terms of other
    observable variables
- [ ]   [Queuing theory](https://en.wikipedia.org/wiki/Queuing_theory)
    - [ ]   [Buzen's algorithm](https://en.wikipedia.org/wiki/Buzen's_algorithm): an algorithm
        for calculating the normalization constant G(K) in the
        [Gordon–Newell theorem](https://en.wikipedia.org/wiki/Gordon–Newell_theorem)
- [ ]   [RANSAC](https://en.wikipedia.org/wiki/RANSAC) (an abbreviation for “RANdom SAmple
    Consensus”): an iterative method to estimate parameters of a
    mathematical model from a set of observed data which contains
    outliers
- [ ]   [Scoring algorithm](https://en.wikipedia.org/wiki/Scoring_algorithm): is a form of
    [Newton's method](https://en.wikipedia.org/wiki/Newton's_method) used to solve [maximum
    likelihood](https://en.wikipedia.org/wiki/maximum_likelihood) equations numerically
- [ ]   [Yamartino method](https://en.wikipedia.org/wiki/Yamartino_method): calculate an
    approximation to the standard deviation σθ of wind direction θ
    during a single pass through the incoming data
- [ ]   [Ziggurat algorithm](https://en.wikipedia.org/wiki/Ziggurat_algorithm): generate random
    numbers from a non-uniform distribution

Computer science
----------------

### Computer architecture

- [ ]   [Tomasulo algorithm](https://en.wikipedia.org/wiki/Tomasulo_algorithm): allows
    sequential instructions that would normally be stalled due to
    certain dependencies to execute non-sequentially

### Computer graphics

- [ ]   [Clipping](Clipping_https://en.wikipedia.org/wiki/(computer_graphics))
    - [ ]   [Line clipping](https://en.wikipedia.org/wiki/Line_clipping)
        - [ ]   [Cohen–Sutherland](https://en.wikipedia.org/wiki/Cohen–Sutherland)
        - [ ]   [Cyrus–Beck](https://en.wikipedia.org/wiki/Cyrus–Beck)
        - [ ]   [Fast-clipping](https://en.wikipedia.org/wiki/Fast_clipping)
        - [ ]   [Liang–Barsky](https://en.wikipedia.org/wiki/Liang–Barsky)
        - [ ]   [Nicholl–Lee–Nicholl](https://en.wikipedia.org/wiki/Nicholl–Lee–Nicholl)
    - [ ]   Polygon clipping
        - [ ]   [Sutherland–Hodgman](https://en.wikipedia.org/wiki/Sutherland–Hodgman)
        - [ ]   [Vatti](https://en.wikipedia.org/wiki/Vatti_clipping_algorithm)
        - [ ]   [Weiler–Atherton](https://en.wikipedia.org/wiki/Weiler–Atherton)
- [ ]   [Contour lines](https://en.wikipedia.org/wiki/Contour_line) and
    [Isosurfaces](https://en.wikipedia.org/wiki/Isosurface)
    - [ ]   [Marching cubes](https://en.wikipedia.org/wiki/Marching_cubes): extract a polygonal
        mesh of an isosurface from a three-dimensional scalar field
        (sometimes called voxels)
    - [ ]   [Marching squares](https://en.wikipedia.org/wiki/Marching_squares): generate
        contour lines for a two-dimensional scalar field
    - [ ]   [Marching tetrahedrons](https://en.wikipedia.org/wiki/Marching_tetrahedrons): an
        alternative to [Marching cubes](https://en.wikipedia.org/wiki/Marching_cubes)
- [ ]   [Discrete Green's Theorem](https://en.wikipedia.org/wiki/Discrete_Green's_Theorem): is
    an algorithm for computing double integral over a generalized
    rectangular domain in constant time. It is a natural extension to
    the summed area table algorithm
- [ ]   [Flood fill](https://en.wikipedia.org/wiki/Flood_fill): fills a connected region of a
    multi-dimensional array with a specified symbol
- [ ]   [Global illumination](https://en.wikipedia.org/wiki/Global_illumination) algorithms:
    Considers direct illumination and reflection from other objects.
    - [ ]   [Ambient occlusion](https://en.wikipedia.org/wiki/Ambient_occlusion)
    - [ ]   [Beam tracing](https://en.wikipedia.org/wiki/Beam_tracing)
    - [ ]   [Cone tracing](https://en.wikipedia.org/wiki/Cone_tracing)
    - [ ]   [Image-based lighting](https://en.wikipedia.org/wiki/Image-based_lighting)
    - [ ]   [Metropolis light
        transport](https://en.wikipedia.org/wiki/Metropolis_light_transport)
    - [ ]   [Path tracing](https://en.wikipedia.org/wiki/Path_tracing)
    - [ ]   [Photon mapping](https://en.wikipedia.org/wiki/Photon_mapping)
    - [ ]   [Radiosity](Radiosity_https://en.wikipedia.org/wiki/(3D_computer_graphics))
    - [ ]   [Ray tracing](Ray_tracing_https://en.wikipedia.org/wiki/(graphics))
- [ ]   [Hidden surface removal](https://en.wikipedia.org/wiki/Hidden_surface_determination) or
    [Visual surface
    determination](https://en.wikipedia.org/wiki/Hidden_surface_determination)
    - [ ]   [Newell's algorithm](https://en.wikipedia.org/wiki/Newell's_algorithm): eliminate
        polygon cycles in the depth sorting required in hidden surface
        removal
    - [ ]   [Painter's algorithm](https://en.wikipedia.org/wiki/Painter's_algorithm): detects
        visible parts of a 3-dimensional scenery
    - [ ]   [Scanline rendering](https://en.wikipedia.org/wiki/Scanline_rendering): constructs
        an image by moving an imaginary line over the image
    - [ ]   [Warnock algorithm](https://en.wikipedia.org/wiki/Warnock_algorithm)
- [ ]   [Line Drawing](https://en.wikipedia.org/wiki/Line_drawing_algorithm): graphical
    algorithm for approximating a line segment on discrete graphical
    media.
    - [ ]   [Bresenham's line
        algorithm](https://en.wikipedia.org/wiki/Bresenham's_line_algorithm): plots points
        of a 2-dimensional array to form a straight line between 2
        specified points (uses decision variables)
    - [ ]   [DDA line
        algorithm](Digital_Differential_Analyzer_https://en.wikipedia.org/wiki/(graphics_algorithm)):
        plots points of a 2-dimensional array to form a straight line
        between 2 specified points (uses floating-point math)
    - [ ]   [Xiaolin Wu's line
        algorithm](https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm): algorithm
        for line antialiasing.
- [ ]   [Midpoint circle algorithm](https://en.wikipedia.org/wiki/Midpoint_circle_algorithm):
    an algorithm used to determine the points needed for drawing a
    circle
- [ ]   [Ramer–Douglas–Peucker
    algorithm](https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm): Given a
    'curve' composed of line segments to find a curve not too dissimilar
    but that has fewer points
- [ ]   [Shading](https://en.wikipedia.org/wiki/Shading)
    - [ ]   [Gouraud shading](https://en.wikipedia.org/wiki/Gouraud_shading): an algorithm to
        simulate the differing effects of light and colour across the
        surface of an object in 3D computer graphics
    - [ ]   [Phong shading](https://en.wikipedia.org/wiki/Phong_shading): an algorithm to
        interpolate surface normal-vectors for surface shading in 3D
        computer graphics
- [ ]   [Slerp](https://en.wikipedia.org/wiki/Slerp) (spherical linear interpolation):
    quaternion interpolation for the purpose of animating 3D rotation
- [ ]   [Summed area table](https://en.wikipedia.org/wiki/Summed_area_table) (also known as an
    integral image): an algorithm for computing the sum of values in a
    rectangular subset of a grid in constant time

### Cryptography

- [ ]   [Asymmetric (public key)
    encryption](https://en.wikipedia.org/wiki/Asymmetric_key_algorithm):
    - [ ]   [ElGamal](https://en.wikipedia.org/wiki/ElGamal_encryption)
    - [ ]   [Elliptic curve
        cryptography](https://en.wikipedia.org/wiki/Elliptic_curve_cryptography)
    - [ ]   [MAE1](https://en.wikipedia.org/wiki/Matei_Array_Encreption_1)
    - [ ]   [NTRUEncrypt](https://en.wikipedia.org/wiki/NTRUEncrypt)
    - [ ]   [RSA](RSA_https://en.wikipedia.org/wiki/(cryptosystem))
- [ ]   [Digital signatures](https://en.wikipedia.org/wiki/Digital_signature) (asymmetric
    authentication):
    - [ ]   [DSA](https://en.wikipedia.org/wiki/Digital_Signature_Algorithm), and its variants:
        - [ ]   [ECDSA](https://en.wikipedia.org/wiki/ECDSA) and [Deterministic
            ECDSA](https://tools.ietf.org/html/rfc6979)
        - [ ]   [EdDSA](https://en.wikipedia.org/wiki/EdDSA) (Ed25519)
    - [ ]   [RSA](RSA_https://en.wikipedia.org/wiki/(cryptosystem))
- [ ]   [Cryptographic hash
    functions](https://en.wikipedia.org/wiki/Cryptographic_hash_function) (see also the
    section on message authentication codes):
    - [ ]   [BLAKE](BLAKE_https://en.wikipedia.org/wiki/(hash_function))
    - [ ]   [MD5](https://en.wikipedia.org/wiki/MD5) – Note that there is now a method of
        generating collisions for MD5
    - [ ]   [RIPEMD-160](https://en.wikipedia.org/wiki/RIPEMD-160)
    - [ ]   [SHA-1](https://en.wikipedia.org/wiki/SHA-1) – Note that there is now a method of
        generating collisions for SHA-1
    - [ ]   [SHA-2](https://en.wikipedia.org/wiki/SHA-2) (SHA-224, SHA-256, SHA-384, SHA-512)
    - [ ]   [SHA-3](https://en.wikipedia.org/wiki/SHA-3) (SHA3-224, SHA3-256, SHA3-384,
        SHA3-512, SHAKE128, SHAKE256)
    - [ ]   [Tiger](Tiger_https://en.wikipedia.org/wiki/(hash)) (TTH), usually used in [Tiger
        tree hashes](Hash_tree_https://en.wikipedia.org/wiki/(persistent_data_structure))
    - [ ]   [WHIRLPOOL](https://en.wikipedia.org/wiki/WHIRLPOOL)
- [ ]   [Cryptographically secure pseudo-random number
    generators](https://en.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator)
    - [ ]   [Blum Blum Shub](https://en.wikipedia.org/wiki/Blum_Blum_Shub) - based on the
        hardness of [factorization](https://en.wikipedia.org/wiki/integer_factorization)
    - [ ]   [Fortuna](Fortuna_https://en.wikipedia.org/wiki/(PRNG)), intended as an improvement
        on [Yarrow algorithm](https://en.wikipedia.org/wiki/Yarrow_algorithm)
    - [ ]   [Linear-feedback shift
        register](https://en.wikipedia.org/wiki/Linear-feedback_shift_register) (note: many
        LFSR-based algorithms are weak or have been broken)
    - [ ]   [Yarrow algorithm](https://en.wikipedia.org/wiki/Yarrow_algorithm)
- [ ]   [Key exchange](https://en.wikipedia.org/wiki/Key_exchange)
    - [ ]   [Diffie–Hellman key
        exchange](https://en.wikipedia.org/wiki/Diffie–Hellman_key_exchange)
    - [ ]   [Elliptic-curve
        Diffie-Hellman](https://en.wikipedia.org/wiki/Elliptic-curve_Diffie-Hellman) (ECDH)
- [ ]   [Key derivation functions](https://en.wikipedia.org/wiki/Key_derivation_function),
    often used for [password hashing](https://en.wikipedia.org/wiki/password_hashing) and
    [key stretching](https://en.wikipedia.org/wiki/key_stretching)
    - [ ]   [bcrypt](https://en.wikipedia.org/wiki/bcrypt)
    - [ ]   [PBKDF2](https://en.wikipedia.org/wiki/PBKDF2)
    - [ ]   [scrypt](https://en.wikipedia.org/wiki/scrypt)
    - [ ]   [Argon2](https://en.wikipedia.org/wiki/Argon2)
- [ ]   [Message authentication
    codes](https://en.wikipedia.org/wiki/Message_authentication_code) (symmetric
    authentication algorithms, which take a key as a parameter):
    - [ ]   [HMAC](https://en.wikipedia.org/wiki/keyed-hash_message_authentication_code):
        keyed-hash message authentication
    - [ ]   [Poly1305](https://en.wikipedia.org/wiki/Poly1305)
    - [ ]   [SipHash](https://en.wikipedia.org/wiki/SipHash)
- [ ]   [Secret sharing](https://en.wikipedia.org/wiki/Secret_sharing), Secret Splitting, Key
    Splitting, M of N algorithms
    - [ ]   Blakey's Scheme
    - [ ]   [Shamir's Scheme](https://en.wikipedia.org/wiki/Shamir's_Secret_Sharing)
- [ ]   [Symmetric (secret key)
    encryption](https://en.wikipedia.org/wiki/symmetric_key_algorithm):
    - [ ]   [Advanced Encryption
        Standard](https://en.wikipedia.org/wiki/Advanced_Encryption_Standard) (AES), winner
        of [NIST](https://en.wikipedia.org/wiki/NIST) competition, also known as
        [Rijndael](https://en.wikipedia.org/wiki/Rijndael)
    - [ ]   [Blowfish](Blowfish_https://en.wikipedia.org/wiki/(cipher))
    - [ ]   [Twofish](https://en.wikipedia.org/wiki/Twofish)
    - [ ]   [Threefish](https://en.wikipedia.org/wiki/Threefish)
    - [ ]   [Data Encryption Standard](https://en.wikipedia.org/wiki/Data_Encryption_Standard)
        (DES), sometimes DE Algorithm, winner of NBS selection
        competition, replaced by AES for most purposes
    - [ ]   [IDEA](https://en.wikipedia.org/wiki/International_Data_Encryption_Algorithm)
    - [ ]   [RC4 (cipher)](RC4_https://en.wikipedia.org/wiki/(cipher))
    - [ ]   [Tiny Encryption
        Algorithm](https://en.wikipedia.org/wiki/Tiny_Encryption_Algorithm) (TEA)
    - [ ]   [Salsa20](https://en.wikipedia.org/wiki/Salsa20), and its updated variant
        [ChaCha20](https://en.wikipedia.org/wiki/Salsa20#ChaCha_variant)
- [ ]   [Post-quantum cryptography](https://en.wikipedia.org/wiki/Post-quantum_cryptography)
- [ ]   [Proof-of-work algorithms](https://en.wikipedia.org/wiki/Proof-of-work_system)

### Digital logic

- [ ]   Boolean minimization
    - [ ]   [Quine–McCluskey
        algorithm](https://en.wikipedia.org/wiki/Quine–McCluskey_algorithm): Also called as
        Q-M algorithm, programmable method for simplifying the boolean
        equations.
    - [ ]   [Petrick's method](https://en.wikipedia.org/wiki/Petrick's_method): Another
        algorithm for boolean simplification.
    - [ ]   [Espresso heuristic logic
        minimizer](https://en.wikipedia.org/wiki/Espresso_heuristic_logic_minimizer): Fast
        algorithm for boolean function minimization.

### Machine learning and statistical classification

- [ ]   [ALOPEX](https://en.wikipedia.org/wiki/ALOPEX): a correlation-based [machine-learning
    algorithm](https://en.wikipedia.org/wiki/Machine_learning)
- [ ]   [Association rule learning](https://en.wikipedia.org/wiki/Association_rule_learning):
    discover interesting relations between variables, used in [data
    mining](https://en.wikipedia.org/wiki/data_mining)
    - [ ]   [Apriori algorithm](https://en.wikipedia.org/wiki/Apriori_algorithm)
    - [ ]   [Eclat algorithm](https://en.wikipedia.org/wiki/Eclat_algorithm)
    - [ ]   [FP-growth
        algorithm](https://en.wikipedia.org/wiki/Association_rule_learning#FP-growth_algorithm)
    - [ ]   [One-attribute rule](https://en.wikipedia.org/wiki/One-attribute_rule)
    - [ ]   [Zero-attribute
        rule](https://en.wikipedia.org/wiki/Association_rule_learning#Zero-attribute_rule)
- [ ]   [Boosting (meta-algorithm)](Boosting_https://en.wikipedia.org/wiki/(meta-algorithm)):
    Use many weak learners to boost effectiveness
    - [ ]   [AdaBoost](https://en.wikipedia.org/wiki/AdaBoost): adaptive boosting
    - [ ]   [BrownBoost](https://en.wikipedia.org/wiki/BrownBoost):a boosting algorithm that
        may be robust to noisy datasets
    - [ ]   [LogitBoost](https://en.wikipedia.org/wiki/LogitBoost): [logistic
        regression](https://en.wikipedia.org/wiki/logistic_regression) boosting
    - [ ]   [LPBoost](https://en.wikipedia.org/wiki/LPBoost): [linear
        programming](https://en.wikipedia.org/wiki/linear_programming) boosting
- [ ]   [Bootstrap aggregating](https://en.wikipedia.org/wiki/Bootstrap_aggregating) (bagging):
    technique to improve stability and classification accuracy
- [ ]   [Computer Vision](https://en.wikipedia.org/wiki/Computer_Vision)
    - [ ]   [Grabcut](https://en.wikipedia.org/wiki/Grabcut) based on [Graph
        cuts](https://en.wikipedia.org/wiki/Graph_cuts_in_computer_vision)
- [ ]   [Decision Trees](https://en.wikipedia.org/wiki/Decision_tree_learning)
    - [ ]   [C4.5 algorithm](https://en.wikipedia.org/wiki/C4.5_algorithm): an extension to ID3
    - [ ]   [ID3 algorithm](https://en.wikipedia.org/wiki/ID3_algorithm) (Iterative
        Dichotomiser 3): Use heuristic to generate small decision trees
- [ ]   [Clustering](https://en.wikipedia.org/wiki/Cluster_analysis): Class of [unsupervised
    learning](https://en.wikipedia.org/wiki/unsupervised_learning) algorithms for grouping
    and bucketing related input vector.
    - [ ]   [k-nearest neighbors](https://en.wikipedia.org/wiki/k-nearest_neighbors) (k-NN): a
        method for classifying objects based on closest training
        examples in the [feature space](https://en.wikipedia.org/wiki/feature_space)
- [ ]   [Linde–Buzo–Gray algorithm](https://en.wikipedia.org/wiki/Linde–Buzo–Gray_algorithm): a
    vector quantization algorithm used to derive a good codebook
- [ ]   [Locality-sensitive hashing](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)
    (LSH): a method of performing probabilistic dimension reduction of
    high-dimensional data
- [ ]   [Neural Network](https://en.wikipedia.org/wiki/Artificial_neural_network)
    - [ ]   [Backpropagation](https://en.wikipedia.org/wiki/Backpropagation): A [supervised
        learning](https://en.wikipedia.org/wiki/supervised_learning) method which requires
        a teacher that knows, or can calculate, the desired output for
        any given input
    - [ ]   [Hopfield net](https://en.wikipedia.org/wiki/Hopfield_net): a [Recurrent neural
        network](https://en.wikipedia.org/wiki/Recurrent_neural_network) in which all
        connections are symmetric
    - [ ]   [Perceptron](https://en.wikipedia.org/wiki/Perceptron): the simplest kind of
        feedforward neural network: a [linear
        classifier](https://en.wikipedia.org/wiki/linear_classifier).
    - [ ]   [Pulse-coupled neural
        networks](https://en.wikipedia.org/wiki/Pulse-coupled_neural_networks) (PCNN):
        [Neural models](https://en.wikipedia.org/wiki/Artificial_neural_network) proposed
        by modeling a cat's [visual cortex](https://en.wikipedia.org/wiki/visual_cortex)
        and developed for high-performance
        [biomimetic](https://en.wikipedia.org/wiki/Bionics) image processing.
    - [ ]   [Radial basis function
        network](https://en.wikipedia.org/wiki/Radial_basis_function_network): an
        artificial neural network that uses radial [basis
        functions](https://en.wikipedia.org/wiki/basis_function) as activation functions
    - [ ]   [Self-organizing map](https://en.wikipedia.org/wiki/Self-organizing_map): an
        unsupervised network that produces a low-dimensional
        representation of the input space of the training samples
- [ ]   [Random forest](https://en.wikipedia.org/wiki/Random_forest): classify using many
    decision trees
- [ ]   [Reinforcement Learning](https://en.wikipedia.org/wiki/Reinforcement_Learning):
    - [ ]   [Q-learning](https://en.wikipedia.org/wiki/Q-learning): learn an action-value
        function that gives the expected utility of taking a given
        action in a given state and following a fixed policy thereafter
    - [ ]   [State-Action-Reward-State-Action](https://en.wikipedia.org/wiki/State-Action-Reward-State-Action)
        (SARSA): learn a [Markov decision
        process](https://en.wikipedia.org/wiki/Markov_decision_process) policy
    - [ ]   [Temporal difference
        learning](https://en.wikipedia.org/wiki/Temporal_difference_learning)
- [ ]   [Relevance Vector Machine](https://en.wikipedia.org/wiki/Relevance_Vector_Machine)
    (RVM): similar to SVM, but provides probabilistic classification
- [ ]   [Supervised Learning](https://en.wikipedia.org/wiki/Supervised_Learning): Learning by
    examples (labelled data-set split into training-set and test-set)
- [ ]   [Support Vector Machines](https://en.wikipedia.org/wiki/Support_Vector_Machines) (SVM):
    a set of methods which divide multidimensional data by finding a
    dividing hyperplane with the maximum margin between the two sets
    - [ ]   [Structured SVM](https://en.wikipedia.org/wiki/Structured_SVM): allows training of
        a classifier for general structured output labels.
- [ ]   [Winnow algorithm](https://en.wikipedia.org/wiki/Winnow_algorithm): related to the
    perceptron, but uses a [multiplicative weight-update
    scheme](https://en.wikipedia.org/wiki/Multiplicative_Weight_Update_Method)

### Programming language theory

- [ ]   [C3 linearization](https://en.wikipedia.org/wiki/C3_linearization): an algorithm used
    primarily to obtain a consistent linearization of a multiple
    inheritance hierarchy in object-oriented programming
- [ ]   [Chaitin's algorithm](https://en.wikipedia.org/wiki/Chaitin's_algorithm): a bottom-up,
    graph coloring register allocation algorithm that uses cost/degree
    as its spill metric
- [ ]   [Hindley–Milner type inference
    algorithm](https://en.wikipedia.org/wiki/Hindley-Milner_type_inference)
- [ ]   [Rete algorithm](https://en.wikipedia.org/wiki/Rete_algorithm): an efficient pattern
    matching algorithm for implementing [production
    rule](Start_symbol_https://en.wikipedia.org/wiki/(formal_languages)) systems
- [ ]   [Sethi-Ullman algorithm](https://en.wikipedia.org/wiki/Sethi-Ullman_algorithm):
    generate optimal code for arithmetic expressions

#### Parsing

- [ ]   [CYK algorithm](https://en.wikipedia.org/wiki/CYK_algorithm): An O(n^3^) algorithm for
    parsing [context-free grammars](https://en.wikipedia.org/wiki/context-free_grammar) in
    [Chomsky normal form](https://en.wikipedia.org/wiki/Chomsky_normal_form)
- [ ]   [Earley parser](https://en.wikipedia.org/wiki/Earley_parser): Another O(n^3^) algorithm
    for parsing any [context-free
    grammar](https://en.wikipedia.org/wiki/context-free_grammar)
- [ ]   [GLR parser](https://en.wikipedia.org/wiki/GLR_parser):An algorithm for parsing any
    [context-free grammar](https://en.wikipedia.org/wiki/context-free_grammar) by [Masaru
    Tomita](https://en.wikipedia.org/wiki/Masaru_Tomita). It is tuned for deterministic
    grammars, on which it performs almost [linear
    time](https://en.wikipedia.org/wiki/linear_time) and O(n^3^) in worst case.
- [ ]   [Inside-outside algorithm](https://en.wikipedia.org/wiki/Inside-outside_algorithm): An
    O(n^3^) algorithm for re-estimating production probabilities in
    [probabilistic context-free
    grammars](https://en.wikipedia.org/wiki/probabilistic_context-free_grammar)
- [ ]   [LL parser](https://en.wikipedia.org/wiki/LL_parser): A relatively simple [linear
    time](https://en.wikipedia.org/wiki/linear_time) parsing algorithm for a limited class
    of [context-free grammars](https://en.wikipedia.org/wiki/context-free_grammar)
- [ ]   [LR parser](https://en.wikipedia.org/wiki/LR_parser): A more complex [linear
    time](https://en.wikipedia.org/wiki/linear_time) parsing algorithm for a larger class
    of [context-free grammars](https://en.wikipedia.org/wiki/context-free_grammar).
    Variants:
    - [ ]   [Canonical LR parser](https://en.wikipedia.org/wiki/Canonical_LR_parser)
    - [ ]   [LALR (Look-ahead LR) parser](https://en.wikipedia.org/wiki/Look-ahead_LR_parser)
    - [ ]   [Operator-precedence
        parser](https://en.wikipedia.org/wiki/Operator-precedence_parser)
    - [ ]   [SLR (Simple LR) parser](https://en.wikipedia.org/wiki/Simple_LR_parser)
    - [ ]   [Simple precedence parser](https://en.wikipedia.org/wiki/Simple_precedence_parser)
- [ ]   [Packrat parser](https://en.wikipedia.org/wiki/Packrat_parser): A [linear
    time](https://en.wikipedia.org/wiki/linear_time) parsing algorithm supporting some
    [context-free grammars](https://en.wikipedia.org/wiki/context-free_grammar) and
    [parsing expression grammars](https://en.wikipedia.org/wiki/parsing_expression_grammar)
- [ ]   [Recursive descent parser](https://en.wikipedia.org/wiki/Recursive_descent_parser): A
    [top-down parser](https://en.wikipedia.org/wiki/top-down_parsing) suitable for LL(*k*)
    grammars
- [ ]   [Shunting yard algorithm](https://en.wikipedia.org/wiki/Shunting_yard_algorithm):
    convert an infix-notation math expression to postfix
- [ ]   [Pratt parser](https://en.wikipedia.org/wiki/Pratt_parser)
- [ ]   [Lexical analysis](https://en.wikipedia.org/wiki/Lexical_analysis)

### Quantum algorithms

- [ ]   [Deutsch-Jozsa algorithm](https://en.wikipedia.org/wiki/Deutsch-Jozsa_algorithm):
    criterion of balance for Boolean function
- [ ]   [Grover's algorithm](https://en.wikipedia.org/wiki/Grover's_algorithm): provides
    quadratic speedup for many search problems
- [ ]   [Shor's algorithm](https://en.wikipedia.org/wiki/Shor's_algorithm): provides
    [exponential](https://en.wikipedia.org/wiki/exponential_function) speedup (relative to
    currently known non-quantum algorithms) for factoring a number
- [ ]   [Simon's algorithm](https://en.wikipedia.org/wiki/Simon's_algorithm): provides a
    provably [exponential](https://en.wikipedia.org/wiki/exponential_function) speedup
    (relative to any non-quantum algorithm) for a black-box problem

### Theory of computation and automata

- [ ]   [Hopcroft's
    algorithm](https://en.wikipedia.org/wiki/DFA_minimization#Hopcroft's_algorithm),
    [Moore's algorithm](https://en.wikipedia.org/wiki/DFA_minimization#Moore's_algorithm),
    and [Brzozowski's
    algorithm](https://en.wikipedia.org/wiki/DFA_minimization#Brzozowski's_algorithm):
    algorithms for [minimizing the number of states in a deterministic
    finite automaton](https://en.wikipedia.org/wiki/DFA_minimization)
- [ ]   [Powerset construction](https://en.wikipedia.org/wiki/Powerset_construction): Algorithm
    to convert nondeterministic automaton to [deterministic
    automaton](https://en.wikipedia.org/wiki/deterministic_automaton).
- [ ]   [Tarski–Kuratowski
    algorithm](https://en.wikipedia.org/wiki/Tarski–Kuratowski_algorithm): a
    [non-deterministic
    algorithm](https://en.wikipedia.org/wiki/non-deterministic_algorithm) which provides an
    upper bound for the complexity of formulas in the [arithmetical
    hierarchy](https://en.wikipedia.org/wiki/arithmetical_hierarchy) and [analytical
    hierarchy](https://en.wikipedia.org/wiki/analytical_hierarchy)

Information theory and signal processing
----------------------------------------

### Coding theory

#### Error detection and correction

- [ ]   [BCH Codes](https://en.wikipedia.org/wiki/BCH_Code)
    - [ ]   [Berlekamp–Massey
        algorithm](https://en.wikipedia.org/wiki/Berlekamp–Massey_algorithm)
    - [ ]   [Peterson–Gorenstein–Zierler
        algorithm](https://en.wikipedia.org/wiki/Peterson–Gorenstein–Zierler_algorithm)
    - [ ]   [Reed–Solomon error
        correction](https://en.wikipedia.org/wiki/Reed–Solomon_error_correction)
- [ ]   [BCJR algorithm](https://en.wikipedia.org/wiki/BCJR_algorithm): decoding of error
    correcting codes defined on trellises (principally convolutional
    codes)
- [ ]   [Forward error correction](https://en.wikipedia.org/wiki/Forward_error_correction)
- [ ]   [Gray code](https://en.wikipedia.org/wiki/Gray_code)
- [ ]   [Hamming codes](https://en.wikipedia.org/wiki/Hamming_code)
    - [ ]   [Hamming(7,4)](Hamminghttps://en.wikipedia.org/wiki/(7,4)): a [Hamming
        code](https://en.wikipedia.org/wiki/Hamming_code) that encodes 4 bits of data into
        7 bits by adding 3 parity bits
    - [ ]   [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance): sum number of
        positions which are different
    - [ ]   [Hamming weight](https://en.wikipedia.org/wiki/Hamming_weight) (population count):
        find the number of 1 bits in a binary word
- [ ]   [Redundancy checks](https://en.wikipedia.org/wiki/Redundancy_check)
    - [ ]   [Adler-32](https://en.wikipedia.org/wiki/Adler-32)
    - [ ]   [Cyclic redundancy check](https://en.wikipedia.org/wiki/Cyclic_redundancy_check)
    - [ ]   [Damm algorithm](https://en.wikipedia.org/wiki/Damm_algorithm)
    - [ ]   [Fletcher's checksum](https://en.wikipedia.org/wiki/Fletcher's_checksum)
    - [ ]   [Longitudinal redundancy
        check](https://en.wikipedia.org/wiki/Longitudinal_redundancy_check) (LRC)
    - [ ]   [Luhn algorithm](https://en.wikipedia.org/wiki/Luhn_algorithm): a method of
        validating identification numbers
    - [ ]   [Luhn mod N algorithm](https://en.wikipedia.org/wiki/Luhn_mod_N_algorithm):
        extension of Luhn to non-numeric characters
    - [ ]   [Parity](https://en.wikipedia.org/wiki/Parity_bit): simple/fast error detection
        technique
    - [ ]   [Verhoeff algorithm](https://en.wikipedia.org/wiki/Verhoeff_algorithm)

#### Lossless compression algorithms

- [ ]   [Burrows–Wheeler transform](https://en.wikipedia.org/wiki/Burrows–Wheeler_transform):
    preprocessing useful for improving [lossless
    compression](https://en.wikipedia.org/wiki/Lossless_data_compression)
- [ ]   [Context tree weighting](https://en.wikipedia.org/wiki/Context_tree_weighting)
- [ ]   [Delta encoding](https://en.wikipedia.org/wiki/Delta_encoding): aid to compression of
    data in which sequential data occurs frequently
- [ ]   [Dynamic Markov compression](https://en.wikipedia.org/wiki/Dynamic_Markov_compression):
    Compression using predictive arithmetic coding
- [ ]   [Dictionary coders](https://en.wikipedia.org/wiki/Dictionary_coder)
    - [ ]   [Byte pair encoding](https://en.wikipedia.org/wiki/Byte_pair_encoding) (BPE)
    - [ ]   [DEFLATE](DEFLATE_https://en.wikipedia.org/wiki/(algorithm))
    - [ ]   [Lempel–Ziv](https://en.wikipedia.org/wiki/Lempel–Ziv)
        - [ ]   [LZ77 and LZ78](https://en.wikipedia.org/wiki/LZ77_and_LZ78)
        - [ ]   [Lempel–Ziv Jeff Bonwick](https://en.wikipedia.org/wiki/LZJB) (LZJB)
        - [ ]   [Lempel–Ziv–Markov chain
            algorithm](https://en.wikipedia.org/wiki/Lempel–Ziv–Markov_chain_algorithm)
            (LZMA)
        - [ ]   [Lempel–Ziv–Oberhumer](https://en.wikipedia.org/wiki/Lempel–Ziv–Oberhumer)
            (LZO): speed oriented
        - [ ]   [Lempel–Ziv–Stac](https://en.wikipedia.org/wiki/Lempel–Ziv–Stac) (LZS)
        - [ ]   [Lempel–Ziv–Storer–Szymanski](https://en.wikipedia.org/wiki/Lempel–Ziv–Storer–Szymanski)
            (LZSS)
        - [ ]   [Lempel–Ziv–Welch](https://en.wikipedia.org/wiki/Lempel–Ziv–Welch) (LZW)
        - [ ]   [LZWL](https://en.wikipedia.org/wiki/LZWL): syllable-based variant
        - [ ]   [LZX](https://en.wikipedia.org/wiki/LZX)
        - [ ]   [Lempel–Ziv Ross Williams](https://en.wikipedia.org/wiki/LZRW) (LZRW)
- [ ]   [Entropy encoding](https://en.wikipedia.org/wiki/Entropy_encoding): coding scheme that
    assigns codes to symbols so as to match code lengths with the
    probabilities of the symbols
    - [ ]   [Arithmetic coding](https://en.wikipedia.org/wiki/Arithmetic_coding): advanced
        [entropy](https://en.wikipedia.org/wiki/entropy) coding
        - [ ]   [Range encoding](https://en.wikipedia.org/wiki/Range_encoding): same as
            [arithmetic coding](https://en.wikipedia.org/wiki/arithmetic_coding), but
            looked at in a slightly different way
    - [ ]   [Huffman coding](https://en.wikipedia.org/wiki/Huffman_coding): simple lossless
        compression taking advantage of relative character frequencies
        - [ ]   [Adaptive Huffman
            coding](https://en.wikipedia.org/wiki/Adaptive_Huffman_coding): [adaptive
            coding](https://en.wikipedia.org/wiki/adaptive_coding) technique based on
            Huffman coding
        - [ ]   [Package-merge
            algorithm](https://en.wikipedia.org/wiki/Package-merge_algorithm): Optimizes
            Huffman coding subject to a length restriction on code
            strings
    - [ ]   [Shannon–Fano coding](https://en.wikipedia.org/wiki/Shannon–Fano_coding)
    - [ ]   [Shannon–Fano–Elias
        coding](https://en.wikipedia.org/wiki/Shannon–Fano–Elias_coding): precursor to
        arithmetic encoding[^1]
- [ ]   [Entropy coding with known entropy
    characteristics](https://en.wikipedia.org/wiki/Entropy_encoding)
    - [ ]   [Golomb coding](https://en.wikipedia.org/wiki/Golomb_coding): form of entropy
        coding that is optimal for alphabets following geometric
        distributions
    - [ ]   [Rice coding](https://en.wikipedia.org/wiki/Rice_coding): form of entropy coding
        that is optimal for alphabets following geometric distributions
    - [ ]   [Truncated binary
        encoding](https://en.wikipedia.org/wiki/Truncated_binary_encoding)
    - [ ]   [Unary coding](https://en.wikipedia.org/wiki/Unary_coding): code that represents a
        number n with n ones followed by a zero
    - [ ]   [Universal codes](Universal_code_https://en.wikipedia.org/wiki/(data_compression)):
        encodes positive integers into binary code words
        - [ ]   Elias [delta](https://en.wikipedia.org/wiki/Elias_delta_coding),
            [gamma](https://en.wikipedia.org/wiki/Elias_gamma_coding), and
            [omega](https://en.wikipedia.org/wiki/Elias_omega_coding) coding
        - [ ]   [Exponential-Golomb
            coding](https://en.wikipedia.org/wiki/Exponential-Golomb_coding)
        - [ ]   [Fibonacci coding](https://en.wikipedia.org/wiki/Fibonacci_coding)
        - [ ]   [Levenshtein coding](https://en.wikipedia.org/wiki/Levenshtein_coding)
- [ ]   [Fast Efficient & Lossless Image Compression
    System](https://en.wikipedia.org/wiki/FELICS) (FELICS): a lossless image compression
    algorithm
- [ ]   [Incremental encoding](https://en.wikipedia.org/wiki/Incremental_encoding): delta
    encoding applied to sequences of strings
- [ ]   [Prediction by partial
    matching](https://en.wikipedia.org/wiki/PPM_compression_algorithm) (PPM): an adaptive
    statistical data compression technique based on context modeling and
    prediction
- [ ]   [Run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding): lossless data
    compression taking advantage of strings of repeated characters
- [ ]   [SEQUITUR algorithm](https://en.wikipedia.org/wiki/SEQUITUR_algorithm): lossless
    compression by incremental grammar inference on a string

#### Lossy compression algorithms

- [ ]   [3Dc](https://en.wikipedia.org/wiki/3Dc): a lossy data compression algorithm for
    [normal maps](https://en.wikipedia.org/wiki/Normal_mapping)
- [ ]   [Audio](https://en.wikipedia.org/wiki/Audio_data_compression) and
    [Speech](https://en.wikipedia.org/wiki/speech_encoding) compression
    - [ ]   [A-law algorithm](https://en.wikipedia.org/wiki/A-law_algorithm): standard
        companding algorithm
    - [ ]   [Code-excited linear
        prediction](https://en.wikipedia.org/wiki/Code-excited_linear_prediction) (CELP):
        low bit-rate speech compression
    - [ ]   [Linear predictive coding](https://en.wikipedia.org/wiki/Linear_predictive_coding)
        (LPC): lossy compression by representing the [spectral
        envelope](https://en.wikipedia.org/wiki/spectral_envelope) of a digital signal of
        speech in compressed form
    - [ ]   [Mu-law algorithm](https://en.wikipedia.org/wiki/Mu-law_algorithm): standard analog
        signal compression or companding algorithm
    - [ ]   [Warped Linear Predictive
        Coding](https://en.wikipedia.org/wiki/Warped_Linear_Predictive_Coding) (WLPC)
- [ ]   [Image Compression](https://en.wikipedia.org/wiki/Image_Compression)
    - [ ]   [Block Truncation Coding](https://en.wikipedia.org/wiki/Block_Truncation_Coding)
        (BTC): a type of lossy image compression technique for greyscale
        images
    - [ ]   [Embedded Zerotree
        Wavelet](https://en.wikipedia.org/wiki/Embedded_Zerotree_Wavelet) (EZW)
    - [ ]   [Fast Cosine Transform
        algorithms](https://en.wikipedia.org/wiki/Fast_Cosine_Transform) (FCT algorithms):
        compute Discrete Cosine Transform (DCT) efficiently
    - [ ]   [Fractal compression](https://en.wikipedia.org/wiki/Fractal_compression): method
        used to compress images using fractals
    - [ ]   [Set Partitioning in Hierarchical
        Trees](https://en.wikipedia.org/wiki/Set_Partitioning_in_Hierarchical_Trees)
        (SPIHT)
    - [ ]   [Wavelet compression](https://en.wikipedia.org/wiki/Wavelet_compression): form of
        data compression well suited for [image
        compression](https://en.wikipedia.org/wiki/image_compression) (sometimes also video
        compression and audio compression)
- [ ]   [Transform coding](https://en.wikipedia.org/wiki/Transform_coding): type of data
    compression for “natural” data like audio signals or photographic
    images
- [ ]   [Video compression](https://en.wikipedia.org/wiki/Video_compression)
- [ ]   [Vector quantization](https://en.wikipedia.org/wiki/Vector_quantization): technique
    often used in lossy data compression

### Digital signal processing

- [ ]   [Adaptive-additive
    algorithm](https://en.wikipedia.org/wiki/Adaptive-additive_algorithm) (AA algorithm):
    find the spatial frequency phase of an observed wave source
- [ ]   [Discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform):
    determines the frequencies contained in a (segment of a) signal
    - [ ]   [Bluestein's FFT
        algorithm](https://en.wikipedia.org/wiki/Bluestein's_FFT_algorithm)
    - [ ]   [Bruun's FFT algorithm](https://en.wikipedia.org/wiki/Bruun's_FFT_algorithm)
    - [ ]   [Cooley&ndash;Tukey FFT
        algorithm](https://en.wikipedia.org/wiki/Cooley&ndash;Tukey_FFT_algorithm)
    - [ ]   [Fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform)
    - [ ]   [Prime-factor FFT
        algorithm](https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm)
    - [ ]   [Rader's FFT algorithm](https://en.wikipedia.org/wiki/Rader's_FFT_algorithm)
- [ ]   [Fast folding algorithm](https://en.wikipedia.org/wiki/Fast_folding_algorithm): an
    efficient algorithm for the detection of approximately periodic
    events within time series data
- [ ]   [Gerchberg–Saxton algorithm](https://en.wikipedia.org/wiki/Gerchberg–Saxton_algorithm):
    Phase retrieval algorithm for optical planes
- [ ]   [Goertzel algorithm](https://en.wikipedia.org/wiki/Goertzel_algorithm): identify a
    particular frequency component in a signal. Can be used for
    [DTMF](https://en.wikipedia.org/wiki/DTMF) digit decoding.
- [ ]   [Karplus-Strong string
    synthesis](https://en.wikipedia.org/wiki/Karplus-Strong_string_synthesis): physical
    modelling synthesis to simulate the sound of a hammered or plucked
    string or some types of percussion

#### Image processing

- [ ]   Contrast Enhancement
    - [ ]   [Histogram equalization](https://en.wikipedia.org/wiki/Histogram_equalization): use
        histogram to improve image contrast
    - [ ]   [Adaptive histogram
        equalization](https://en.wikipedia.org/wiki/Adaptive_histogram_equalization):
        histogram equalization which adapts to local changes in contrast
- [ ]   [Connected-component
    labeling](https://en.wikipedia.org/wiki/Connected-component_labeling): find and label
    disjoint regions
- [ ]   [Dithering](https://en.wikipedia.org/wiki/Dithering) and
    [half-toning](https://en.wikipedia.org/wiki/half-toning)
    - [ ]   [Error diffusion](https://en.wikipedia.org/wiki/Error_diffusion)
    - [ ]   [Floyd–Steinberg
        dithering](https://en.wikipedia.org/wiki/Floyd–Steinberg_dithering)
    - [ ]   [Ordered dithering](https://en.wikipedia.org/wiki/Ordered_dithering)
    - [ ]   [Riemersma dithering](https://en.wikipedia.org/wiki/Riemersma_dithering)
- [ ]   Elser [difference-map
    algorithm](https://en.wikipedia.org/wiki/difference-map_algorithm): a search algorithm
    for general constraint satisfaction problems. Originally used for
    [X-Ray diffraction](https://en.wikipedia.org/wiki/X-ray_crystallography) microscopy
- [ ]   [Feature detection](Feature_detection_https://en.wikipedia.org/wiki/(computer_vision))
    - [ ]   [Canny edge detector](https://en.wikipedia.org/wiki/Canny_edge_detector): detect a
        wide range of edges in images
    - [ ]   [Generalised Hough
        transform](https://en.wikipedia.org/wiki/Generalised_Hough_transform)
    - [ ]   [Hough transform](https://en.wikipedia.org/wiki/Hough_transform)
    - [ ]   [Marr–Hildreth algorithm](https://en.wikipedia.org/wiki/Marr–Hildreth_algorithm):
        an early [edge detection](https://en.wikipedia.org/wiki/edge_detection) algorithm
    - [ ]   [SIFT](https://en.wikipedia.org/wiki/Scale-invariant_feature_transform)
        (Scale-invariant feature transform): is an algorithm to detect
        and describe local features in images.
    - [ ]   [SURF (Speeded Up Robust
        Features)](SURF_https://en.wikipedia.org/wiki/(Speeded_Up_Robust_Features)): is a
        robust local feature detector, first presented by Herbert Bay et
        al. in 2006, that can be used in computer vision tasks like
        object recognition or 3D reconstruction. It is partly inspired
        by the SIFT descriptor. The standard version of SURF is several
        times faster than SIFT and claimed by its authors to be more
        robust against different image transformations than
        SIFT.[^2][^3]
- [ ]   [Richardson–Lucy
    deconvolution](https://en.wikipedia.org/wiki/Richardson–Lucy_deconvolution): image
    de-blurring algorithm
- [ ]   [Blind deconvolution](https://en.wikipedia.org/wiki/Blind_deconvolution): image
    de-blurring algorithm when [point spread
    function](https://en.wikipedia.org/wiki/point_spread_function) is unknown.
- [ ]   [Median filtering](https://en.wikipedia.org/wiki/Median_filtering)
- [ ]   [Seam carving](https://en.wikipedia.org/wiki/Seam_carving): content-aware image
    resizing algorithm
- [ ]   [Segmentation](Segmentation_https://en.wikipedia.org/wiki/(image_processing)):
    partition a digital image into two or more regions
    - [ ]   [GrowCut algorithm](https://en.wikipedia.org/wiki/GrowCut_algorithm): an
        interactive segmentation algorithm
    - [ ]   [Random walker algorithm](https://en.wikipedia.org/wiki/Random_walker_algorithm)
    - [ ]   [Region growing](https://en.wikipedia.org/wiki/Region_growing)
    - [ ]   [Watershed transformation](Watershed_https://en.wikipedia.org/wiki/(algorithm)): a
        class of algorithms based on the watershed analogy

Software engineering
--------------------

- [ ]   [Cache algorithms](https://en.wikipedia.org/wiki/Cache_algorithms)
- [ ]   [CHS conversion](https://en.wikipedia.org/wiki/CHS_conversion): converting between disk
    addressing systems
- [ ]   [Double dabble](https://en.wikipedia.org/wiki/Double_dabble): Convert binary numbers to
    BCD
- [ ]   [Hash Function](https://en.wikipedia.org/wiki/Hash_Function): convert a large, possibly
    variable-sized amount of data into a small datum, usually a single
    integer that may serve as an index into an array
    - [ ]   [Fowler–Noll–Vo hash
        function](https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function): fast with
        low collision rate
    - [ ]   [Pearson hashing](https://en.wikipedia.org/wiki/Pearson_hashing): computes 8 bit
        value only, optimized for 8 bit computers
    - [ ]   [Zobrist hashing](https://en.wikipedia.org/wiki/Zobrist_hashing): used in the
        implementation of [transposition
        tables](https://en.wikipedia.org/wiki/transposition_table)
- [ ]   [Unicode Collation
    Algorithm](https://en.wikipedia.org/wiki/Unicode_Collation_Algorithm)
- [ ]   [Xor swap algorithm](https://en.wikipedia.org/wiki/Xor_swap_algorithm): swaps the
    values of two variables without using a buffer

Database algorithms
-------------------

- [ ]   [Algorithms for Recovery and Isolation Exploiting
    Semantics](https://en.wikipedia.org/wiki/Algorithms_for_Recovery_and_Isolation_Exploiting_Semantics)
    (ARIES): [transaction](transaction_https://en.wikipedia.org/wiki/(database)) recovery
- [ ]   [Join algorithms](Join_https://en.wikipedia.org/wiki/(SQL))
    - [ ]   [Block nested loop](https://en.wikipedia.org/wiki/Block_nested_loop)
    - [ ]   [Hash join](https://en.wikipedia.org/wiki/Hash_join)
    - [ ]   [Nested loop join](https://en.wikipedia.org/wiki/Nested_loop_join)
    - [ ]   [Sort-Merge Join](https://en.wikipedia.org/wiki/Sort-Merge_Join)

Distributed systems algorithms
------------------------------

- [ ]   [Bully algorithm](https://en.wikipedia.org/wiki/Bully_algorithm): a method for
    dynamically selecting a coordinator
- [ ]   [Byzantine fault tolerance](https://en.wikipedia.org/wiki/Byzantine_fault_tolerance):
    good [fault tolerance](https://en.wikipedia.org/wiki/Fault-tolerant_system).
- [ ]   [Clock synchronization](https://en.wikipedia.org/wiki/Clock_synchronization)
    - [ ]   [Berkeley algorithm](https://en.wikipedia.org/wiki/Berkeley_algorithm)
    - [ ]   [Cristian's algorithm](https://en.wikipedia.org/wiki/Cristian's_algorithm)
    - [ ]   [Intersection algorithm](https://en.wikipedia.org/wiki/Intersection_algorithm)
    - [ ]   [Marzullo's algorithm](https://en.wikipedia.org/wiki/Marzullo's_algorithm)
- [ ]   Detection of Process Termination
    - [ ]   [Dijkstra-Scholten
        algorithm](https://en.wikipedia.org/wiki/Dijkstra-Scholten_algorithm)
    - [ ]   [Huang's algorithm](https://en.wikipedia.org/wiki/Huang's_algorithm)
- [ ]   [Lamport ordering](https://en.wikipedia.org/wiki/Lamport_ordering): a [partial
    ordering](https://en.wikipedia.org/wiki/partial_order) of events based on the
    *happened-before* relation
- [ ]   [Mutual exclusion](https://en.wikipedia.org/wiki/Mutual_exclusion)
    - [ ]   [Lamport's Distributed Mutual Exclusion
        Algorithm](https://en.wikipedia.org/wiki/Lamport's_Distributed_Mutual_Exclusion_Algorithm)
    - [ ]   [Naimi-Trehel's log(n)
        Algorithm](https://en.wikipedia.org/wiki/Naimi-Trehel's_log(n)_Algorithm)
    - [ ]   [Maekawa's Algorithm](https://en.wikipedia.org/wiki/Maekawa's_Algorithm)
    - [ ]   [Raymond's Algorithm](https://en.wikipedia.org/wiki/Raymond's_Algorithm)
    - [ ]   [Ricart-Agrawala
        Algorithm](https://en.wikipedia.org/wiki/Ricart-Agrawala_Algorithm)
- [ ]   [Paxos algorithm](https://en.wikipedia.org/wiki/Paxos_algorithm): a family of protocols
    for solving consensus in a network of unreliable processors
- [ ]   [Snapshot algorithm](https://en.wikipedia.org/wiki/Snapshot_algorithm): record a
    consistent global state for an asynchronous system
    - [ ]   [Chandy-Lamport algorithm](https://en.wikipedia.org/wiki/Chandy-Lamport_algorithm)
- [ ]   [Vector clocks](https://en.wikipedia.org/wiki/Vector_clocks): generate a [partial
    ordering](https://en.wikipedia.org/wiki/partial_ordering) of events in a distributed
    system and detect [causality](https://en.wikipedia.org/wiki/causality) violations

### Memory allocation and deallocation algorithms

- [ ]   [Buddy memory allocation](https://en.wikipedia.org/wiki/Buddy_memory_allocation):
    Algorithm to allocate memory such that fragmentation is less.
- [ ]   [Garbage
    collectors](Garbage_collection_https://en.wikipedia.org/wiki/(computer_science))
    - [ ]   [Cheney's algorithm](https://en.wikipedia.org/wiki/Cheney's_algorithm): An
        improvement on the [Semi-space
        collector](https://en.wikipedia.org/wiki/Semi-space_collector)
    - [ ]   [Generational garbage
        collector](garbage_collection_https://en.wikipedia.org/wiki/(computer_science)):
        Fast garbage collectors that segregate memory by age
    - [ ]   [Mark-compact algorithm](https://en.wikipedia.org/wiki/Mark-compact_algorithm): a
        combination of the [mark-sweep
        algorithm](https://en.wikipedia.org/wiki/Mark_and_sweep) and [Cheney's copying
        algorithm](https://en.wikipedia.org/wiki/Cheney's_algorithm)
    - [ ]   [Mark and sweep](https://en.wikipedia.org/wiki/Mark_and_sweep)
    - [ ]   [Semi-space collector](https://en.wikipedia.org/wiki/Semi-space_collector): An
        early copying collector
- [ ]   [Reference counting](https://en.wikipedia.org/wiki/Reference_counting)

Networking
----------

- [ ]   [Karn's algorithm](https://en.wikipedia.org/wiki/Karn's_algorithm): addresses the
    problem of getting accurate estimates of the round-trip time for
    messages when using TCP
- [ ]   [Luleå algorithm](https://en.wikipedia.org/wiki/Luleå_algorithm): a technique for
    storing and searching internet routing tables efficiently
- [ ]   [Network congestion](https://en.wikipedia.org/wiki/Network_congestion)
    - [ ]   [Exponential backoff](https://en.wikipedia.org/wiki/Exponential_backoff)
    - [ ]   [Nagle's algorithm](https://en.wikipedia.org/wiki/Nagle's_algorithm): improve the
        efficiency of TCP/IP networks by coalescing packets
    - [ ]   [Truncated binary exponential
        backoff](https://en.wikipedia.org/wiki/Truncated_binary_exponential_backoff)

Operating systems algorithms
----------------------------

- [ ]   [Banker's algorithm](https://en.wikipedia.org/wiki/Banker's_algorithm): Algorithm used
    for deadlock avoidance.
- [ ]   [Page replacement
    algorithms](https://en.wikipedia.org/wiki/Page_replacement_algorithms): Selecting the
    victim page under low memory conditions.
    - [ ]   [Adaptive replacement
        cache](https://en.wikipedia.org/wiki/Adaptive_replacement_cache): better
        performance than LRU
    - [ ]   [Clock with Adaptive
        Replacement](https://en.wikipedia.org/wiki/Clock_with_Adaptive_Replacement) (CAR):
        is a page replacement algorithm that has performance comparable
        to [Adaptive replacement
        cache](https://en.wikipedia.org/wiki/Adaptive_replacement_cache)

### Process synchronization

- [ ]   [Dekker's algorithm](https://en.wikipedia.org/wiki/Dekker's_algorithm)
- [ ]   [Lamport's Bakery algorithm](https://en.wikipedia.org/wiki/Lamport's_Bakery_algorithm)
- [ ]   [Peterson's algorithm](https://en.wikipedia.org/wiki/Peterson's_algorithm)

### Scheduling

- [ ]   [Earliest deadline first
    scheduling](https://en.wikipedia.org/wiki/Earliest_deadline_first_scheduling)
- [ ]   [Fair-share scheduling](https://en.wikipedia.org/wiki/Fair-share_scheduling)
- [ ]   [Least slack time
    scheduling](https://en.wikipedia.org/wiki/Least_slack_time_scheduling)
- [ ]   [List scheduling](https://en.wikipedia.org/wiki/List_scheduling)
- [ ]   [Multi level feedback queue](https://en.wikipedia.org/wiki/Multi_level_feedback_queue)
- [ ]   [Rate-monotonic scheduling](https://en.wikipedia.org/wiki/Rate-monotonic_scheduling)
- [ ]   [Round-robin scheduling](https://en.wikipedia.org/wiki/Round-robin_scheduling)
- [ ]   [Shortest job next](https://en.wikipedia.org/wiki/Shortest_job_next)
- [ ]   [Shortest remaining time](https://en.wikipedia.org/wiki/Shortest_remaining_time)
- [ ]   [Top-nodes algorithm](https://en.wikipedia.org/wiki/Top-nodes_algorithm): resource
    calendar management

#### Disk scheduling

- [ ]   [Elevator algorithm](https://en.wikipedia.org/wiki/Elevator_algorithm): Disk scheduling
    algorithm that works like an elevator.
- [ ]   [Shortest seek first](https://en.wikipedia.org/wiki/Shortest_seek_first): Disk
    scheduling algorithm to reduce [seek time](https://en.wikipedia.org/wiki/seek_time).
