# Simulating Bond Percolation

## Background

[Percolation theory](https://en.wikipedia.org/wiki/Percolation_theory) is the study of the behaviour of the random subgraphs obtained when adding edges to a lattice with a given probability $p$. More precisely, what we have just described is _bond percolation_ - in the alternate model of _site percolation_ nodes are occupied with a given probability, with edges then implicitly drawn to all neighbouring occupied nodes. Unless otherwise specified, we will always be referring to bond percolation from now on.

We are typically interested in studying the behaviour of the connected components of the percolation as $p$ varies. For the case of infinite lattices, we ask the question: Are there infinite clusters? If so, how many, and how does this depend on $p$?

Obviously, for $p=0$, every cluster has size 1, and for $p=1$, there is one single infinite cluster. By Kolmogorov's zero–one law, the probability of there being an infinite cluster for a given $p$ is either 0 or 1. It is not hard to show that this is an increasing function of $p$. We therefore obtain that there must be a critical value $p_c$ such that for every $p<p_c$, there are no infinite clusters, and for $p>p_c$ there is always an infinite cluster.


## Simulation

Without referring to any prior research in the field, two obvious ways of simulating this behaviour stood out to me.

The first would be a way of simulating the growth of a single cluster: Starting at a given point, add all of the neighbours (along a valid edge in the lattice) to a queue each with probability $p$, then add the point to the visited set. Continue the algorithm by popping a neighbour off the queue, adding all unvisited neighbours to the back of the queue with probablity $p$, then adding the current point to the visited set and so on. If the cluster was finite, the algorithm would terminate when the queue of neighbours was empty.

With this algorithm we obtain $O(N)$ runtime and memory usage.

The second would be to use a [disjoint set forest](https://en.wikipedia.org/wiki/Disjoint-set_data_structure) data structure: We iterate over every possible edge in the lattice and merge the two clusters at each point with probability $p$.

With this method (using union by size and path halving) we obtain $O(N)$ memory usage and $O(N\alpha(N))$ runtime (where $\alpha$ is the inverse Ackermann function, which can be considered essentially constant for our purposes).

The second method has the obvious advantage of simulating the entire system, not just a single cluster, therefore we shall choose the latter as both offer comparable performance.
