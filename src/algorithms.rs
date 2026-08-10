//! Algorithms using GBWT and GBZ.

use crate::{GBZ, Orientation};
use crate::support::{self, Chains};

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::fmt::{Display, Formatter};
use std::{cmp, fmt};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Returns the longest common subsequence of integer sequences `a` and `b`, weighted by the given function.
///
/// The subsequence is returned as pairs of positions, and the second return value is the total weight of the LCS.
/// Weights are applied to the elements of the sequences, and they are assumed to be non-zero.
/// This version of the algorithm uses naive dynamic programming.
/// It is not suitable for long sequences.
/// See [`fast_weighted_lcs`] for an implementation based on Myers' O(nd) algorithm.
pub fn naive_weighted_lcs<F: Fn(usize) -> usize>(a: &[usize], b: &[usize], weight: F) -> (Vec<(usize, usize)>, usize) {
    let mut dp = vec![vec![0; b.len() + 1]; a.len() + 1];
    for (i, a_val) in a.iter().enumerate() {
        for (j, b_val) in b.iter().enumerate() {
            dp[i + 1][j + 1] = cmp::max(dp[i + 1][j], dp[i][j + 1]);
            if a_val == b_val {
                dp[i + 1][j + 1] = cmp::max(dp[i + 1][j + 1], dp[i][j] + weight(a[i]));
            }
        }
    }

    let mut result = Vec::new();
    let mut i = a.len();
    let mut j = b.len();
    while i > 0 && j > 0 {
        if a[i - 1] == b[j - 1] {
            result.push((i - 1, j - 1));
            i -= 1;
            j -= 1;
        } else if dp[i][j] == dp[i - 1][j] {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    result.reverse();
    (result, dp[a.len()][b.len()])
}

//-----------------------------------------------------------------------------

// TODO: We do not need to store the weight. Due to the invariant, we can derive it from
// edits and the prefix sums.
// TODO: We could make this more space-efficient by using 32-bit integers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct DPPoint {
    // Twice the total weight of the LCS up to this point.
    // This guarantees an invariant: `weight + edits = a_prefix_sum + b_prefix_sum`.
    weight: usize,

    // Non-weighted offset in the first sequence.
    a_offset: usize,

    // Non-weighted offset in the second sequence.
    b_offset: usize,

    // Length of the non-weighted run of matches.
    matches: usize,
}

impl DPPoint {
    fn new(weight: usize, a_offset: usize, b_offset: usize) -> DPPoint {
        DPPoint {
            weight, a_offset, b_offset, matches: 0,
        }
    }
}

impl Display for DPPoint {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "(weight {}, a {}, b {}, matches {})", self.weight, self.a_offset, self.b_offset, self.matches)
    }
}

struct DPMatrix<'a> {
    // First sequence.
    a: &'a [usize],

    // Second sequence.
    b: &'a [usize],

    // Prefix sum of weights on the first sequence.
    a_weights: Vec<usize>,

    // Prefix sum of weights on the second sequence.
    b_weights: Vec<usize>,

    // Furthest point for the given number of weighted edits on the weighted diagonal.
    // A diagonal is a_prefix_sum - b_prefix_sum.
    // The point refers to the first unprocessed pair of offsets.
    points: BTreeMap<(usize, isize), DPPoint>,
}

impl<'a> DPMatrix<'a> {
    fn prefix_sum<F: Fn(usize) -> usize>(sequence: &[usize], weight: F) -> Vec<usize> {
        let mut result = Vec::with_capacity(sequence.len() + 1);
        result.push(0);
        for i in 0..sequence.len() {
            result.push(result[i] + weight(sequence[i]));
        }
        result
    }

    fn new<F: Fn(usize) -> usize>(a: &'a [usize], b: &'a [usize], weight: F) -> Self {
        let a_weights = Self::prefix_sum(a, &weight);
        let b_weights = Self::prefix_sum(b, &weight);
        let mut points = BTreeMap::new();
        points.insert((0, 0), DPPoint::new(0, 0, 0));
        DPMatrix {
            a, b, a_weights, b_weights, points
        }
    }

    fn diagonals_for(&self, edits: usize) -> Vec<isize> {
        let mut result = Vec::new();
        for ((_, diagonal), _) in self.points.range((edits, isize::MIN)..(edits + 1, isize::MIN)) {
            result.push(*diagonal);
        }
        result
    }

    // Inserts the point if it is better than the existing one.
    fn try_insert(&mut self, edits: usize, diagonal: isize, point: DPPoint) {
        if let Some(existing) = self.points.get_mut(&(edits, diagonal)) {
            if point.weight > existing.weight {
                *existing = point;
            }
        } else {
            self.points.insert((edits, diagonal), point);
        }
    }

    fn a_weight(&self, offset: usize) -> usize {
        self.a_weights[offset + 1] - self.a_weights[offset]
    }

    fn b_weight(&self, offset: usize) -> usize {
        self.b_weights[offset + 1] - self.b_weights[offset]
    }

    // Extends matches over all diagonals for the given number of weighted edits.
    // Also adds successor points reachable with a single edit from each extension.
    // Returns the point if we reached the end.
    fn extend(&mut self, edits: usize) -> Option<DPPoint> {
        let diagonals = self.diagonals_for(edits);
        for diagonal in diagonals.into_iter() {
            let point = self.points.get(&(edits, diagonal)).copied();
            if point.is_none() {
                continue;
            }
            let mut point = point.unwrap();
            while point.a_offset < self.a.len() && point.b_offset < self.b.len() && self.a[point.a_offset] == self.b[point.b_offset] {
                point.weight += 2 * self.a_weight(point.a_offset);
                point.a_offset += 1;
                point.b_offset += 1;
                point.matches += 1;
            }
            if point.matches > 0 {
                self.points.insert((edits, diagonal), point);
            }
            if point.a_offset == self.a.len() && point.b_offset == self.b.len() {
                return Some(point);
            }
            if point.a_offset < self.a.len() {
                let weight = self.a_weight(point.a_offset);
                let new_point = DPPoint::new(point.weight, point.a_offset + 1, point.b_offset);
                self.try_insert(edits + weight, diagonal + (weight as isize), new_point);
            }
            if point.b_offset < self.b.len() {
                let weight = self.b_weight(point.b_offset);
                let new_point = DPPoint::new(point.weight, point.a_offset, point.b_offset + 1);
                self.try_insert(edits + weight, diagonal - (weight as isize), new_point);
            }
        }
        None
    }

    // Returns the next possible number of edits after the given number.
    // This assumes that `extend` has been called for the given number of edits.
    fn next_edits(&self, edits: usize) -> Option<usize> {
        if let Some(((value, _), _)) = self.points.range((edits + 1, isize::MIN)..).next() {
            Some(*value)
        } else {
            None
        }
    }

    // Returns the predecessor point and number of edits for the given point and number of edits.
    fn predecessor(&self, a_offset: usize, b_offset: usize, edits: usize) -> Option<(DPPoint, usize)> {
        let diagonal = (self.a_weights[a_offset] as isize) - (self.b_weights[b_offset] as isize);
        let prev = if a_offset > 0 && self.a_weight(a_offset - 1) <= edits {
            let weight = self.a_weight(a_offset - 1);
            self.points.get(&(edits - weight, diagonal - (weight as isize))).copied()
        } else {
            None
        };
        let next = if b_offset > 0 && self.b_weight(b_offset - 1) <= edits {
            let weight = self.b_weight(b_offset - 1);
            self.points.get(&(edits - weight, diagonal + (weight as isize))).copied()
        } else {
            None
        };
        match (prev, next) {
            (Some(p), Some(n)) => {
                if p.weight > n.weight {
                    Some((p, edits - self.a_weight(a_offset - 1)))
                } else {
                    Some((n, edits - self.b_weight(b_offset - 1)))
                }
            },
            (Some(p), None) => Some((p, edits - self.a_weight(a_offset - 1))),
            (None, Some(n)) => Some((n, edits - self.b_weight(b_offset - 1))),
            (None, None) => None,
        }
    }
}

/// Returns the longest common subsequence of integer sequences `a` and `b`, weighted by the given function.
///
/// The subsequence is returned as pairs of positions, and the second return value is the total weight of the LCS.
/// Weights are applied to the elements of the sequences, and they are assumed to be non-zero.
/// This version of the algorithm is based on Myers' O(nd) algorithm.
/// It is efficient with long sequences, as long as they are similar.
/// See [`naive_weighted_lcs`] for an implementation using naive dynamic programming.
///
/// The following specializations are available:
///
/// * [`lcs`] for unweighted LCS.
/// * [`path_lcs`] for paths in a graph, using sequence lengths as weights.
pub fn fast_weighted_lcs<F: Fn(usize) -> usize>(a: &[usize], b: &[usize], weight: F) -> (Vec<(usize, usize)>, usize) {
    if a.is_empty() || b.is_empty() {
        return (Vec::new(), 0);
    }

    // Find the furthest point on each diagonal with each possible number of edits, until we reach the end.
    let mut matrix = DPMatrix::new(a, b, &weight);
    let mut edits = 0;
    let mut point = DPPoint::new(0, 0, 0);
    loop {
        if let Some(next_point) = matrix.extend(edits) {
            point = next_point;
            break;
        }
        if let Some(next_edits) = matrix.next_edits(edits) {
            edits = next_edits;
        } else {
            // This should not happen.
            break;
        }
    }

    // Trace back the LCS.
    let mut result = Vec::new();
    let final_weight = point.weight / 2;
    loop {
        for _ in 0..point.matches {
            point.a_offset -= 1;
            point.b_offset -= 1;
            result.push((point.a_offset, point.b_offset));
        }
        if let Some((p, e)) = matrix.predecessor(point.a_offset, point.b_offset, edits) {
            point = p;
            edits = e;
        } else {
            break;
        }
    }
    result.reverse();
    (result, final_weight)
}

//-----------------------------------------------------------------------------

/// Returns the longest common subsequence of `a` and `b`.
///
/// The return value consists of pairs of positions in the input vectors.
///
/// # Examples
///
/// ```
/// use gbz::algorithms::lcs;
///
/// let a = vec![1, 2, 3, 4, 5];
/// let b = vec![2, 4, 6, 8, 10];
/// let truth = vec![(1, 0), (3, 1)];
/// assert_eq!(lcs(&a, &b), truth);
/// ``````
pub fn lcs(a: &[usize], b: &[usize]) -> Vec<(usize, usize)> {
    fast_weighted_lcs(a, b, |_| 1).0
}

/// Returns the longest common subsequence of paths `a` and `b` in the graph.
///
/// The LCS is weighted by the length of the node and returned as a vector of pairs of positions.
/// The second return value is the total weight of the LCS.
/// If a node is not found in the graph, its length is assumed to be zero.
/// In such cases, the LCS may not be meaningful.
///
/// # Examples
///
/// ```
/// use gbz::{GBZ, Orientation, algorithms};
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// fn get_path(graph: &GBZ, path_id: usize) -> Vec<usize> {
///     graph.path(path_id, Orientation::Forward)
///         .unwrap()
///         .map(|(id, o)| support::encode_node(id, o))
///         .collect()
/// }
///
/// // (1: GA), (2: T), (3: T), (5: CA), (6: G), (9: A), (11: TA)
/// let a = get_path(&gbz, 1);
///
/// // (1: GA), (2: T), (4: A), (5: CA), (6: G), (10: T), (11: TA)
/// let b = get_path(&gbz, 2);
///
/// let truth = vec![
///     (0, 0), (1, 1), (3, 3), (4, 4), (6, 6)
/// ];
/// let len = 8;
///
/// assert_eq!(algorithms::path_lcs(&a, &b, &gbz), (truth, len));
/// ```
pub fn path_lcs(a: &[usize], b: &[usize], graph: &GBZ) -> (Vec<(usize, usize)>, usize) {
    let weight = |handle| graph.sequence_len(support::node_id(handle)).unwrap_or(0);
    fast_weighted_lcs(a, b, weight)
}

//-----------------------------------------------------------------------------

// TODO: chains.rs with the `Chains` struct and related code.
// TODO: We could have a version that returns the actual chains (including the non-boundary nodes?).
// We could also compute the minimum distance (in bp) from the start of the chain to each node.
// That would be good enough to replace the distance index in haplotype sampling.
// TODO: Chain as a vector of (handle, is_snarl_start, is_snarl_end, distance_from_start, distance_to_end) tuples?
// TODO: Process the components in parallel?

/// Finds top-level chains in the graph.
///
/// A top-level chain is a sequence of boundary nodes bordering snarls in a weakly connected component.
/// The result contains links between consecutive boundary nodes in each chain.
///
/// Each component must have exactly two tip nodes, and there must be a directed path between them.
/// If these assumptions do not hold, the component is skipped and not included in the result.
/// In more complex cases, the chains can be extracted from a snarl decomposition or a distance index using `vg chains`.
///
/// The algorithm uses 32-bit integers internally to save memory.
/// This is not a serious limitation, as GBWTs in the vg ecosystem also use 32-bit node identifiers in many places.
///
/// NOTE: The integrated tests using `cargo test` only cover small instances.
/// If the construction logic is modified, it should be tested with the `build-chains` binary against existing chains files.
///
/// # Examples
///
/// ```
/// use gbz::{GBZ, Orientation};
/// use gbz::{algorithms, support};
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// let chains = algorithms::find_chains(&gbz);
/// assert_eq!(chains.len(), 1);
/// assert_eq!(chains.components(), Some(1));
/// assert_eq!(chains.links(), 3);
/// assert_eq!(
///     chains.next(support::encode_node(2, Orientation::Forward)),
///     Some(support::encode_node(5, Orientation::Forward))
/// );
/// assert_eq!(
///     chains.next(support::encode_node(5, Orientation::Forward)),
///     Some(support::encode_node(6, Orientation::Forward))
/// );
/// assert_eq!(
///    chains.next(support::encode_node(6, Orientation::Forward)),
///    Some(support::encode_node(11, Orientation::Forward))
/// );
/// ```
pub fn find_chains(graph: &GBZ) -> Chains {
    let mut chains = Chains::new();
    let components = graph.weakly_connected_components();

    let mut trivial_chains = 0;
    for component in components.iter() {
        let tips = find_tips(graph, component);
        if tips.len() != 2 {
            continue;
        }
        let from_handle = tips[0];
        let to_handle = support::flip_node(tips[1]);
        let bridges = find_bridge_nodes(graph, component);

        // Determine all handles in the top-level chain and which of them are snarl entry points.
        let result = walk_chain(graph, from_handle, to_handle, &bridges);
        if result.is_none() {
            continue;
        }
        let chain_handles = result.unwrap();
        drop(bridges);

        // Returns `true` if the handle is a snarl entry point.
        // Either the outdegree is > 1 or the indegree of the only successor is > 1.
        let is_snarl_entry = |handle: u32| -> bool {
            let mut outdegree = 0;
            let mut first_successor = None;
            if let Some(iter) = graph.successors(support::node_id(handle as usize), support::node_orientation(handle as usize)) {
                for (next_id, next_o) in iter {
                    outdegree += 1;
                    if outdegree > 1 {
                        return true;
                    }
                    first_successor = Some((next_id, next_o));
                }
            }
            if let Some((id, orientation)) = first_successor {
                match graph.indegree(id, orientation) {
                    Some(degree) => degree > 1,
                    None => false,
                }
            } else {
                false
            }
        };

        // Returns `true` if the handle is a boundary node of a snarl.
        let is_boundary_node = |handle: u32| -> bool {
            is_snarl_entry(handle) || is_snarl_entry(support::flip_node(handle as usize) as u32)
        };

        // Add links between consecutive boundary nodes.
        let boundary: Vec<u32> = chain_handles.iter().copied()
            .filter(|handle| is_boundary_node(*handle))
            .collect();
        drop(chain_handles);
        if boundary.len() < 2 {
            trivial_chains += 1;
        }
        for pair in boundary.windows(2) {
            chains.add_link(pair[0] as usize, pair[1] as usize);
        }
    }

    chains.count_chains();
    chains.set_trivial_chains(Some(trivial_chains));
    chains.set_components(Some(components.len()));
    chains
}

// Returns inward-oriented handles for tip nodes in the component.
// A tip node has successors on one side but not the other.
fn find_tips(graph: &GBZ, component: &[usize]) -> Vec<usize> {
    let mut tips = Vec::new();
    for &node_id in component {
        let fwd_succ = graph.successors(node_id, Orientation::Forward)
            .map_or(0, |iter| iter.count());
        let fwd_pred = graph.predecessors(node_id, Orientation::Forward)
            .map_or(0, |iter| iter.count());

        if fwd_succ == 0 && fwd_pred == 0 {
            continue;
        }
        if fwd_succ > 0 && fwd_pred > 0 {
            continue;
        }
        // Tip: has edges on one side only. Inward handle is the orientation with successors.
        if fwd_succ > 0 {
            tips.push(support::encode_node(node_id, Orientation::Forward));
        } else {
            tips.push(support::encode_node(node_id, Orientation::Reverse));
        }
    }
    tips
}

// Finds bridge nodes using a variant of Tarjan's algorithm.
// Returns node ids (not handles).
fn find_bridge_nodes(graph: &GBZ, component: &[usize]) -> HashSet<u32> {
    // Map from node ids to indices in the component vector.
    let mut node_to_index: HashMap<u32, u32> = HashMap::new();
    for (i, &node_id) in component.iter().enumerate() {
        node_to_index.insert(node_id as u32, i as u32);
    }
    let node_to_index = node_to_index;

    // These arrays are indexed by encoded node sides, except that we use indexes in
    // the component vector in place of node ids. Each node side is effectively
    // considered a separate node.
    let n: usize = 2 * component.len();
    let mut visited = vec![false; n];
    let mut preorder_rank = vec![0u32; n];
    let mut lowest_reachable = vec![0u32; n];
    let mut parent = vec![u32::MAX; n];
    let mut result = HashSet::new();
    let mut timer = 0u32;

    // Returns the `i`-th neighbor of an encoded node side. The first neighbor (`i == 0`)
    // is the other side of the node, and the rest are the actual neighbors of the side.
    // Because we have to decompress the edges for each call, we may just as well iterate
    // over them.
    let get_neighbor = |encoded_side: u32, i: u32| -> Option<u32> {
        if i == 0 {
            // The first neighbor is the other side of the node.
            return Some(support::flip_node_side(encoded_side as usize) as u32);
        }

        let (node_index, side) = support::decode_node_side(encoded_side as usize);
        let node_id = component[node_index];
        let orientation = support::exit_orientation(side);
        if let Some(iter) = graph.successors(node_id, orientation) {
            for (j, (next_id, next_o)) in iter.enumerate() {
                if j + 1 == i as usize {
                    let next_index = node_to_index[&(next_id as u32)];
                    let next_side = support::entry_side(next_o);
                    let next = support::encode_node_side(next_index as usize, next_side) as u32;
                    return Some(next);
                }
            }
        }
        None
    };

    // Iterative Tarjan's algorithm for bridges. We only report the "internal edges" connecting the sides
    // of the same node. If there is also an external edge connecting the sides, the node is not a bridge.
    // Stack entries: (encoded side, next neighbor).
    let mut stack: Vec<(u32, u32)> = Vec::new();
    visited[0] = true;
    preorder_rank[0] = timer;
    lowest_reachable[0] = timer;
    timer += 1;
    stack.push((0, 0));

    while let Some((encoded_from, next_neighbor)) = stack.last_mut() {
        if let Some(encoded_to) = get_neighbor(*encoded_from, *next_neighbor) {
            *next_neighbor += 1;
            if !visited[encoded_to as usize] {
                visited[encoded_to as usize] = true;
                preorder_rank[encoded_to as usize] = timer;
                lowest_reachable[encoded_to as usize] = timer;
                timer += 1;
                parent[encoded_to as usize] = *encoded_from;
                stack.push((encoded_to, 0));
            } else {
                // We update `lowest_reachable` if this is a back edge. If the edge is
                // to the parent, it can still be a back edge if it's to the other side
                // of the same node using an external edge (in the actual edge list).
                let is_parent = encoded_to == parent[*encoded_from as usize];
                let is_other_side = encoded_to == support::flip_node_side(*encoded_from as usize) as u32;
                let is_external_edge = *next_neighbor > 1;
                if !is_parent || (is_other_side && is_external_edge) {
                    lowest_reachable[*encoded_from as usize] = cmp::min(
                        lowest_reachable[*encoded_from as usize], preorder_rank[encoded_to as usize]
                    );
                }
            }
        } else {
            // Done processing node from, pop and update parent.
            let encoded_from = *encoded_from;
            stack.pop();
            if let Some((encoded_parent, _)) = stack.last() {
                lowest_reachable[*encoded_parent as usize] = cmp::min(
                    lowest_reachable[*encoded_parent as usize], lowest_reachable[encoded_from as usize]
                );
                if lowest_reachable[encoded_from as usize] > preorder_rank[*encoded_parent as usize]
                    && *encoded_parent == (support::flip_node_side(encoded_from as usize) as u32) {
                        // The first condition means that this is a bridge edge. The second
                        // means that the edge is actually a bridge node in the bidirected graph.
                        let node_index = support::node_side_id(encoded_from as usize);
                        result.insert(component[node_index] as u32);
                }
            }
        }
    }

    result
}

// Returns all handles in the top-level chain with `from_handle` and `to_handle` as the tips.
// Returns `None` if there is no directed path between the tips.
fn walk_chain(
    graph: &GBZ,
    from_handle: usize, to_handle: usize,
    bridges: &HashSet<u32>
) -> Option<Vec<u32>> {
    let shortest_path = shortest_unweighted_path(graph, from_handle, to_handle)?;

    // If we visit a bridge node twice, it is a hairpin-like structure instead of
    // a node in a top-level chain.
    let mut visited_bridges = HashSet::new();
    let mut duplicate_bridges = HashSet::new();
    for &handle in shortest_path.iter() {
        let node_id = support::node_id(handle as usize) as u32;
        if bridges.contains(&node_id) {
            let first_visit = visited_bridges.insert(node_id);
            if !first_visit {
                duplicate_bridges.insert(node_id);
            }
        }
    }
    drop(visited_bridges);

    // Now determine the handles in the top-level chain and whether they are
    // entry points to a snarl.
    let mut chain_handles: Vec<u32> = Vec::new();
    for &handle in shortest_path.iter() {
        let is_tip = handle == from_handle as u32 || handle == to_handle as u32;
        let node_id = support::node_id(handle as usize) as u32;
        let is_unique_bridge = bridges.contains(&node_id) && !duplicate_bridges.contains(&node_id);
        if is_tip || is_unique_bridge {
            chain_handles.push(handle);
        }
    }

    Some(chain_handles)
}

// Finds the shortest directed path between the given handles.
// Each distance is 1 instead of the length of the node.
fn shortest_unweighted_path(graph: &GBZ, start_handle: usize, end_handle: usize) -> Option<Vec<u32>> {
    let mut queue: VecDeque<u32> = VecDeque::new();
    let mut visited: HashSet<u32> = HashSet::new();
    let mut parent: HashMap<u32, u32> = HashMap::new();

    visited.insert(start_handle as u32);
    queue.push_back(start_handle as u32);

    while let Some(handle) = queue.pop_front() {
        if handle == end_handle as u32 {
            // Reconstruct path.
            let mut path = Vec::new();
            let mut current = handle;
            while current != start_handle as u32 {
                path.push(current);
                current = parent[&current];
            }
            path.push(start_handle as u32);
            path.reverse();
            return Some(path);
        }

        let (node_id, orientation) = support::decode_node(handle as usize);
        if let Some(iter) = graph.successors(node_id, orientation) {
            for (next_id, next_o) in iter {
                let next_handle = support::encode_node(next_id, next_o) as u32;
                if !visited.contains(&next_handle) {
                    visited.insert(next_handle);
                    parent.insert(next_handle, handle);
                    queue.push_back(next_handle);
                }
            }
        }
    }

    None
}

//-----------------------------------------------------------------------------
