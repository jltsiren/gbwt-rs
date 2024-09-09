//! Algorithms using GBWT and GBZ.

use crate::{GBZ, support};

use std::collections::BTreeMap;
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
    for i in 0..a.len() {
        for j in 0..b.len() {
            dp[i + 1][j + 1] = cmp::max(dp[i + 1][j], dp[i][j + 1]);
            if a[i] == b[j] {
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
                return Some(point.clone());
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
        if prev.is_some() && next.is_some() {
            let prev = prev.unwrap();
            let next = next.unwrap();
            if prev.weight > next.weight {
                Some((prev, edits - self.a_weight(a_offset - 1)))
            } else {
                Some((next, edits - self.b_weight(b_offset - 1)))
            }
        } else if let Some(point) = prev {
            Some((point, edits - self.a_weight(a_offset - 1)))
        } else if let Some(point) = next {
            Some((point, edits - self.b_weight(b_offset - 1)))
        } else {
            None
        }
    }
}

// FIXME tests against the naive algorithm
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
/// use gbwt::algorithms::lcs;
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
/// use gbwt::{GBZ, Orientation, algorithms};
/// use gbwt::support;
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
