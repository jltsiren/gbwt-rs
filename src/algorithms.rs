//! Algorithms using GBWT and GBZ.

use crate::{GBZ, support};

use std::cmp;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

// FIXME replace with the O(nd) algorithm
/// Returns the longest common subsequence of integer sequences `a` and `b`, weighted by the given function.
///
/// The subsequence is returned as pairs of positions, and the second return value is the total weight of the LCS.
/// See [`lcs`] for an unweighted version and [`path_lcs`] for finding the LCS of two paths in a graph.
pub fn weighted_lcs<F: Fn(usize) -> usize>(a: &[usize], b: &[usize], weight: F) -> (Vec<(usize, usize)>, usize) {
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
    weighted_lcs(a, b, |_| 1).0
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
    weighted_lcs(a, b, weight)
}

//-----------------------------------------------------------------------------
