//! Fill-reducing ordering: nested dissection by coordinate bisection.
//!
//! Nested dissection (A. George, *Nested dissection of a regular finite element mesh*, SIAM J.
//! Numer. Anal. 10(2):345–363, 1973) orders a separator after the two parts it splits, so the
//! Cholesky factor of a 3D mesh with n nodes fills O(n^{4/3}) instead of the O(n^{5/3}) of a band
//! or envelope ordering. The separator here is geometric: the nodes are split at the median of
//! the longest axis of their bounding box, and the nodes of the smaller side that touch the other
//! side form the separator (the coordinate nested dissection of H. D. Simon, *Partitioning of
//! unstructured problems for parallel processing*, Computing Systems in Engineering 2:135–148,
//! 1991, with vertex separators). Ties are broken by node index, so the ordering is
//! deterministic.

/// Parts this small are ordered as they come.
const LEAF: usize = 16;

/// A permutation of the nodes: `order[k]` is the node eliminated k-th.
pub fn nested_dissection(xyz: &[[f64; 3]], adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = xyz.len();
    let mut order = Vec::with_capacity(n);
    // 0 = not in the current part, 1 = left, 2 = right.
    let mut side = vec![0u8; n];
    // An explicit stack instead of recursion: (part, separator to emit after it).
    enum Task {
        Split(Vec<usize>),
        Emit(Vec<usize>),
    }
    let mut stack = vec![Task::Split((0..n).collect())];
    while let Some(task) = stack.pop() {
        let mut part = match task {
            Task::Emit(sep) => {
                order.extend(sep);
                continue;
            }
            Task::Split(p) => p,
        };
        if part.len() <= LEAF {
            order.extend(part);
            continue;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &v in &part {
            let p = xyz.get(v).copied().unwrap_or_default();
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let axis = (0..3).fold(0, |best, k| if hi[k] - lo[k] > hi[best] - lo[best] { k } else { best });
        let key = |v: usize| xyz.get(v).map_or(0.0, |p| p[axis]);
        part.sort_by(|&a, &b| key(a).total_cmp(&key(b)).then(a.cmp(&b)));
        let half = part.len() / 2;
        for (i, &v) in part.iter().enumerate() {
            if let Some(s) = side.get_mut(v) {
                *s = if i < half { 1 } else { 2 };
            }
        }
        let touches = |v: usize, other: u8| adjacency.get(v).is_some_and(|row| row.iter().any(|w| side.get(*w) == Some(&other)));
        let (left, right) = part.split_at(half);
        let left_boundary: Vec<usize> = left.iter().copied().filter(|&v| touches(v, 2)).collect();
        let right_boundary: Vec<usize> = right.iter().copied().filter(|&v| touches(v, 1)).collect();
        let separator = if left_boundary.len() <= right_boundary.len() { left_boundary } else { right_boundary };
        for &v in &part {
            if let Some(s) = side.get_mut(v) {
                *s = 0;
            }
        }
        // Mark the separator so the parts leave it out.
        for &v in &separator {
            if let Some(s) = side.get_mut(v) {
                *s = 3;
            }
        }
        let a: Vec<usize> = left.iter().copied().filter(|v| side.get(*v) != Some(&3)).collect();
        let b: Vec<usize> = right.iter().copied().filter(|v| side.get(*v) != Some(&3)).collect();
        for &v in &separator {
            if let Some(s) = side.get_mut(v) {
                *s = 0;
            }
        }
        // Order: a, then b, then the separator (the stack runs in reverse).
        stack.push(Task::Emit(separator));
        stack.push(Task::Split(b));
        stack.push(Task::Split(a));
    }
    order
}

/// A permutation of the nodes by nested dissection on the graph alone: each part's separator is
/// a middle level of the breadth-first level structure rooted at a pseudo-peripheral node, so
/// separators follow the mesh rather than a coordinate plane (A. George, J. W. H. Liu, *An
/// automatic nested dissection algorithm for irregular finite element problems*, SIAM J. Numer.
/// Anal. 15(5):1053–1069, 1978; the pseudo-peripheral node of N. E. Gibbs, W. G. Poole,
/// P. K. Stockmeyer, SIAM J. Numer. Anal. 13(2):236–250, 1976). A thin curved wall (a shell, a
/// bellows) is cut across its surface instead of many times by a plane. Disconnected parts are
/// ordered one after the other; ties go by node index, so the ordering is deterministic.
pub fn graph_dissection(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let mut order = Vec::with_capacity(n);
    // The part each node belongs to (a token per split), and its level in the current search.
    let mut owner = vec![0usize; n];
    let mut level = vec![usize::MAX; n];
    let mut token = 1usize;
    enum Task {
        Split(Vec<usize>),
        Emit(Vec<usize>),
    }
    let mut stack = vec![Task::Split((0..n).collect())];
    // The level structure of the nodes of part `token` reached from `root`.
    let levels = |root: usize, token: usize, owner: &[usize], level: &mut Vec<usize>| -> Vec<Vec<usize>> {
        let mut out: Vec<Vec<usize>> = vec![vec![root]];
        if let Some(l) = level.get_mut(root) {
            *l = 0;
        }
        loop {
            let depth = out.len();
            let mut next = Vec::new();
            for &v in out.last().into_iter().flatten() {
                for &w in adjacency.get(v).into_iter().flatten() {
                    if owner.get(w) == Some(&token) && level.get(w) == Some(&usize::MAX) {
                        if let Some(l) = level.get_mut(w) {
                            *l = depth;
                        }
                        next.push(w);
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            out.push(next);
        }
        out
    };
    let degree = |v: usize, token: usize, owner: &[usize]| adjacency.get(v).map_or(0, |row| row.iter().filter(|w| owner.get(**w) == Some(&token)).count());
    while let Some(task) = stack.pop() {
        let part = match task {
            Task::Emit(sep) => {
                order.extend(sep);
                continue;
            }
            Task::Split(p) => p,
        };
        if part.len() <= LEAF {
            order.extend(part);
            continue;
        }
        token += 1;
        for &v in &part {
            if let Some(o) = owner.get_mut(v) {
                *o = token;
            }
        }
        let reset = |level: &mut Vec<usize>, nodes: &[usize]| {
            for &v in nodes {
                if let Some(l) = level.get_mut(v) {
                    *l = usize::MAX;
                }
            }
        };
        // A pseudo-peripheral root: from the first node, the end of the deepest level structure
        // with the fewest neighbours, while the depth grows.
        let mut root = part.first().copied().unwrap_or(0);
        let mut structure = levels(root, token, &owner, &mut level);
        for _ in 0..PERIPHERAL_SEARCHES {
            let reached: Vec<usize> = structure.iter().flatten().copied().collect();
            reset(&mut level, &reached);
            let Some(candidate) = structure.last().and_then(|last| last.iter().copied().min_by_key(|v| (degree(*v, token, &owner), *v))) else { break };
            let deeper = levels(candidate, token, &owner, &mut level);
            if deeper.len() <= structure.len() {
                let reached: Vec<usize> = deeper.iter().flatten().copied().collect();
                reset(&mut level, &reached);
                structure = levels(root, token, &owner, &mut level);
                break;
            }
            (root, structure) = (candidate, deeper);
        }
        let reached: Vec<usize> = structure.iter().flatten().copied().collect();
        reset(&mut level, &reached);
        if reached.len() < part.len() {
            // Disconnected: this component, then the rest.
            let mut inside = vec![false; 0];
            inside.resize(n, false);
            for &v in &reached {
                if let Some(x) = inside.get_mut(v) {
                    *x = true;
                }
            }
            let rest: Vec<usize> = part.iter().copied().filter(|v| !inside.get(*v).copied().unwrap_or(false)).collect();
            let mut component = reached;
            component.sort_unstable();
            stack.push(Task::Split(rest));
            stack.push(Task::Split(component));
            continue;
        }
        if structure.len() < 3 {
            order.extend(part);
            continue;
        }
        // The middle level: the first whose levels up to it hold half the part.
        let mut count = 0;
        let mut middle = 1;
        for (i, l) in structure.iter().enumerate().skip(1) {
            count += l.len();
            if 2 * count >= part.len() {
                middle = i.min(structure.len() - 2);
                break;
            }
        }
        // Its nodes that touch the next level separate the levels before it from those after.
        for (i, l) in structure.iter().enumerate() {
            for &v in l {
                if let Some(x) = level.get_mut(v) {
                    *x = i;
                }
            }
        }
        let touches_next = |v: usize, level: &[usize]| adjacency.get(v).is_some_and(|row| row.iter().any(|w| owner.get(*w) == Some(&token) && level.get(*w) == Some(&(middle + 1))));
        let mut separator = Vec::new();
        let mut before: Vec<usize> = structure.iter().take(middle).flatten().copied().collect();
        for &v in structure.get(middle).into_iter().flatten() {
            if touches_next(v, &level) { separator.push(v) } else { before.push(v) }
        }
        let after: Vec<usize> = structure.iter().skip(middle + 1).flatten().copied().collect();
        reset(&mut level, &part);
        before.sort_unstable();
        separator.sort_unstable();
        let mut after = after;
        after.sort_unstable();
        stack.push(Task::Emit(separator));
        stack.push(Task::Split(after));
        stack.push(Task::Split(before));
    }
    order
}

/// Searches for a pseudo-peripheral node.
const PERIPHERAL_SEARCHES: usize = 8;

/// The entries of the Cholesky factor of the node graph `adjacency` eliminated in `order` (one
/// per node pair): the fill an ordering gives.
pub fn fill(adjacency: &[Vec<usize>], order: &[usize]) -> usize {
    let n = adjacency.len();
    let mut position = vec![0usize; n];
    for (k, &v) in order.iter().enumerate() {
        if let Some(p) = position.get_mut(v) {
            *p = k;
        }
    }
    let mut ptr = vec![0usize];
    let mut idx = Vec::new();
    for &v in order {
        let k = position.get(v).copied().unwrap_or(0);
        let mut rows: Vec<usize> = adjacency.get(v).into_iter().flatten().filter_map(|w| position.get(*w).copied()).filter(|r| *r <= k).collect();
        rows.push(k);
        rows.sort_unstable();
        rows.dedup();
        idx.extend(rows);
        ptr.push(idx.len());
    }
    let val = vec![1.0; idx.len()];
    crate::cholesky::factor_entries(&crate::cholesky::UpperColumns { n, ptr, idx, val })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_permutation_with_separators_last() {
        // A 12 × 12 grid graph.
        let m = 12;
        let id = |i: usize, j: usize| i * m + j;
        let xyz: Vec<[f64; 3]> = (0..m * m).map(|v| [(v / m) as f64, (v % m) as f64, 0.0]).collect();
        let mut adj = vec![Vec::new(); m * m];
        for i in 0..m {
            for j in 0..m {
                for (di, dj) in [(0i64, 1i64), (1, 0), (0, -1), (-1, 0), (0, 0)] {
                    let (a, b) = (i as i64 + di, j as i64 + dj);
                    if (0..m as i64).contains(&a) && (0..m as i64).contains(&b) {
                        adj[id(i, j)].push(id(a as usize, b as usize));
                    }
                }
            }
        }
        let order = nested_dissection(&xyz, &adj);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..m * m).collect::<Vec<_>>());
        // The top separator (last 12 entries) is one grid line.
        let last: Vec<usize> = order[m * m - m..].to_vec();
        assert!(last.iter().all(|v| v / m == last[0] / m) || last.iter().all(|v| v % m == last[0] % m), "{last:?}");
        assert_eq!(order, nested_dissection(&xyz, &adj), "deterministic");
    }
}
