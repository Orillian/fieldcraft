//! Fan and porous data appended to the cell-flags buffer, so no kernel binds another buffer.
//!
//! The flags keep their own layout (fluid bit, wall-model row or on-node mark; `native.rs`);
//! zone ownership is a separate field. With zones the buffer holds, in words:
//!
//! - `[0, n)`: the flags, unchanged;
//! - `[n, 2n)`: each cell's zone plus one (0 outside every zone);
//! - [`WORDS`] per zone: its header (layout below);
//! - each fan's curve, its points as pairs of f32 bits;
//! - the fan table: each fan's zone index, in zone order;
//! - each fan's fluid cells, padded with [`NONE`] to whole workgroups of [`WG`], so a workgroup
//!   of `fan_partials` sums the cells of one fan only.
//!
//! Header words: 0 fan (1) or porous (0); 1 axis; 2–4 lo and 5–7 hi (cell indices); 8 α and
//! 9 β (porous resistance, f32 bits); 10 the fan's direction sign; 11 its relaxation; 12 its
//! current pressure rise (the only word a kernel writes: `fan_update`); 13 its curve's first
//! word and 14 its point count; 15 the layers along the axis; 16 its first workgroup in the cell
//! lists (which is also its first slot in `partials`) and 17 its number of workgroups; 18–19
//! unused. All offsets are absolute word indices.

use fieldcraft_lbm::{Error, Flow, Result, ZoneKind};

pub(super) const WORDS: usize = fieldcraft_lbm::memory::DEVICE_ZONE_WORDS;

/// Cells per workgroup (`WG` in lbm.wgsl).
pub(super) const WG: usize = 256;

/// A padding entry of a fan's cell list.
pub(super) const NONE: u32 = u32::MAX;

/// Where the zones are in the flags buffer, for the kernels' parameters and dispatches.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Layout {
    /// Fan and porous zones.
    pub zones: u32,
    /// Fans among them.
    pub fans: u32,
    /// The fan table's first word; the cell lists follow it.
    pub fan_table: u32,
    /// Workgroups over every fan's padded cell list (slots of `partials` that `fan_partials`
    /// fills).
    pub fan_groups: u32,
}

fn index(n: usize) -> Result<u32> {
    u32::try_from(n).ok().filter(|v| *v != NONE).ok_or_else(|| Error::Limit("too much fan/porous data for the GPU".into()))
}

fn value(name: &str, label: &str, v: f64) -> Result<u32> {
    let stored = v as f32;
    if !stored.is_finite() || (v != 0.0 && stored == 0.0) {
        return Err(Error::Invalid(format!("zone `{name}`: {label} cannot be represented in GPU single precision")));
    }
    Ok(stored.to_bits())
}

/// Add the zones of `flow` to `flags` (one per cell, prepared by the caller), keeping any fan
/// rise a warmed-up CPU flow has reached. Refuses curves that collapse in single precision before
/// any device work. Leaves `flags` unchanged without zones.
pub(super) fn append(flow: &Flow, flags: &mut Vec<u32>) -> Result<Layout> {
    let zones = flow.zones();
    if zones.is_empty() {
        return Ok(Layout::default());
    }
    let n = flow.grid().len();
    if flags.len() != n {
        return Err(Error::Invalid("GPU zone upload needs one flag per cell".into()));
    }
    let cells = flow.zone_cells();
    let too_much = || Error::Limit("too much fan/porous data for the GPU".into());
    let fans: Vec<usize> = zones.iter().enumerate().filter(|(_, z)| z.is_fan()).map(|(k, _)| k).collect();
    let mut extra = n.checked_add(zones.len().checked_mul(WORDS).ok_or_else(too_much)?).ok_or_else(too_much)?;
    for &k in &fans {
        let curve = match zones.get(k).map(|z| &z.kind) {
            Some(ZoneKind::Fan { curve, .. }) => curve.len(),
            _ => 0,
        };
        let list = cells.get(k).map_or(0, |c| c.len().div_ceil(WG) * WG);
        extra = extra.checked_add(curve.checked_mul(2).ok_or_else(too_much)?.checked_add(1 + list).ok_or_else(too_much)?).ok_or_else(too_much)?;
    }
    index(n.checked_add(extra).ok_or_else(too_much)?)?;
    flags.try_reserve_exact(extra).map_err(|_| Error::Limit("not enough memory for GPU zones".into()))?;
    // Zone ownership, then room for the headers.
    flags.resize(2 * n, 0);
    for (k, list) in cells.iter().enumerate() {
        let name = zones.get(k).map_or("", |z| z.name.as_str());
        for &cell in list {
            let slot = flags
                .get_mut(n + cell as usize)
                .filter(|_| (cell as usize) < n)
                .ok_or_else(|| Error::Invalid(format!("zone `{name}` has a cell outside the grid")))?;
            *slot = index(k + 1)?;
        }
    }
    let headers = flags.len();
    flags.resize(headers + zones.len() * WORDS, 0);
    let samples = flow.zone_samples();
    let mut data = vec![[0u32; WORDS]; zones.len()];
    for (k, (zone, d)) in zones.iter().zip(data.iter_mut()).enumerate() {
        d[1] = index(zone.axis)?;
        for axis in 0..3 {
            d[2 + axis] = index(zone.lo[axis])?;
            d[5 + axis] = index(zone.hi[axis])?;
        }
        d[10] = 1.0f32.to_bits();
        d[15] = index(zone.layers())?;
        match &zone.kind {
            ZoneKind::Porous { alpha, beta } => {
                d[8] = value(&zone.name, "Darcy resistance", *alpha)?;
                d[9] = value(&zone.name, "inertial resistance", *beta)?;
            }
            ZoneKind::Fan { sign, curve, relaxation } => {
                d[0] = 1;
                d[10] = value(&zone.name, "direction", *sign)?;
                d[11] = value(&zone.name, "relaxation", *relaxation)?;
                d[12] = value(&zone.name, "current pressure rise", samples.get(k).map_or(0.0, |s| s.fan_rise))?;
                d[13] = index(flags.len())?;
                d[14] = index(curve.len())?;
                let mut previous = None;
                for point in curve {
                    let x = value(&zone.name, "curve flow", point[0])?;
                    let y = value(&zone.name, "curve rise", point[1])?;
                    let at = f32::from_bits(x);
                    if previous.is_some_and(|p| at <= p) {
                        return Err(Error::Invalid(format!("zone `{}`: fan curve flows are indistinguishable in GPU single precision", zone.name)));
                    }
                    previous = Some(at);
                    flags.extend([x, y]);
                }
            }
        }
    }
    let fan_table = index(flags.len())?;
    for &k in &fans {
        flags.push(index(k)?);
    }
    let mut group = 0usize;
    for &k in &fans {
        let list = cells.get(k).map_or(&[][..], |c| &c[..]);
        let groups = list.len().div_ceil(WG);
        if let Some(d) = data.get_mut(k) {
            d[16] = index(group)?;
            d[17] = index(groups)?;
        }
        flags.extend_from_slice(list);
        flags.resize(flags.len() + groups * WG - list.len(), NONE);
        group += groups;
    }
    for (k, d) in data.iter().enumerate() {
        let start = headers + k * WORDS;
        flags.get_mut(start..start + WORDS).ok_or_else(|| Error::Invalid("missing GPU zone header".into()))?.copy_from_slice(d);
    }
    Ok(Layout { zones: index(zones.len())?, fans: index(fans.len())?, fan_table, fan_groups: index(group)? })
}
