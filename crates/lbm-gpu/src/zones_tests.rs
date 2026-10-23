//! Upload checks that need no GPU device.
#[path = "../../lbm/tests/common/mod.rs"]
#[allow(clippy::unwrap_used)]
mod common;
use super::zones;
use common::*;
use fieldcraft_model::{BoxFace, FaceCondition, Fan, FanDirection};

#[test]
fn curves_that_collapse_in_single_precision_are_refused_before_device_work() {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    d.fans = vec![Fan {
        name: "Indistinguishable".into(),
        min: [2.0, 0.0, 0.0],
        max: [4.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.002], [0.1, 0.001], [0.1 + 1e-10, 0.0]],
        relaxation: None,
    }];
    let m = model(1.0, 0.1, d, step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    let mut flags = vec![1; flow.grid().len()];
    let e = zones::append(&flow, &mut flags).unwrap_err().to_string();
    assert!(e.contains("Indistinguishable") && e.contains("single precision"), "{e}");
}

#[test]
fn zones_leave_the_flags_alone_and_offsets_are_absolute() {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    let m = model(1.0, 0.1, d.clone(), step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    let mut flags = vec![1; 32];
    assert_eq!(zones::append(&flow, &mut flags).unwrap(), zones::Layout::default());
    assert_eq!(flags, vec![1; 32]);
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [2.0, 0.0, 0.0],
        max: [4.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.002], [0.1, 0.0]],
        relaxation: None,
    }];
    let m = model(1.0, 0.1, d, step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    // Flag words a wall model would have set stay as they are: ownership is its own field.
    let mut flags: Vec<u32> = (0..32).map(|c| 1 | (c as u32) << 1).collect();
    let before = flags.clone();
    let layout = zones::append(&flow, &mut flags).unwrap();
    assert_eq!(flags[..32], before[..]);
    let owned: Vec<u32> = (0..32u32).filter(|c| flags[32 + *c as usize] == 1).collect();
    assert_eq!(owned.len(), 8);
    assert!(flags[32..64].iter().all(|z| *z <= 1));
    let header = 64;
    assert_eq!(flags[header + 13] as usize, header + zones::WORDS);
    let curve = flags[header + 13] as usize;
    assert_eq!(f32::from_bits(flags[curve]), 0.0);
    assert_eq!(f32::from_bits(flags[curve + 1]), 0.002f32);
    assert_eq!(layout, zones::Layout { zones: 1, fans: 1, fan_table: (curve + 4) as u32, fan_groups: 1 });
    assert_eq!(flags[layout.fan_table as usize], 0);
    assert_eq!((flags[header + 16], flags[header + 17]), (0, 1));
    let list = &flags[layout.fan_table as usize + 1..];
    assert_eq!(list.len(), zones::WG);
    assert_eq!(list[..8], owned[..]);
    assert!(list[8..].iter().all(|c| *c == zones::NONE));
}
