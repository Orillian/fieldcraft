//! Initial fields: the state a structural analysis starts from, given at the nodes of the mesh,
//! usually mapped from another analysis's results (a moulding or casting simulation, a thermal
//! analysis, a coarser or finer mesh of the part):
//! - an initial **temperature**: the temperature each node starts at, from which thermal strains
//!   are measured, ε_th = α (T − T₀) (a node it leaves out starts at its material's reference
//!   temperature); a static step's temperature loads then strain the part by the change, so a
//!   part ejected hot and cooled to the room temperature warps by its own temperature field;
//! - an initial **stress** (xx, yy, zz, xy, xz, yz): a residual stress the part carries before
//!   the step, σ = D (ε − ε_th) + σ₀; what of it is not in equilibrium with the supports is
//!   released into deformation (warpage after ejection) and the stress that remains.
//!
//! The fields act in every linear static and coupled temperature–displacement step ([`acts_in`]). They
//! are tied to node labels: a remeshed part keeps them only where the labels survive, so a field
//! is mapped again after remeshing.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, MAX_NODES, Model, Procedure, Result, Step, context, insert_new, remove_existing};

/// Whether initial fields act in a step of this procedure: linear static and coupled
/// temperature–displacement steps (other procedures leave them out and say so).
pub fn acts_in(procedure: &Procedure) -> bool {
    matches!(procedure, Procedure::Static {} | Procedure::CoupledTemperatureDisplacement {})
}

/// The warning for a step that leaves the model's initial fields out (every solver gives it,
/// the explicit one included); `None` when the step takes them or there are none.
pub fn left_out(model: &Model, step: &Step) -> Option<String> {
    (!model.initial_fields().is_empty() && !acts_in(&step.procedure)).then(|| {
        format!(
            "the initial fields act in static and coupled temperature–displacement steps only; step `{}` starts without them (at rest, at the reference temperature)",
            step.name
        )
    })
}

/// Longest description of where a field came from.
pub const MAX_SOURCE_LEN: usize = 1024;

/// An initial field: one value per node, by node label.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum InitialField {
    /// The temperature each node starts at.
    Temperature {
        #[serde(with = "by_node")]
        values: BTreeMap<u32, f64>,
        /// Where the values came from (a file and step), for the person reading the model.
        #[serde(default, skip_serializing_if = "String::is_empty")]
        source: String,
    },
    /// The stress at each node, xx, yy, zz, xy, xz, yz.
    Stress {
        #[serde(with = "by_node")]
        values: BTreeMap<u32, [f64; 6]>,
        #[serde(default, skip_serializing_if = "String::is_empty")]
        source: String,
    },
}

/// Values by node, written as `[[label, value], …]` in label order (a map's keys would be
/// strings, which a tagged enum cannot read back as labels); a label listed twice is refused.
mod by_node {
    use std::collections::BTreeMap;

    use serde::de::Error as _;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<V: Serialize, S: Serializer>(values: &BTreeMap<u32, V>, s: S) -> Result<S::Ok, S::Error> {
        s.collect_seq(values.iter())
    }

    pub fn deserialize<'de, V: Deserialize<'de>, D: Deserializer<'de>>(d: D) -> Result<BTreeMap<u32, V>, D::Error> {
        let pairs = Vec::<(u32, V)>::deserialize(d)?;
        let mut out = BTreeMap::new();
        for (node, v) in pairs {
            if out.insert(node, v).is_some() {
                return Err(D::Error::custom(format!("node {node} is listed twice")));
            }
        }
        Ok(out)
    }
}

impl InitialField {
    /// `temperature` or `stress`.
    pub fn kind(&self) -> &'static str {
        match self {
            InitialField::Temperature { .. } => "temperature",
            InitialField::Stress { .. } => "stress",
        }
    }

    /// The number of nodes it gives a value.
    pub fn len(&self) -> usize {
        match self {
            InitialField::Temperature { values, .. } => values.len(),
            InitialField::Stress { values, .. } => values.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Where its values came from.
    pub fn source(&self) -> &str {
        match self {
            InitialField::Temperature { source, .. } | InitialField::Stress { source, .. } => source,
        }
    }

    /// The node labels it gives a value, ascending.
    pub fn nodes(&self) -> Vec<u32> {
        match self {
            InitialField::Temperature { values, .. } => values.keys().copied().collect(),
            InitialField::Stress { values, .. } => values.keys().copied().collect(),
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.len() > MAX_NODES {
            return Err(Error::Limit(format!("an initial field gives at most {MAX_NODES} nodes a value")));
        }
        if self.source().len() > MAX_SOURCE_LEN {
            return Err(Error::Invalid(format!("the source description is longer than {MAX_SOURCE_LEN} bytes")));
        }
        let bad = match self {
            InitialField::Temperature { values, .. } => values.iter().find(|(_, v)| !v.is_finite()).map(|(n, _)| *n),
            InitialField::Stress { values, .. } => values.iter().find(|(_, v)| v.iter().any(|c| !c.is_finite())).map(|(n, _)| *n),
        };
        if let Some(node) = bad {
            return Err(Error::Invalid(format!("node {node}: the {} must be finite", self.kind())));
        }
        Ok(())
    }
}

impl Model {
    /// The initial fields, by name.
    pub fn initial_fields(&self) -> &BTreeMap<String, InitialField> {
        &self.initial_fields
    }

    pub fn initial_field(&self, name: &str) -> Option<&InitialField> {
        self.initial_fields.get(name)
    }

    /// Add an initial field, or replace the one of that name (mapping it again).
    pub fn set_initial_field(&mut self, name: &str, field: InitialField) -> Result<()> {
        field.validate().map_err(|e| context("initial field", name, e))?;
        match self.initial_fields.get_mut(name) {
            Some(slot) => {
                *slot = field;
                Ok(())
            }
            None => insert_new(&mut self.initial_fields, "initial field", name, field),
        }
    }

    pub fn remove_initial_field(&mut self, name: &str) -> Result<InitialField> {
        remove_existing(&mut self.initial_fields, "initial field", name)
    }

    /// The initial temperature of every node one is given, and the initial stress likewise;
    /// fails when two fields give one node different values. Labels the mesh no longer has are
    /// kept (the solver skips them).
    pub fn initial_state(&self) -> Result<(BTreeMap<u32, f64>, BTreeMap<u32, [f64; 6]>)> {
        let (mut temperature, mut stress) = (BTreeMap::new(), BTreeMap::new());
        fn merge<V: Copy + PartialEq + std::fmt::Debug>(out: &mut BTreeMap<u32, V>, name: &str, values: &BTreeMap<u32, V>) -> Result<()> {
            for (node, v) in values {
                match out.insert(*node, *v) {
                    Some(old) if old != *v => {
                        return Err(Error::Invalid(format!(
                            "node {node} is given both {old:?} and {v:?} (initial field `{name}` and another of its kind)"
                        )));
                    }
                    _ => {}
                }
            }
            Ok(())
        }
        for (name, f) in &self.initial_fields {
            match f {
                InitialField::Temperature { values, .. } => merge(&mut temperature, name, values)?,
                InitialField::Stress { values, .. } => merge(&mut stress, name, values)?,
            }
        }
        Ok((temperature, stress))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temperature(values: &[(u32, f64)]) -> InitialField {
        InitialField::Temperature { values: values.iter().copied().collect(), source: String::new() }
    }

    #[test]
    fn fields_validate_round_trip_and_merge() {
        let mut m = Model::default();
        m.set_initial_field("hot", temperature(&[(1, 400.0), (2, 410.0)])).unwrap();
        let stress = InitialField::Stress { values: [(2, [1e6, 0.0, 0.0, 0.0, 0.0, -2e5])].into(), source: "mould.fieldcraft, step Cool".into() };
        m.set_initial_field("residual", stress.clone()).unwrap();
        assert_eq!(m.initial_field("residual"), Some(&stress));
        // Replacing maps it again.
        m.set_initial_field("hot", temperature(&[(1, 400.0), (2, 420.0)])).unwrap();
        let (t, s) = m.initial_state().unwrap();
        assert_eq!((t.len(), s.len(), t[&2]), (2, 1, 420.0));
        // Two fields agreeing at a shared node are fine; disagreeing ones are refused.
        m.set_initial_field("more", temperature(&[(2, 420.0), (3, 300.0)])).unwrap();
        assert_eq!(m.initial_state().unwrap().0.len(), 3);
        m.set_initial_field("clash", temperature(&[(1, 0.0)])).unwrap();
        assert!(m.initial_state().unwrap_err().to_string().contains("node 1"));
        m.remove_initial_field("clash").unwrap();
        assert!(m.remove_initial_field("clash").is_err());
        // Through JSON and back.
        let json = serde_json::to_string(&m).unwrap();
        let back: Model = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);
        // A model without fields writes nothing for them.
        assert!(!serde_json::to_string(&Model::default()).unwrap().contains("initial_fields"));
    }

    #[test]
    fn the_check_names_clashes_and_nodes_the_mesh_lacks() {
        let mut m = Model::default();
        m.set_node(1, [0.0; 3]).unwrap();
        m.set_initial_field("hot", temperature(&[(1, 400.0), (7, 400.0)])).unwrap();
        m.set_initial_field("cold", temperature(&[(1, 300.0)])).unwrap();
        let found = m.check();
        let said = |what: &str| found.iter().any(|d| format!("{} {}", d.subject, d.message).contains(what));
        assert!(said("initial field `hot` 1 of its 2 nodes are not in the mesh"), "{found:?}");
        assert!(said("node 1 is given both"), "{found:?}");
    }

    #[test]
    fn hostile_fields_are_refused() {
        let mut m = Model::default();
        assert!(m.set_initial_field("bad", temperature(&[(1, f64::NAN)])).is_err());
        assert!(
            m.set_initial_field(
                "bad",
                InitialField::Stress { values: [(1, [0.0, f64::INFINITY, 0.0, 0.0, 0.0, 0.0])].into(), source: String::new() }
            )
            .is_err()
        );
        assert!(m.set_initial_field("bad", InitialField::Temperature { values: BTreeMap::new(), source: "x".repeat(MAX_SOURCE_LEN + 1) }).is_err());
        assert!(m.set_initial_field("1 bad", temperature(&[])).is_err());
        assert!(serde_json::from_str::<Model>(r#"{"initial_fields": {"f": {"type": "temperature", "values": [[1, 1.0]], "extra": 1}}}"#).is_err());
        assert!(serde_json::from_str::<Model>(r#"{"initial_fields": {"f": {"type": "stress", "values": [[1, [1.0]]]}}}"#).is_err());
        assert!(serde_json::from_str::<Model>(r#"{"initial_fields": {"f": {"type": "temperature", "values": [[1, 1.0], [1, 2.0]]}}}"#).is_err());
        assert!(serde_json::from_str::<Model>(r#"{"initial_fields": {"f": {"type": "temperature", "values": [[-1, 1.0]]}}}"#).is_err());
        assert!(m.initial_fields().is_empty());
    }
}
