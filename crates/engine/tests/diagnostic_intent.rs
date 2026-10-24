//! The engine as a dependency (not cfg(test)): a session no host gave a folder records no solver
//! diagnostic, wherever the environment points the crash folder.
use fieldcraft_engine::{EngineError, Session};
use serde_json::json;

#[test]
fn child_embedded_session() {
    if std::env::var_os("FIELDCRAFT_DIAGNOSTIC_INTENT_CHILD").is_none() {
        return;
    }
    let mut session = Session::default();
    session.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
    // Without its clamp the static solve fails (singular): a solver failure.
    session.execute("model.step.edit", &json!({"name": "Load", "boundary_conditions": []})).unwrap();
    let error = session.execute("solve.static", &json!({"step": "Load"})).unwrap_err();
    assert!(matches!(error, EngineError::Solver(fieldcraft_solver::Error::Singular(_))), "{error}");
    assert_eq!(session.take_solver_diagnostic(), None);
}

#[test]
fn a_dependency_session_does_not_write_reports_without_host_intent() {
    let directory = std::env::temp_dir().join(format!("fieldcraft-diagnostic-intent-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "child_embedded_session", "--nocapture"])
        .env("FIELDCRAFT_DIAGNOSTIC_INTENT_CHILD", "1")
        .env("LOCALAPPDATA", &directory)
        .env("HOME", &directory)
        .env("XDG_DATA_HOME", &directory)
        .env("FIELDCRAFT_CRASH_DIR", directory.join("crashes"))
        .env_remove("CI")
        .output()
        .unwrap();
    assert!(output.status.success(), "{}\n{}", String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 0);
    let _ = std::fs::remove_dir_all(directory);
}
