//! Python binding for batch entity decisions.
use pyo3::prelude::*;
mod precomputed;

#[pymodule]
fn _omnipath_resolver(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(precomputed::resolve_precomputed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(precomputed::resolve_molecular_batch, m)?)?;
    Ok(())
}
