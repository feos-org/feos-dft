use feos_dft::python::feos_dft;
use pyo3::prelude::*;

#[pymodule]
pub fn build_wheel(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    feos_dft(py, m)
}
