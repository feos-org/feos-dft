use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::PyInit_quantity;

mod adsorption;
mod fundamental_measure_theory;
mod interface;
mod profile;
mod solvation;
mod solver;

pub use adsorption::{PyExternalPotential, PyGeometry};
use fundamental_measure_theory::*;
pub use solver::PyDFTSolver;

#[pymodule]
pub fn feos_dft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyGeometry>()?;
    m.add_class::<PyDFTSolver>()?;

    m.add_class::<PyFMTVersion>()?;
    m.add_class::<PyFMTFunctional>()?;

    m.add_class::<PyState>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyPairCorrelation>()?;
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;

    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_dft.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
