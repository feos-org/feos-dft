#[macro_export]
macro_rules! impl_pair_correlation {
    ($func:ty) => {
        /// A one-dimensional initial density profile of a fluid in a pore.
        ///
        /// Parameters
        /// ----------
        /// bulk : State
        ///     The bulk state in equilibrium with the profile.
        /// n_grid : int
        ///     The number of grid points.
        /// width: SINumber
        ///     The width of the system.
        ///
        /// Returns
        /// -------
        /// PairCorrelation
        ///
        #[pyclass(name = "PairCorrelation", unsendable)]
        #[pyo3(text_signature = "(bulk, n_grid, width)")]
        pub struct PyPairCorrelation(PairCorrelation<SIUnit, $func>);

        impl_1d_profile!(PyPairCorrelation, [get_r]);

        #[pymethods]
        impl PyPairCorrelation {
            #[new]
            fn new(bulk: PyState, n_grid: usize, width: PySINumber) -> PyResult<Self> {
                Ok(PyPairCorrelation(PairCorrelation::new(
                    &bulk.0,
                    n_grid,
                    width.into(),
                )?))
            }
        }

        #[pymethods]
        impl PyPairCorrelation {
            #[getter]
            fn get_pair_correlation_function<'py>(
                &self,
                py: Python<'py>,
            ) -> Option<&'py PyArray2<f64>> {
                self.0
                    .pair_correlation_function
                    .as_ref()
                    .map(|g| g.view().to_pyarray(py))
            }

            #[getter]
            fn get_self_solvation_free_energy(&self) -> Option<PySINumber> {
                self.0.self_solvation_free_energy.map(PySINumber::from)
            }

            #[getter]
            fn get_structure_factor(&self) -> Option<f64> {
                self.0.structure_factor
            }
        }
    };
}
