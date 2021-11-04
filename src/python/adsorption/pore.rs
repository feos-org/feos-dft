#[macro_export]
macro_rules! impl_pore {
    ($func:ty, $py_func:ty) => {
        /// A one-dimensional density profile of a fluid in a pore.
        ///
        /// Parameters
        /// ----------
        /// functional : $func
        ///     The Helmholtz energy functional.
        /// geometry : Geometry
        ///     The pore geometry.
        /// pore_size : SINumber
        ///     The width of the slit pore.
        /// potential : ExternalPotential
        ///     The potential used to model wall-fluid interactions.
        /// n_grid : int, optional
        ///     The number of grid points.
        /// potential_cutoff : float, optional
        ///     Maximum value for the external potential.
        ///
        /// Returns
        /// -------
        /// Pore1D
        ///
        #[pyclass(name = "Pore1D", unsendable)]
        #[pyo3(text_signature = "(functional, geometry, pore_size, potential, n_grid=None, potential_cutoff=None)")]
        pub struct PyPore1D(Pore1D<SIUnit, $func>);

        #[pyclass(name = "PoreProfile1D", unsendable)]
        pub struct PyPoreProfile1D(PoreProfile1D<SIUnit, $func>);

        impl_1d_profile!(PyPoreProfile1D, [get_r, get_z]);

        #[pymethods]
        impl PyPore1D {
            #[new]
            fn new(
                functional: &$py_func,
                geometry: PyGeometry,
                pore_size: PySINumber,
                potential: PyExternalPotential,
                n_grid: Option<usize>,
                potential_cutoff: Option<f64>,
            ) -> Self {
                Self(Pore1D::new(
                    &functional.0,
                    geometry.0,
                    pore_size.into(),
                    potential.0,
                    n_grid,
                    potential_cutoff,
                ))
            }

            /// Initialize the pore for the given bulk state
            ///
            /// Parameters
            /// ----------
            /// bulk : State
            ///     The bulk state in equilibrium with the pore.
            /// external_potential : Array2, optional
            ///     The external potential in the pore. Used to
            ///     save computation time in the case of costly
            ///     evaluations of external potentials.
            ///
            /// Returns
            /// -------
            /// PoreProfile1D
            #[pyo3(text_signature = "(bulk, external_potential=None)")]
            fn initialize(
                &self,
                bulk: &PyState,
                external_potential: Option<&PyArray2<f64>>,
            ) -> PyResult<PyPoreProfile1D> {
                Ok(PyPoreProfile1D(self.0.initialize(
                    &bulk.0,
                    external_potential.map(|e| e.to_owned_array()).as_ref(),
                )?))
            }
        }

        #[pymethods]
        impl PyPoreProfile1D {
            #[getter]
            fn get_grand_potential(&self) -> Option<PySINumber> {
                self.0.grand_potential.map(PySINumber::from)
            }

            #[getter]
            fn get_interfacial_tension(&self) -> Option<PySINumber> {
                self.0.interfacial_tension.map(PySINumber::from)
            }
        }

        /// A three-dimensional density profile of a fluid in a pore.
        ///
        /// Parameters
        /// ----------
        /// functional : $func
        ///     The Helmholtz energy functional.
        /// system_size : [SINumber; 3]
        ///     The size of the unit cell.
        /// n_grid : [int; 3]
        ///     The number of grid points in each direction.
        /// coordinates : SIArray2
        ///     The positions of all interaction sites in the solid.
        /// sigma_ss : Array1
        ///     The size parameters of all interaction sites.
        /// epsilon_k_ss : Array1
        ///     The energy parameter of all interaction sites.
        /// potential_cutoff: float, optional
        ///     Maximum value for the external potential.
        /// cutoff_radius: SINumber, optional
        ///     The cutoff radius for the calculation of solid-fluid interactions.
        ///
        /// Returns
        /// -------
        /// Pore3D
        ///
        #[pyclass(name = "Pore3D", unsendable)]
        #[pyo3(text_signature = "(functional, system_size, n_grid, coordinates, sigma_ss, epsilon_k_ss, potential_cutoff=None, cutoff_radius=None)")]
        pub struct PyPore3D(Pore3D<SIUnit, $func>);

        #[pyclass(name = "PoreProfile3D", unsendable)]
        pub struct PyPoreProfile3D(PoreProfile3D<SIUnit, $func>);

        impl_3d_profile!(PyPoreProfile3D, get_x, get_y, get_z);

        #[pymethods]
        impl PyPore3D {
            #[new]
            fn new(
                functional: &$py_func,
                system_size: [PySINumber; 3],
                n_grid: [usize; 3],
                coordinates: &PySIArray2,
                sigma_ss: &PyArray1<f64>,
                epsilon_k_ss: &PyArray1<f64>,
                potential_cutoff: Option<f64>,
                cutoff_radius: Option<PySINumber>,
            ) -> Self {
                Self(Pore3D::new(
                    &functional.0,
                    [system_size[0].into(), system_size[1].into(), system_size[2].into()],
                    n_grid,
                    coordinates.clone().into(),
                    sigma_ss.to_owned_array(),
                    epsilon_k_ss.to_owned_array(),
                    potential_cutoff,
                    cutoff_radius.map(|c| c.into()),
                ))
            }

            /// Initialize the pore for the given bulk state
            ///
            /// Parameters
            /// ----------
            /// bulk : State
            ///     The bulk state in equilibrium with the pore.
            /// external_potential : Array4, optional
            ///     The external potential in the pore. Used to
            ///     save computation time in the case of costly
            ///     evaluations of external potentials.
            ///
            /// Returns
            /// -------
            /// PoreProfile3D
            #[pyo3(text_signature = "(bulk, external_potential=None)")]
            fn initialize(
                &self,
                bulk: &PyState,
                external_potential: Option<&PyArray4<f64>>,
            ) -> PyResult<PyPoreProfile3D> {
                Ok(PyPoreProfile3D(self.0.initialize(
                    &bulk.0,
                    external_potential.map(|e| e.to_owned_array()).as_ref(),
                )?))
            }
        }

        #[pymethods]
        impl PyPoreProfile3D {
            #[getter]
            fn get_grand_potential(&self) -> Option<PySINumber> {
                self.0.grand_potential.map(PySINumber::from)
            }

            #[getter]
            fn get_interfacial_tension(&self) -> Option<PySINumber> {
                self.0.interfacial_tension.map(PySINumber::from)
            }
        }
    };
}
