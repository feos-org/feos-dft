# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- `HelmholtzEnergyFunctional`s can now overwrite the `ideal_gas` method to provide a non-default ideal gas contribution, that is accounted for in the calculation of the entropy, the internal energy and other properties. [#10](https://github.com/feos-org/feos-core/pull/10)

### Changed
- Removed the `functional` field in `Pore1D` and `Pore3D`. [#9](https://github.com/feos-org/feos-core/pull/9)

### Fixed
- Fixed the units of default values for adsorption isotherms. [#8](https://github.com/feos-org/feos-core/pull/8)

## [0.1.0] - 2021-12-22
### Added
- Initial release
