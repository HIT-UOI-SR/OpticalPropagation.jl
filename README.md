# OpticalPropagation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HIT-UOI-SR.github.io/OpticalPropagation.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HIT-UOI-SR.github.io/OpticalPropagation.jl/dev)
[![Build Status](https://github.com/HIT-UOI-SR/OpticalPropagation.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/HIT-UOI-SR/OpticalPropagation.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/HIT-UOI-SR/OpticalPropagation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HIT-UOI-SR/OpticalPropagation.jl)

A package of useful optical propagation methods and simple 2-dimensional light field data types.

Requires Julia 1.10 or later. CI tests the latest Julia LTS and stable releases on
Linux, macOS, and Windows, with an additional experimental nightly run on Linux.

## Development

Instantiate the package environment and run the tests:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Build the documentation with the local package:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Set `GKSwstype=100` when building documentation without a display. Documentation
deployment runs separately in CI on pushes to `master` and version tags.

Package, test, and documentation environments declare their own compatibility
bounds. Their generated manifests stay local so each supported Julia version can
resolve compatible dependencies. Dependabot checks Julia dependencies and GitHub
Actions weekly.
