# PiccoloQuantumObjects

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kestrelquantum.github.io/PiccoloQuantumObjects.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kestrelquantum.github.io/PiccoloQuantumObjects.jl/dev/)
[![Build Status](https://github.com/kestrelquantum/PiccoloQuantumObjects.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kestrelquantum/PiccoloQuantumObjects.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kestrelquantum/PiccoloQuantumObjects.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kestrelquantum/PiccoloQuantumObjects.jl)

## Description

PiccoloQuantumObjects.jl is a Julia package for working with quantum objects, providing tools for constructing and manipulating quantum states and operators. It is designed to be used with other packages in the Piccolo ecosystem, such as [QuantumCollocation.jl](https://github.com/kestrelquantum/QuantumCollocation.jl) and [NamedTrajectories.jl](https://github.com/kestrelquantum/NamedTrajectories.jl).

## Installation

This package is registered! To install, enter the Julia REPL, type `]` to enter pkg mode, and then run:
```julia
pkg> add PiccoloQuantumObjects
```

## Usage

The following example demonstrates how to create a quantum state, create a quantum operator, and apply the operator to the state:

```julia
using PiccoloQuantumObjects

# Create a quantum state
state = ket_from_string("g", [2])

# Create a quantum operator
operator = PAULIS[:X]

# Apply the operator to the state
new_state = operator * state
```

## License

PiccoloQuantumObjects.jl is licensed under the MIT License.
