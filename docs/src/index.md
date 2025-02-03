```@raw html
<div align="center">

<a href="https://github.com/kestrelquantum/Piccolo.jl">
  <img src="assets/logo.svg" alt="Piccolo.jl" width="25%"/>
</a> 

<div style="display: table; width: 100%;">
  <div style="display: table-row;">
    <div style="display: table-cell; text-align: center;"><b>Documentation</b></div>
    <div style="display: table-cell; text-align: center;"><b>Build Status</b></div>
    <div style="display: table-cell; text-align: center;"><b>Support</b></div>
  </div>
  <div style="display: table-row;">
    <div style="display: table-cell; text-align: center;">
      <a href="https://kestrelquantum.github.io/PiccoloQuantumObjects.jl/stable/">
        <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"/>
      </a>
      <a href="https://kestrelquantum.github.io/PiccoloQuantumObjects.jl/dev/">
        <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"/>
      </a>
    </div>
    <div style="display: table-cell; text-align: center;">
      <a href="https://github.com/kestrelquantum/PiccoloQuantumObjects.jl/actions/workflows/CI.yml?query=branch%3Amain">
        <img src="https://github.com/kestrelquantum/PiccoloQuantumObjects.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"/>
      </a>
      <a href="https://codecov.io/gh/kestrelquantum/PiccoloQuantumObjects.jl">
        <img src="https://codecov.io/gh/kestrelquantum/PiccoloQuantumObjects.jl/branch/main/graph/badge.svg" alt="Coverage"/>
      </a>
    </div>
    <div style="display: table-cell; text-align: center;">
      <a href="https://unitary.fund">
        <img src="https://img.shields.io/badge/Supported%20By-Unitary%20Fund-FFFF00.svg" alt="Unitary Fund"/>
      </a>
    </div>
  </div>
</div>

<br>
<i> Make simple transformation of complex objects for quantum numerics </i>
<br>

</div>
```

# PiccoloQuantumObjects

**PiccoloQuantumObjects.jl** is a Julia package for working with quantum objects. It provides tools for constructing and manipulating quantum states and operators. It is designed to be used with other packages in the [Piccolo.jl](https://github.com/kestrelquantum/Piccolo.jl) ecosystem, such as [QuantumCollocation.jl](https://github.com/kestrelquantum/QuantumCollocation.jl) and [NamedTrajectories.jl](https://github.com/kestrelquantum/NamedTrajectories.jl).

## Installation

This package is registered! To install, enter the Julia REPL, type `]` to enter pkg mode, and then run:
```julia
pkg> add PiccoloQuantumObjects
```

## Usage

The following example demonstrates how to create a quantum state, create a quantum operator, and apply the operator to the state:

```Julia
using PiccoloQuantumObjects

# Create a quantum state
state = ket_from_string("g", [2])

# Create a quantum operator
operator = PAULIS.X

# Apply the operator to the state
new_state = operator * state
```
