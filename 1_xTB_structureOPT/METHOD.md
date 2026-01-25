## xTB Method Selection

xTB implements a family of semi-empirical tight-binding methods collectively referred to as **GFN-xTB** (Geometry, Frequency, Noncovalent).

In this work, **GFN2-xTB** is used throughout unless otherwise stated.

### Available GFN Methods

| Option    | Method       | Notes                                                                                                 |
| --------- | ------------ | ----------------------------------------------------------------------------------------------------- |
| `--gfn 1` | GFN1-xTB     | Older parameterization; faster but less accurate                                                      |
| `--gfn 2` | **GFN2-xTB** | **Default and recommended**; improved treatment of geometry, noncovalent interactions, and energetics |

### Rationale for Using GFN2-xTB

GFN2-xTB provides:

* Improved accuracy for equilibrium geometries
* Better description of noncovalent interactions
* More reliable relative energies for conformational analysis
* Robust performance for organic and biomolecular systems

As a result, GFN2-xTB is widely adopted in automated workflows and is the default method used in xTB-based pipelines, including STEPs-related structure preparation.

### Default Behavior

If no method is explicitly specified, xTB defaults to **GFN2-xTB**:

```bash
xtb structure.xyz
```

Explicitly specifying the method (e.g., `--gfn 2`) is recommended in scripts and workflows to ensure reproducibility.
