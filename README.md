# Ada SPARK Implementation: Neural Adaptive Tension for Multi-Geometry Curve Subdivision

**Ada SPARK port of the geometric core and subdivision operators** from:
> **"Neural Adaptive Tension for Multi-Geometry Curve Subdivision: A Unified Approach"**
> H. Ugail and N. Howard
> [Original Repository](https://github.com/ugail/neural_subdivision)

---

## 🎯 Overview

This repository provides an **Ada SPARK implementation** of the **geometric foundations** and **subdivision operators** for interpolatory curve subdivision across three constant-curvature geometries:
- **Euclidean plane (E²)**
- **Unit two-sphere (S²)**
- **Hyperbolic plane (H², Poincaré disk model)**

### ⚠️ Scope
This is a **partial reimplementation** focusing on:
✅ **Deterministic geometric kernels** (points, vectors, metrics, exterior angles)
✅ **Classical subdivision schemes** (4-point, 6-point)
✅ **Angle-based insertion operators** (with G¹ safety bounds)
✅ **Formal verification** (GNATPROVE-compatible SPARK code)

❌ **Not included** (replaced by lookup tables or heuristics):
- Neural network training (Python/TensorFlow)
- Dynamic weight prediction (replaced by **pre-trained lookup tables** from the original model)

The original neural predictor (140K parameters) is **not implemented in Ada SPARK** (not feasible for formal verification).
Instead, we use **pre-trained insertion angles** (extracted from `best_model.pt`) as static lookup tables.

---

## 📦 Repository Structure

```
spark-neural-subdivision/
├── src/
│   ├── geometry_core.ads      # SPARK: Types, specifications, and invariants
│   ├── geometry_core.adb      # SPARK: Implementations (E², S², H²)
│   ├── math_utils.ads         # C bindings for numerical stability (e.g., Möbius addition)
│   ├── alpha_table.ads        # Lookup tables for insertion angles (from best_model.pt)
│   └── main.adb               # Example usage
├── c_src/
│   └── math_utils.c           # C implementations (e.g., hyperbolic operations)
├── python_utils/
│   └── extract_weights.py      # Convert PyTorch weights to Ada-compatible binary
├── gnatprove/
│   └── geometry_subdivision.gpr # GNATPROVE project file
├── outputs/                   # (Optional) Pre-trained lookup tables
│   └── alpha_table.bin        # Binary file for insertion angles
└── README.md
```

---

## 🛠️ Requirements

### For Ada SPARK Development
   Tool | Version | Purpose |
 |------|---------|---------|
 | [GNAT Community Edition](https://www.adacore.com/download) | 2023+ | Ada compiler + SPARK |
 | [GNATPROVE](https://docs.adacore.com/live/wave/gprbuild/html/gprbuild_ug/en/gprbuild.html) | Latest | Formal verification |
 | [Alire](https://alire.ada.dev/) | Optional | Dependency management |

Install GNAT (Linux):
```bash
sudo apt-get install gnat gnatprove
```

### For Lookup Table Generation (Optional)
 | Package | Minimum Version |
 |---------|-----------------|
 | Python  | 3.9             |
 | PyTorch | 2.0             |
 | NumPy   | 1.24            |

Install dependencies:
```bash
pip install torch numpy
```

---

## 🚀 Quick Start

### 1. Clone the Repository
```bash
git clone https://github.com/BiosphereB/spark-neural-subdivision.git
cd spark-neural-subdivision
```

### 2. Generate Lookup Tables (Optional)
If you want to use the **pre-trained insertion angles** from the original model:
1. Download `best_model.pt` from the [original repository](https://github.com/ugail/neural_subdivision/tree/main/outputs).
2. Run the conversion script:
   ```bash
   python python_utils/extract_weights.py
   ```
   This creates `outputs/alpha_table.bin` for use in Ada.

**Alternative:** Use a **heuristic** (e.g., linear adaptive tension) instead of lookup tables.

### 3. Compile & Verify with GNATPROVE
```bash
cd gnatprove
gnatprove -P geometry_subdivision.gpr
```
- **Expected:** Some warnings for non-linear math (e.g., `sin`, `cos`).
  Add [SPARK lemmas](#spark-lemmas) as needed.

### 4. Run Tests
```bash
gnatmake -P geometry_subdivision.gpr
./obj/main
```

---
## 📊 Ada SPARK Implementation Details

### Geometry Types
 | Type | Description | SPARK Annotations |
 |------|-------------|-------------------|
 | `Point_2D` | 2D point (for E² and H²) | `X, Y : Real` |
 | `Point_3D` | 3D point (for S²) | `X, Y, Z : Real` with `X² + Y² + Z² = 1.0` |
 | `Universal_Point` | Discriminated union for all geometries | `G : Geometry_Type` |
 | `Polygon` | Cyclic sequence of points | `N : Polygon_Index` |

### Key Functions
 | Function | Description | Geometry | SPARK Verified |
 |----------|-------------|----------|----------------|
 | `Euclidean_Distance` | Distance between 2D points | E² | ✅ |
 | `Spherical_Distance` | Great-circle distance | S² | ✅ |
 | `Hyperbolic_Distance` | Poincaré disk distance | H² | ✅ (with C binding) |
 | `Exterior_Angle_E2` | Exterior angle at vertex | E² | ✅ |
 | `Exterior_Angle_S2` | Exterior angle at vertex | S² | ⚠️ (Needs lemmas) |
 | `Exterior_Angle_H2` | Exterior angle at vertex | H² | ⚠️ (Needs lemmas) |
 | `Classical_4Point_Subdivide` | Classical 4-point scheme | All | ✅ |
 | `Angle_Based_Insert` | Insert vertex using angle α | All | ✅ |
 | `Clamp_Alpha` | Enforce G¹ safety bound (Theorem 1) | All | ✅ |

### SPARK Lemmas
For non-linear math (e.g., trigonometric identities), add lemmas to `geometry_core.ads`:
```ada
-- Example: Pythagorean identity
lemma Sin_Cos_Identity (X : Real) with
  Post => Sin(X)**2 + Cos(X)**2 = 1.0;

-- Example: Non-negativity of square root
lemma Sqrt_Nonnegative (X : Real) with
  Pre  => X >= 0.0,
  Post => Sqrt(X) >= 0.0;
```

---
## 🔄 Relationship to Original Work

### What’s Included
- **Geometric kernels** (E², S², H² metrics and operations).
- **Subdivision operators** (4-point, 6-point, angle-based insertion).
- **Safety guarantees** (G¹ bounds via `Clamp_Alpha`).

### What’s Replaced
 | Original (Python) | This Project (Ada SPARK) |
 |-------------------|--------------------------|
 | Neural network (`NeuralOperator`) | **Lookup tables** (`alpha_table.bin`) |
 | Dynamic weight prediction | **Static angle values** (from `best_model.pt`) |
 | Training loop | **Not applicable** (use pre-trained weights) |

### Performance Notes
- **Ada SPARK is slower** than Python for numerical computations.
  - **Optimization:** Use `pragma Inline` for critical functions (e.g., `Euclidean_Distance`).
  - **C Bindings:** Offload numerically unstable ops (e.g., Möbius addition) to C.
- **Formal verification** adds overhead but ensures correctness.

---
## 📚 Original Paper & Citation

If you use this code in research, please cite the original paper:
```bibtex
@article{ugail2026neural,
  title={A Neural Tension Operator for Curve Subdivision across Constant Curvature Geometries},
  author={Ugail, Hassan and Howard, Newton},
  journal={arXiv preprint},
  year={2026},
  note={arXiv:2604.XXXX}
}
```

**Original Repository:** [ugail/neural_subdivision](https://github.com/ugail/neural_subdivision)

---
## 🤝 Acknowledgments

- **Hassan Ugail & Newton Howard**: For the original research and open-sourcing the Python implementation.
- **AdaCore**: For GNAT and SPARK, enabling formal verification of geometric algorithms.

---
## 📜 License

This project is **for personal/hobby use only**.
The original code is released for **research and academic use** (see [ugail/neural_subdivision](https://github.com/ugail/neural_subdivision)).
**Do not use this Ada SPARK port for commercial purposes without permission.**
