# Neural Adaptive Tension for Multi-Geometry Curve Subdivision

Code repository for:

> **"Neural Adaptive Tension for Multi-Geometry Curve Subdivision: A Unified Approach"**  
> H. Ugail and N. Howard  
> 

---

## Overview

This repository provides a self-contained implementation of a shared learned tension predictor for interpolatory curve subdivision across three constant-curvature geometries: the Euclidean plane **E²**, the unit two-sphere **S²**, and the hyperbolic plane **H²** (Poincaré disk model). A single 140,505-parameter residual network predicts per-edge insertion angles from local intrinsic features and a trainable geometry embedding, replacing the fixed global tension parameter of classical four- and six-point schemes.

The method achieves markedly lower bending energy and G¹ roughness than all fixed-tension baselines, and generalises to out-of-distribution data as demonstrated on the ISS orbital ground-track example.

---

## Requirements

| Package | Minimum version |
|---------|----------------|
| Python  | 3.9            |
| PyTorch | 2.0            |
| NumPy   | 1.24           |
| Matplotlib | 3.7         |
| SciPy   | 1.10           |

Install dependencies:

```bash
pip install torch numpy matplotlib scipy
```

The notebook is designed to run on **Google Colab** with an A100 GPU. CPU execution is supported but training will take considerably longer (estimated 28 hours on A100; scale accordingly on CPU).

---

## Repository Structure

```
Neural_Subdivision.ipynb   Main notebook (all code, training, evaluation)
outputs/
  best_model.pt            Pre-trained checkpoint (epoch 200)
  fig1_qualitative_E2.png  Qualitative results — Euclidean
  fig2_qualitative_S2.png  Qualitative results — Spherical
  fig3_qualitative_H2.png  Qualitative results — Hyperbolic
  fig4_pareto.png          Fidelity–smoothness Pareto frontier
  fig5_realworld_iss.png   ISS orbital ground-track case study
  fig6_ablation.png        Geometry embedding and loss ablations
```

The notebook is structured as a single sequential document. All cells should be run top to bottom in one session. The pre-trained checkpoint at `outputs/best_model.pt` allows evaluation and inference to be run independently of training (see **Skipping Training** below).

---

## Notebook Structure

| Section | Content |
|---------|---------|
| Install / check dependencies | Version verification, Google Drive mount |
| Imports & device setup | Packages, random seeds, output directory |
| Hyperparameters | All configurable settings in `CFG` dict |
| Differentiable geometry (E², S², H²) | Edge lengths, exterior angles, subdivision steps for all three spaces |
| Neural operator architecture | `NeuralOperator`: geometry embedding, input projection, 4 residual blocks, G¹-safe sigmoid head |
| Synthetic dataset | 400 curves per geometry; Euclidean ellipses / Fourier / superellipses; spherical Lissajous / polar-Fourier / perturbed great circles; hyperbolic circles / offset ellipses |
| Subdivision engine & loss functions | Classical warm-up, neural steps, Chamfer / smoothness / bending / equivariance losses |
| Evaluation metrics | Mean nearest-neighbour, Hausdorff, bending energy, G¹ proxy |
| Early stopping & LR scheduler | Two-level strategy: LR halving + full stop |
| Training loop | Geometry-grouped batch processing, gradient clipping, checkpoint saving |
| Build, baseline, train | Model construction, baseline evaluation, full training run |
| Evaluate & print results table | Aggregate metrics across all geometries and methods |
| Qualitative figures | Six multi-panel publication figures |
| ISS real-world case study (S²) | Out-of-distribution test on closed sinusoidal orbital ground track |
| OOD evaluation | Systematic out-of-distribution test across all three geometries |
| Custom polygon inference | Subdivide your own control polygon with the trained model |

---

## Hyperparameters

All settings are consolidated in the `CFG` dictionary. Key parameters:

```python
CFG = dict(
    # Data
    n_per_geometry = 400,     # curves per geometry (E², S², H²)
    n_ctrl         = 12,      # control polygon vertices
    n_gt           = 1000,    # ground-truth curve resolution

    # Model
    hidden         = 128,
    n_res_blocks   = 4,
    dropout        = 0.05,

    # Training
    n_epochs       = 300,
    lr             = 1e-3,
    batch_size     = 8,
    grad_clip      = 0.5,
    save_every     = 50,

    # Loss weights
    lambda_chamfer = 1.00,
    lambda_bending = 1e-4,
    lambda_equiv   = 0.10,
    lambda_smooth_E = 0.05,   # Euclidean smoothness weight
    lambda_smooth_S = 0.15,   # Spherical  smoothness weight (elevated)
    lambda_smooth_H = 0.05,   # Hyperbolic smoothness weight
    ...
)
```

The elevated spherical smoothness weight (`lambda_smooth_S = 0.15`) compensates for the additional angular variation from geodesic curvature on S², which would otherwise under-penalise smoothness relative to the other geometries.

---

## Skipping Training

To use the pre-trained model directly, place `best_model.pt` in the `outputs/` directory, then run only the following sections in order:

1. Install / check dependencies  
2. Imports & device setup  
3. Hyperparameters  
4. Differentiable geometry  
5. Neural operator architecture  
6. Subdivision engine & loss functions  
7. Evaluation metrics  
8. Evaluate & print results table (and/or any figures section)

The checkpoint was saved at epoch 200 with a cross-geometry mean nearest-neighbour distance of 0.02325.

---

## Custom Polygon Inference

The final cell allows subdivision of an arbitrary control polygon. Set `geometry` to `'euclidean'`, `'spherical'`, or `'hyperbolic'` and supply a control polygon in the appropriate format:

```python
# Euclidean: (N, 2) float array
# Spherical: (N, 3) float array of unit vectors
# Hyperbolic: (N, 2) float array with ||z|| < 1  (Poincaré disk)

my_polygon = torch.tensor([
    [1., 0.], [.7, .7], [0., 1.], [-.7, .7],
    [-1., 0.], [-.7, -.7], [0., -1.], [.7, -.7],
], dtype=torch.float32, device=DEVICE)

geometry = 'euclidean'
n_levels = 5

with torch.no_grad():
    refined = neural_subdivide(my_polygon, model, geometry, n_iters=n_levels)
```

---

## Baselines

Three comparison methods are implemented alongside the neural predictor:

| Method | Description |
|--------|-------------|
| **4-point DGL** (μ = 0) | Classical Dyn–Gregory–Levin scheme, C² limit curves |
| **6-point** (μ = −0.25) | Six-point interpolatory scheme, C² limit curves |
| **Log-exp lifts** | Riemannian manifold lifts of both classical schemes via Wallner–Dyn proximity construction |
| **Best-μ oracle** | Grid search over 21 tension values on the validation set |
| **Linear adaptive heuristic** | Curvature-proportional per-edge tension, no learned parameters |

---

## Results Summary

On 240 held-out validation curves (80 per geometry):

| Method | Mean-NN ↓ | Bending energy ↓ | G¹ proxy ↓ |
|--------|-----------|-----------------|-----------|
| 4-point DGL | baseline | baseline | baseline |
| 6-point | similar | higher | higher |
| Log-exp lift | lower | higher | higher |
| **Neural predictor** | modest increase | **−88% vs DGL** | **−91% vs DGL** |

The neural predictor occupies a distinct position on the fidelity–smoothness Pareto frontier: it accepts a modest mean-NN increase relative to manifold lifts in exchange for substantially smoother curvature profiles.

On the out-of-distribution ISS ground-track example (S²): bending energy −41%, G¹ proxy −68%, Hausdorff distance +69% vs four-point baseline. The absolute Hausdorff increase is 0.00285, less than 0.3% of the sphere circumference.

---

## Citation


```

---

## Licence

This code is released for research and academic use. Please cite the above manuscript when using or building on this work.
