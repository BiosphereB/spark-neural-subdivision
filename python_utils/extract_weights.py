#!/usr/bin/env python3
"""
extract_weights.py
===================
Converts pre-trained PyTorch model weights from 'best_model.pt' into a binary file
for use in Ada SPARK as a lookup table for insertion angles.

Usage:
    1. Download best_model.pt from https://github.com/ugail/neural_subdivision/tree/main/outputs
    2. Place it in the outputs/ directory
    3. Run: python python_utils/extract_weights.py

Output:
    - outputs/alpha_table.bin: Binary file with float64 values for Ada
    - outputs/alpha_table.txt: Human-readable version (for debugging)
"""

import os
import torch
import numpy as np
from pathlib import Path

# --- Configuration ---
MODEL_PATH = Path("outputs/best_model.pt")  # Path to pre-trained model
OUTPUT_BIN = Path("outputs/alpha_table.bin")  # Binary output for Ada
OUTPUT_TXT = Path("outputs/alpha_table.txt")  # Human-readable output
FEATURE_DIM = 7  # Input feature dimension (delta_{j-1..j+2}, e_j/ebar, e_{j+1}/ebar, kappa)
N_EDGES = 10000  # Number of entries in lookup table (adjust as needed)

def load_model(model_path: Path) -> torch.nn.Module:
    """Load the pre-trained model from best_model.pt"""
    try:
        # Try loading as state_dict (common for PyTorch models)
        state_dict = torch.load(model_path, map_location=torch.device('cpu'))
        if isinstance(state_dict, torch.nn.Module):
            return state_dict

        # If it's a state_dict, we need to reconstruct the model architecture
        # Based on the paper: 140,505 parameters, 4 residual blocks, etc.
        class NeuralOperator(torch.nn.Module):
            def __init__(self):
                super().__init__()
                # Geometry embedding: 3 geometries -> 8-dim embedding
                self.geometry_embedding = torch.nn.Embedding(3, 8)

                # Input projection: 7 (features) + 8 (geometry) = 15 -> 128
                self.input_proj = torch.nn.Linear(15, 128)

                # 4 residual blocks (each: 128 -> 128)
                self.res_blocks = torch.nn.ModuleList([
                    torch.nn.Sequential(
                        torch.nn.Linear(128, 128),
                        torch.nn.GELU(),
                        torch.nn.LayerNorm(128),
                        torch.nn.Dropout(0.05),
                        torch.nn.Linear(128, 128)
                    )
                    for _ in range(4)
                ])

                # Output head: 128 -> 32 -> 1
                self.output_head = torch.nn.Sequential(
                    torch.nn.Linear(128, 32),
                    torch.nn.GELU(),
                    torch.nn.Linear(32, 1)
                )

            def forward(self, x, kappa):
                # x: (batch, 7) features
                # kappa: (batch,) geometry indices (-1, 0, +1 for H2, E2, S2)
                g = self.geometry_embedding(kappa + 1)  # +1 to map [-1,0,1] -> [0,1,2]
                h = torch.cat([x, g], dim=-1)
                h = self.input_proj(h)
                for block in self.res_blocks:
                    h = h + block(h)  # Residual connection
                alpha = self.output_head(h)
                # Apply sigmoid and rescale to (alpha_min, alpha_max)
                alpha = torch.sigmoid(alpha) * (np.pi/2 - 0.04) + (-np.pi/4 + 0.02)
                return alpha

        model = NeuralOperator()
        model.load_state_dict(state_dict)
        model.eval()
        return model

    except Exception as e:
        print(f"Error loading model: {e}")
        print("Make sure best_model.pt is in outputs/ and is a valid PyTorch model.")
        raise

def generate_lookup_table(model: torch.nn.Module, n_entries: int = N_EDGES) -> np.ndarray:
    """
    Generate a lookup table of insertion angles by sampling random inputs.
    This approximates the model's behavior for Ada SPARK.
    """
    # Create random input features (normalized exterior angles + edge lengths)
    rng = np.random.RandomState(42)  # Reproducible

    # Generate random features: 4 exterior angles (normalized by pi), 2 edge lengths, 1 geometry code
    features = rng.uniform(low=-1.0, high=1.0, size=(n_entries, 7))
    features[:, 4:6] = rng.uniform(low=0.5, high=2.0, size=(n_entries, 2))  # Edge lengths (positive)
    features[:, 6] = rng.choice([-1, 0, 1], size=n_entries)  # Geometry code (H2, E2, S2)

    # Convert to tensors
    x_tensor = torch.tensor(features[:, :7], dtype=torch.float32)
    kappa_tensor = torch.tensor(features[:, 6].astype(int), dtype=torch.long)

    # Predict angles
    with torch.no_grad():
        alpha_tensor = model(x_tensor, kappa_tensor)

    # Convert to numpy and clamp to G1-safe range (Theorem 1)
    alpha_table = alpha_tensor.numpy().flatten()
    alpha_min = -np.pi/4 + 0.02
    alpha_max = np.pi/4 - 0.02
    alpha_table = np.clip(alpha_table, alpha_min, alpha_max)

    return alpha_table.astype(np.float64)  # Ada expects float64

def save_binary(table: np.ndarray, output_path: Path):
    """Save the lookup table as binary file for Ada"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'wb') as f:
        table.tofile(f)
    print(f"Saved binary lookup table to {output_path} ({len(table)} entries)")

def save_text(table: np.ndarray, output_path: Path):
    """Save the lookup table as human-readable text for debugging"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        for i, alpha in enumerate(table):
            f.write(f"{i:6d}: {alpha:.10f}\n")
    print(f"Saved text lookup table to {output_path}")

def main():
    print("=== PyTorch to Ada Lookup Table Converter ===")

    # Check if model exists
    if not MODEL_PATH.exists():
        print(f"Error: Model file not found at {MODEL_PATH}")
        print("Download best_model.pt from:")
        print("https://github.com/ugail/neural_subdivision/tree/main/outputs")
        return

    # Load model
    print("Loading model...")
    model = load_model(MODEL_PATH)

    # Generate lookup table
    print(f"Generating lookup table with {N_EDGES} entries...")
    alpha_table = generate_lookup_table(model, N_EDGES)

    # Save outputs
    save_binary(alpha_table, OUTPUT_BIN)
    save_text(alpha_table, OUTPUT_TXT)

    print("\nDone! You can now use alpha_table.bin in your Ada SPARK code.")
    print("Example Ada usage:")
    print("  type Alpha_Table is array (1 .. 10000) of Real;")
    print("  Lookup_Alpha : constant Alpha_Table := ...")

if __name__ == "__main__":
    main()
