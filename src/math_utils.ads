--  math_utils.ads
--  ==============
--  Ada bindings for C functions in math_utils.c.
--  Provides numerical utilities for hyperbolic (H2), spherical (S2), and Euclidean (E2) geometries.
--
--  Usage:
--     with Math_Utils; use Math_Utils;
--     Result := C_Mobius_Add(Z, W);
--
--  Note: SPARK_Mode is Off because C bindings cannot be formally verified by GNATPROVE.

with Geometry_Subdivision; use Geometry_Subdivision;

package Math_Utils with
   SPARK_Mode => Off  -- C bindings cannot be verified by SPARK
is
   -- ===================================================================
   -- TYPE DEFINITIONS (must match C structs)
   -- ===================================================================

   -- 2D Point type (matches C's Point2D struct)
   type Point_2D is record
      X, Y : Real;
   end record;

   -- 3D Point type (for spherical geometry)
   type Point_3D is record
      X, Y, Z : Real;
   end record;

   -- ===================================================================
   -- HYPERBOLIC GEOMETRY (H2) - Poincaré Disk
   -- ===================================================================

   -- Möbius addition: z ⊕ w (Eq. 3 from the paper)
   function C_Mobius_Add (Z, W : Point_2D) return Point_2D
     with Import,
          Convention => C,
          External_Name => "mobius_add";

   -- Hyperbolic distance: d_D(z, w) = 2 * arctanh(||(-z) ⊕ w||) (Eq. 3)
   function C_Hyperbolic_Distance (Z, W : Point_2D) return Real
     with Import,
          Convention => C,
          External_Name => "hyperbolic_distance";

   -- Unit hyperbolic tangent vector from z to w (for Eq. 8)
   function C_Hyperbolic_Tangent (Z, W : Point_2D) return Point_2D
     with Import,
          Convention => C,
          External_Name => "hyperbolic_tangent";

   -- Clamp a point to the Poincaré disk (||z|| < 0.999)
   function C_Clamp_To_Disk (Z : Point_2D) return Point_2D
     with Import,
          Convention => C,
          External_Name => "clamp_to_disk";

   -- ===================================================================
   -- EUCLIDEAN GEOMETRY (E2)
   -- ===================================================================

   -- Euclidean distance between two points
   function C_Euclidean_Distance (A, B : Point_2D) return Real
     with Import,
          Convention => C,
          External_Name => "euclidean_distance";

   -- Normalize a 2D vector
   function C_Normalize_Vector (V : Point_2D) return Point_2D
     with Import,
          Convention => C,
          External_Name => "normalize_vector";

   -- Rotate a 2D vector by angle alpha (in radians)
   function C_Rotate_Vector (V : Point_2D; Alpha : Real) return Point_2D
     with Import,
          Convention => C,
          External_Name => "rotate_vector";

   -- ===================================================================
   -- SPHERICAL GEOMETRY (S2)
   -- ===================================================================

   -- Spherical distance (great-circle distance) between two unit vectors
   function C_Spherical_Distance (
      AX, AY, AZ : Real;
      BX, BY, BZ : Real) return Real
     with Import,
          Convention => C,
          External_Name => "spherical_distance";

   -- Spherical tangent vector from point P to Q (Eq. 2 from the paper)
   -- Returns the unit tangent vector in the direction from P to Q on S2
   procedure C_Spherical_Tangent (
      PX, PY, PZ : Real;
      QX, QY, QZ : Real;
      TX, TY, TZ : out Real)
     with Import,
          Convention => C,
          External_Name => "spherical_tangent";

   -- ===================================================================
   -- UTILITY TYPES AND CONSTANTS
   -- ===================================================================

   -- Real type matching C's double precision
   type Real is digits 15 range -1.0E308 .. 1.0E308;

   -- Pi constant (for convenience)
   Pi : constant Real := Ada.Numerics.Pi;

private
   -- Private declarations (none for this package)

end Math_Utils;
