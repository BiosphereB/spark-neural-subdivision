--  geometry_core.ads
--  ==================
--  Specification of geometric types and subdivision operators
--  for Euclidean (E2), Spherical (S2), and Hyperbolic (H2) geometries.
--
--  This package provides:
--  1. Type definitions for points, polygons, and geometry types
--  2. Basic geometry operations (distance, angles)
--  3. Subdivision operators (4-point, angle-based insertion)
--  4. Safety guarantees (G1 bounds via Clamp_Alpha)
--
--  Usage:
--     with Geometry_Subdivision; use Geometry_Subdivision;
--
--  SPARK: All functions include Pre/Post conditions for formal verification

with Ada.Numerics; use Ada.Numerics;

package Geometry_Subdivision with
   SPARK_Mode => On
is
   -- ===================================================================
   -- TYPE DEFINITIONS
   -- ===================================================================

   -- Precision for all calculations (matches C's double)
   type Real is digits 15 range -1.0E308 .. 1.0E308;

   -- Geometry types
   type Geometry_Type is (E2, S2, H2);  -- Euclidean, Spherical, Hyperbolic

   -- 2D Point (for E2 and H2)
   type Point_2D is record
      X, Y : Real;
   end record;

   -- 3D Point (for S2)
   type Point_3D is record
      X, Y, Z : Real;
   end record;

   -- Universal point type (discriminated by geometry)
   type Universal_Point (G : Geometry_Type := E2) is record
      case G is
         when E2 | H2 =>
            P2 : Point_2D;
         when S2 =>
            P3 : Point_3D;
      end case;
   end record;

   -- Polygon index type (max 1000 vertices)
   type Polygon_Index is range 1 .. 1000;

   -- Polygon type (cyclic sequence of points)
   type Polygon (G : Geometry_Type; N : Polygon_Index) is
     array (1 .. N) of Universal_Point(G);

   -- ===================================================================
   -- BASIC GEOMETRY OPERATIONS (Specifications)
   -- ===================================================================

   -- Euclidean distance between two 2D points
   function Euclidean_Distance (A, B : Point_2D) return Real
     with Post => Euclidean_Distance(A, B) >= 0.0;

   -- Spherical distance (great-circle distance) between two unit vectors
   function Spherical_Distance (A, B : Point_3D) return Real
     with Pre  => A.X**2 + A.Y**2 + A.Z**2 = 1.0 and
                 B.X**2 + B.Y**2 + B.Z**2 = 1.0,
          Post => Spherical_Distance(A, B) >= 0.0;

   -- Hyperbolic distance in Poincaré disk
   function Hyperbolic_Distance (A, B : Point_2D) return Real
     with Pre  => A.X**2 + A.Y**2 < 1.0 and B.X**2 + B.Y**2 < 1.0,
          Post => Hyperbolic_Distance(A, B) >= 0.0;

   -- Möbius addition (for H2)
   function Mobius_Add (A, B : Point_2D) return Point_2D
     with Pre  => A.X**2 + A.Y**2 < 1.0 and B.X**2 + B.Y**2 < 1.0;

   -- ===================================================================
   -- EXTERIOR ANGLE CALCULATIONS (Paper Section 3.2)
   -- ===================================================================

   -- Exterior angle at vertex J for E2 (Paper Eq. 1)
   function Exterior_Angle_E2 (P : Polygon(E2, N); J : Polygon_Index) return Real
     with Pre => N >= 3 and J in 1 .. N;

   -- Exterior angle at vertex J for S2 (Paper Eq. 2)
   function Exterior_Angle_S2 (P : Polygon(S2, N); J : Polygon_Index) return Real
     with Pre => N >= 3 and J in 1 .. N;

   -- Exterior angle at vertex J for H2 (Paper Eq. 3 + Euclidean formula)
   function Exterior_Angle_H2 (P : Polygon(H2, N); J : Polygon_Index) return Real
     with Pre => N >= 3 and J in 1 .. N;

   -- ===================================================================
   -- SUBDIVISION OPERATORS (Paper Section 3.3)
   -- ===================================================================

   -- Classical 4-point subdivision scheme (μ = 0 by default)
   procedure Classical_4Point_Subdivide (
      P      : in  Polygon(G, N);
      P_New  : out Polygon(G, 2 * N);
      Mu     : Real := 0.0)  -- Tension parameter
     with Pre  => N >= 3,
          Post => P_New'Length = 2 * N;

   -- Angle-based vertex insertion (Paper Eq. 6-8)
   procedure Angle_Based_Insert (
      P      : in  Polygon(G, N);
      J      : in  Polygon_Index;
      Alpha  : in  Real;
      Q      : out Universal_Point(G))
     with Pre => N >= 3 and J in 1 .. N and
                 Alpha > -Pi/4 + 0.02 and Alpha < Pi/4 - 0.02;

   -- ===================================================================
   -- SAFETY GUARANTEES (Theorem 1)
   -- ===================================================================

   -- Clamp alpha to G1-safe range (Paper Eq. 12)
   function Clamp_Alpha (Alpha : Real) return Real
     with Post => Clamp_Alpha(Alpha) > -Pi/4 + 0.02 and
                  Clamp_Alpha(Alpha) < Pi/4 - 0.02;

   -- ===================================================================
   -- UTILITY TYPES AND CONSTANTS
   -- ===================================================================

   -- Feature vector type (Paper Eq. 17)
   -- [delta_{j-1}/pi, delta_j/pi, delta_{j+1}/pi, delta_{j+2}/pi,
   --  e_j/ebar, e_{j+1}/ebar, kappa]
   type Feature_Vector is array (1 .. 7) of Real;

   -- Pi constant
   Pi : constant Real := Ada.Numerics.Pi;

   -- ===================================================================
   -- LEMMAS FOR SPARK VERIFICATION
   -- ===================================================================

   -- Pythagorean identity
   lemma Sin_Cos_Identity (X : Real) with
     Post => Sin(X)**2 + Cos(X)**2 = 1.0;

   -- Non-negativity of square root
   lemma Sqrt_Nonnegative (X : Real) with
     Pre  => X >= 0.0,
     Post => Sqrt(X) >= 0.0;

   -- Range of Arctan2
   lemma Arctan2_Range (Y, X : Real) with
     Post => Arctan2(Y, X) >= -Pi and Arctan2(Y, X) <= Pi;

   -- Range of Arcsin
   lemma Arcsin_Range (X : Real) with
     Pre  => X >= -1.0 and X <= 1.0,
     Post => Arcsin(X) >= -Pi/2 and Arcsin(X) <= Pi/2;

   -- Range of Arccos
   lemma Arccos_Range (X : Real) with
     Pre  => X >= -1.0 and X <= 1.0,
     Post => Arccos(X) >= 0.0 and Arccos(X) <= Pi;

private
   -- Private declarations (none for this package)

end Geometry_Subdivision;
