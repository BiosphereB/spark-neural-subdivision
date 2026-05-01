--  geometry_core.adb
--  ===================
--  Implementation of geometric operations and subdivision operators
--  for Euclidean (E2), Spherical (S2), and Hyperbolic (H2) geometries.
--
--  This file contains the core mathematical operations for curve subdivision,
--  including distance calculations, exterior angles, and vertex insertion.

with Math_Utils; use Math_Utils;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;

package body Geometry_Subdivision with
   SPARK_Mode => On
is

   -- ===================================================================
   -- BASIC GEOMETRY OPERATIONS
   -- ===================================================================

   function Euclidean_Distance (A, B : Point_2D) return Real is
   begin
      return C_Euclidean_Distance(A, B);
   end Euclidean_Distance;

   function Spherical_Distance (A, B : Point_3D) return Real is
   begin
      return C_Spherical_Distance(A.X, A.Y, A.Z, B.X, B.Y, B.Z);
   end Spherical_Distance;

   function Hyperbolic_Distance (A, B : Point_2D) return Real is
   begin
      return C_Hyperbolic_Distance(A, B);
   end Hyperbolic_Distance;

   function Mobius_Add (A, B : Point_2D) return Point_2D is
   begin
      return C_Mobius_Add(A, B);
   end Mobius_Add;

   -- ===================================================================
   -- EXTERIOR ANGLE CALCULATIONS (Paper Section 3.2)
   -- ===================================================================

   function Exterior_Angle_E2 (P : Polygon(E2, N); J : Polygon_Index) return Real is
      -- Tangent vectors (Paper Eq. 1)
      T_In, T_Out : Point_2D;
      Prev_Index, Next_Index : Polygon_Index;
   begin
      Prev_Index := (if J = 1 then N else J - 1);
      Next_Index := (if J = N then 1 else J + 1);

      -- Outgoing tangent: P[J] -> P[J+1]
      T_Out.X := P(Next_Index).P2.X - P(J).P2.X;
      T_Out.Y := P(Next_Index).P2.Y - P(J).P2.Y;

      -- Incoming tangent: P[J-1] -> P[J]
      T_In.X := P(J).P2.X - P(Prev_Index).P2.X;
      T_In.Y := P(J).P2.Y - P(Prev_Index).P2.Y;

      -- Normalize tangents
      declare
         Norm_Out : constant Real := Euclidean_Distance(T_Out, (0.0, 0.0));
         Norm_In  : constant Real := Euclidean_Distance(T_In, (0.0, 0.0));
         T_Out_Norm : Point_2D;
         T_In_Norm  : Point_2D;
      begin
         if Norm_Out > 1e-12 then
            T_Out_Norm := (T_Out.X / Norm_Out, T_Out.Y / Norm_Out);
         else
            T_Out_Norm := (0.0, 0.0);
         end if;

         if Norm_In > 1e-12 then
            T_In_Norm := (T_In.X / Norm_In, T_In.Y / Norm_In);
         else
            T_In_Norm := (0.0, 0.0);
         end if;

         -- atan2(det, dot) for angle between vectors
         return Arctan2(
           T_In_Norm.X * T_Out_Norm.Y - T_In_Norm.Y * T_Out_Norm.X,
           T_In_Norm.X * T_Out_Norm.X + T_In_Norm.Y * T_Out_Norm.Y
         );
      end;
   end Exterior_Angle_E2;

   function Exterior_Angle_S2 (P : Polygon(S2, N); J : Polygon_Index) return Real is
      -- Spherical tangent vectors (Paper Eq. 2)
      V, T : Point_3D;
      Prev_Index, Next_Index : Polygon_Index;
      TX, TY, TZ : Real;
   begin
      Prev_Index := (if J = 1 then N else J - 1);
      Next_Index := (if J = N then 1 else J + 1);

      -- Outgoing tangent: P[J] -> P[J+1]
      C_Spherical_Tangent(
         P(J).P3.X, P(J).P3.Y, P(J).P3.Z,
         P(Next_Index).P3.X, P(Next_Index).P3.Y, P(Next_Index).P3.Z,
         TX, TY, TZ
      );
      T := (TX, TY, TZ);

      -- Incoming tangent: P[J-1] -> P[J]
      C_Spherical_Tangent(
         P(Prev_Index).P3.X, P(Prev_Index).P3.Y, P(Prev_Index).P3.Z,
         P(J).P3.X, P(J).P3.Y, P(J).P3.Z,
         TX, TY, TZ
      );
      V := (TX, TY, TZ);

      -- Exterior angle: sin(delta) = (V × T) · P[J]
      declare
         Cross_X : constant Real := V.Y * T.Z - V.Z * T.Y;
         Cross_Y : constant Real := V.Z * T.X - V.X * T.Z;
         Cross_Z : constant Real := V.X * T.Y - V.Y * T.X;
         Dot_Product : constant Real := Cross_X * P(J).P3.X +
                                       Cross_Y * P(J).P3.Y +
                                       Cross_Z * P(J).P3.Z;
      begin
         return Arcsin(Dot_Product);
      end;
   end Exterior_Angle_S2;

   function Exterior_Angle_H2 (P : Polygon(H2, N); J : Polygon_Index) return Real is
      -- Hyperbolic tangent vectors (using Euclidean formula from Paper Eq. 1)
      V, T : Point_2D;
      Prev_Index, Next_Index : Polygon_Index;
   begin
      Prev_Index := (if J = 1 then N else J - 1);
      Next_Index := (if J = N then 1 else J + 1);

      -- Outgoing tangent: P[J] -> P[J+1]
      T := C_Hyperbolic_Tangent(P(J).P2, P(Next_Index).P2);

      -- Incoming tangent: P[J-1] -> P[J]
      V := C_Hyperbolic_Tangent(P(Prev_Index).P2, P(J).P2);
      V := C_Normalize_Vector((-V.X, -V.Y));  -- Reverse direction

      -- Use Euclidean formula for angle between vectors (Paper Eq. 1)
      declare
         Dot_Product : constant Real := V.X * T.X + V.Y * T.Y;
         Det : constant Real := V.X * T.Y - V.Y * T.X;
      begin
         return Arctan2(Det, Dot_Product);
      end;
   end Exterior_Angle_H2;

   -- ===================================================================
   -- SUBDIVISION OPERATORS (Paper Section 3.3)
   -- ===================================================================

   procedure Classical_4Point_Subdivide (
      P      : in  Polygon(G, N);
      P_New  : out Polygon(G, 2 * N);
      Mu     : Real := 0.0)  -- Tension parameter
   is
      -- Insert new vertices at each edge using classical rule (Paper Eq. 9)
   begin
      case G is
         when E2 =>
            for J in 1 .. N loop
               -- Retain original vertex
               P_New(2 * J - 1) := P(J);

               -- Insert new vertex at edge J (classical 4-point rule)
               declare
                  Prev_Prev_Index : constant Polygon_Index :=
                    (if J-2 < 1 then N + (J-2) else J-2);
                  Prev_Index       : constant Polygon_Index :=
                    (if J-1 < 1 then N else J-1);
                  Next_Index       : constant Polygon_Index :=
                    (if J+1 > N then 1 else J+1);
                  Next_Next_Index  : constant Polygon_Index :=
                    (if J+2 > N then J+2 - N else J+2);

                  -- Exterior angles (simplified for E2)
                  Delta_Jm1 : constant Real := Exterior_Angle_E2(P, Prev_Prev_Index);
                  Delta_J   : constant Real := Exterior_Angle_E2(P, Prev_Index);
                  Delta_Jp1 : constant Real := Exterior_Angle_E2(P, J);
                  Delta_Jp2 : constant Real := Exterior_Angle_E2(P, Next_Index);

                  -- Classical alpha (Paper Eq. 9)
                  Alpha : constant Real :=
                    0.125 * (Mu * (Delta_Jm1 + Delta_Jp2) + (1.0 - Mu) * (Delta_J + Delta_Jp1));
               begin
                  -- Insert new vertex using angle-based rule (E2 version, Paper Eq. 6)
                  Angle_Based_Insert(P, J, Alpha, P_New(2 * J));
               end;
            end loop;

         when S2 =>
            for J in 1 .. N loop
               -- Retain original vertex
               P_New(2 * J - 1) := P(J);

               -- Insert new vertex at edge J (classical 4-point rule)
               declare
                  Prev_Prev_Index : constant Polygon_Index :=
                    (if J-2 < 1 then N + (J-2) else J-2);
                  Prev_Index       : constant Polygon_Index :=
                    (if J-1 < 1 then N else J-1);
                  Next_Index       : constant Polygon_Index :=
                    (if J+1 > N then 1 else J+1);
                  Next_Next_Index  : constant Polygon_Index :=
                    (if J+2 > N then J+2 - N else J+2);

                  -- Exterior angles for S2
                  Delta_Jm1 : constant Real := Exterior_Angle_S2(P, Prev_Prev_Index);
                  Delta_J   : constant Real := Exterior_Angle_S2(P, Prev_Index);
                  Delta_Jp1 : constant Real := Exterior_Angle_S2(P, J);
                  Delta_Jp2 : constant Real := Exterior_Angle_S2(P, Next_Index);

                  -- Classical alpha (Paper Eq. 9)
                  Alpha : constant Real :=
                    0.125 * (Mu * (Delta_Jm1 + Delta_Jp2) + (1.0 - Mu) * (Delta_J + Delta_Jp1));
               begin
                  -- Insert new vertex using angle-based rule (S2 version, Paper Eq. 7)
                  Angle_Based_Insert(P, J, Alpha, P_New(2 * J));
               end;
            end loop;

         when H2 =>
            for J in 1 .. N loop
               -- Retain original vertex
               P_New(2 * J - 1) := P(J);

               -- Insert new vertex at edge J (classical 4-point rule)
               declare
                  Prev_Prev_Index : constant Polygon_Index :=
                    (if J-2 < 1 then N + (J-2) else J-2);
                  Prev_Index       : constant Polygon_Index :=
                    (if J-1 < 1 then N else J-1);
                  Next_Index       : constant Polygon_Index :=
                    (if J+1 > N then 1 else J+1);
                  Next_Next_Index  : constant Polygon_Index :=
                    (if J+2 > N then J+2 - N else J+2);

                  -- Exterior angles for H2
                  Delta_Jm1 : constant Real := Exterior_Angle_H2(P, Prev_Prev_Index);
                  Delta_J   : constant Real := Exterior_Angle_H2(P, Prev_Index);
                  Delta_Jp1 : constant Real := Exterior_Angle_H2(P, J);
                  Delta_Jp2 : constant Real := Exterior_Angle_H2(P, Next_Index);

                  -- Classical alpha (Paper Eq. 9)
                  Alpha : constant Real :=
                    0.125 * (Mu * (Delta_Jm1 + Delta_Jp2) + (1.0 - Mu) * (Delta_J + Delta_Jp1));
               begin
                  -- Insert new vertex using angle-based rule (H2 version, Paper Eq. 8)
                  Angle_Based_Insert(P, J, Alpha, P_New(2 * J));
               end;
            end loop;
      end case;
   end Classical_4Point_Subdivide;

   procedure Angle_Based_Insert (
      P      : in  Polygon(G, N);
      J      : in  Polygon_Index;
      Alpha  : in  Real;
      Q      : out Universal_Point(G))
   is
      Prev_Index, Next_Index : Polygon_Index;
      Clamped_Alpha : constant Real := Clamp_Alpha(Alpha);
   begin
      Prev_Index := (if J = 1 then N else J - 1);
      Next_Index := (if J = N then 1 else J + 1);

      case G is
         when E2 =>
            declare
               T : Point_2D := (P(Next_Index).P2.X - P(J).P2.X,
                               P(Next_Index).P2.Y - P(J).P2.Y);
               Norm_T : constant Real := Euclidean_Distance(T, (0.0, 0.0));
               T_Norm : Point_2D;
               R_T : Point_2D;
               E_New : Real;
            begin
               if Norm_T > 1e-12 then
                  T_Norm := (T.X / Norm_T, T.Y / Norm_T);
               else
                  T_Norm := (0.0, 0.0);
               end if;

               -- Rotation matrix for angle Alpha
               R_T := C_Rotate_Vector(T_Norm, Clamped_Alpha);

               -- New edge length (Paper Eq. 6)
               if abs Clamped_Alpha < Pi/4 - 0.02 then
                  E_New := Sin(abs Clamped_Alpha) / Sin(Pi - 2.0 * abs Clamped_Alpha) * Norm_T;
               else
                  E_New := Norm_T / 2.0;  -- Fallback
               end if;

               Q := (G => E2, P2 => (P(J).P2.X + E_New * R_T.X,
                                      P(J).P2.Y + E_New * R_T.Y));
            end;

         when S2 =>
            declare
               T : Point_3D;
               TX, TY, TZ : Real;
               N : Point_3D;
               E_New : Real;
               Cos_Alpha : constant Real := Cos(Clamped_Alpha);
               Sin_Alpha : constant Real := Sin(Clamped_Alpha);
            begin
               -- Tangent vector from P(J) to P(J+1)
               C_Spherical_Tangent(
                  P(J).P3.X, P(J).P3.Y, P(J).P3.Z,
                  P(Next_Index).P3.X, P(Next_Index).P3.Y, P(Next_Index).P3.Z,
                  TX, TY, TZ
               );
               T := (TX, TY, TZ);

               -- Normal vector (P(J) × T)
               N := (
                  P(J).P3.Y * T.Z - P(J).P3.Z * T.Y,
                  P(J).P3.Z * T.X - P(J).P3.X * T.Z,
                  P(J).P3.X * T.Y - P(J).P3.Y * T.X
               );

               -- New edge length (Paper Eq. 7)
               declare
                  Edge_Length : constant Real := Spherical_Distance(P(J).P3, P(Next_Index).P3);
               begin
                  E_New := Arctan(Tan(Edge_Length / 2.0) / Cos(abs Clamped_Alpha));
               exception
                  when others =>
                     E_New := Edge_Length / 2.0;  -- Fallback
               end;

               -- New vertex (Paper Eq. 7)
               Q := (G => S2, P3 => (
                  Cos(E_New) * P(J).P3.X +
                  Sin(E_New) * (Cos_Alpha * T.X + Sin_Alpha * N.X),

                  Cos(E_New) * P(J).P3.Y +
                  Sin(E_New) * (Cos_Alpha * T.Y + Sin_Alpha * N.Y),

                  Cos(E_New) * P(J).P3.Z +
                  Sin(E_New) * (Cos_Alpha * T.Z + Sin_Alpha * N.Z)
               ));
            end;

         when H2 =>
            declare
               T : Point_2D := C_Hyperbolic_Tangent(P(J).P2, P(Next_Index).P2);
               R_T : Point_2D;
               E_New : Real;
               New_Point : Point_2D;
            begin
               -- Rotate tangent vector by Alpha
               R_T := C_Rotate_Vector(T, Clamped_Alpha);

               -- New edge length (Paper Eq. 8)
               declare
                  Edge_Length : constant Real := Hyperbolic_Distance(P(J).P2, P(Next_Index).P2);
               begin
                  E_New := 2.0 * Arctanh(Tanh(Edge_Length / 2.0) / Cos(abs Clamped_Alpha));
               exception
                  when others =>
                     E_New := Edge_Length / 2.0;  -- Fallback
               end;

               -- Scale the rotated tangent vector
               New_Point := C_Normalize_Vector(R_T);
               New_Point := (New_Point.X * Tanh(E_New / 2.0),
                             New_Point.Y * Tanh(E_New / 2.0));

               -- New vertex (Paper Eq. 8)
               Q := (G => H2, P2 => C_Mobius_Add(P(J).P2, New_Point));

               -- Clamp to Poincaré disk
               Q.P2 := C_Clamp_To_Disk(Q.P2);
            end;
      end case;
   end Angle_Based_Insert;

   -- ===================================================================
   -- SAFETY GUARANTEES (Theorem 1)
   -- ===================================================================

   function Clamp_Alpha (Alpha : Real) return Real is
      Alpha_Min : constant Real := -Pi/4 + 0.02;
      Alpha_Max : constant Real := Pi/4 - 0.02;
   begin
      if Alpha < Alpha_Min then
         return Alpha_Min;
      elsif Alpha > Alpha_Max then
         return Alpha_Max;
      else
         return Alpha;
      end if;
   end Clamp_Alpha;

end Geometry_Subdivision;
