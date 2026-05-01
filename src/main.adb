--  main.adb
--  ==========
--  Main test program for Ada SPARK implementation of curve subdivision.
--  Demonstrates usage of geometry_core and math_utils packages.
--
--  Compile with: gnatmake -P ../gnatprove/geometry_subdivision.gpr
--  Run with:     ./obj/main

with Geometry_Subdivision; use Geometry_Subdivision;
with Math_Utils;            use Math_Utils;
with Alpha_Table;           use Alpha_Table;
with Ada.Text_IO;           use Ada.Text_IO;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;

procedure Main is
   -- ===================================================================
   -- TEST 1: Basic Geometry Operations (E2, S2, H2)
   -- ===================================================================

   procedure Test_Geometry_Operations is
      -- Euclidean (E2)
      P1, P2 : Point_2D := (0.0, 0.0);
      Dist_E2 : Real;

      -- Spherical (S2)
      Q1, Q2 : Point_3D := (1.0, 0.0, 0.0);
      Dist_S2 : Real;

      -- Hyperbolic (H2)
      R1, R2 : Point_2D := (0.0, 0.0);
      Dist_H2 : Real;
      Tangent_H2 : Point_2D;
   begin
      Put_Line("=== TEST 1: Geometry Operations ===");

      -- E2: Euclidean distance
      P2 := (3.0, 4.0);
      Dist_E2 := Euclidean_Distance(P1, P2);
      Put_Line("E2 Distance between (0,0) and (3,4):" & Dist_E2'Image);

      -- S2: Spherical distance (points on unit sphere)
      Q2 := (0.0, 1.0, 0.0);
      Dist_S2 := Spherical_Distance(Q1, Q2);
      Put_Line("S2 Distance between (1,0,0) and (0,1,0):" & Dist_S2'Image);

      -- H2: Hyperbolic distance
      R2 := (0.5, 0.3);
      Dist_H2 := C_Hyperbolic_Distance(R1, R2);
      Put_Line("H2 Distance between (0,0) and (0.5,0.3):" & Dist_H2'Image);

      -- H2: Hyperbolic tangent
      Tangent_H2 := C_Hyperbolic_Tangent(R1, R2);
      Put_Line("H2 Tangent from (0,0) to (0.5,0.3): (" &
              Tangent_H2.X'Image & "," & Tangent_H2.Y'Image & ")");

      New_Line;
   end Test_Geometry_Operations;

   -- ===================================================================
   -- TEST 2: Exterior Angles (E2, S2, H2)
   -- ===================================================================

   procedure Test_Exterior_Angles is
      -- E2: Square polygon
      type Square_Polygon is new Polygon(E2, 4);
      Square : Square_Polygon := (
        (G => E2, P2 => (0.0, 0.0)),
        (G => E2, P2 => (1.0, 0.0)),
        (G => E2, P2 => (1.0, 1.0)),
        (G => E2, P2 => (0.0, 1.0))
      );
      Angle_E2 : Real;

      -- S2: Tetrahedron vertices (approximate)
      type Tetra_Polygon is new Polygon(S2, 4);
      Tetra : Tetra_Polygon := (
        (G => S2, P3 => (1.0, 0.0, 0.0)),
        (G => S2, P3 => (0.0, 1.0, 0.0)),
        (G => S2, P3 => (0.0, 0.0, 1.0)),
        (G => S2, P3 => (0.0, 0.0, -1.0))
      );
      Angle_S2 : Real;

      -- H2: Small polygon near origin
      type H2_Polygon is new Polygon(H2, 4);
      H2_Poly : H2_Polygon := (
        (G => H2, P2 => (0.1, 0.1)),
        (G => H2, P2 => (0.2, 0.1)),
        (G => H2, P2 => (0.2, 0.2)),
        (G => H2, P2 => (0.1, 0.2))
      );
      Angle_H2 : Real;
   begin
      Put_Line("=== TEST 2: Exterior Angles ===");

      -- E2 exterior angle at vertex 1
      Angle_E2 := Exterior_Angle_E2(Square, 1);
      Put_Line("E2 Exterior angle at vertex 1 (square):" & Angle_E2'Image);

      -- S2 exterior angle at vertex 1
      Angle_S2 := Exterior_Angle_S2(Tetra, 1);
      Put_Line("S2 Exterior angle at vertex 1 (tetrahedron):" & Angle_S2'Image);

      -- H2 exterior angle at vertex 1
      Angle_H2 := Exterior_Angle_H2(H2_Poly, 1);
      Put_Line("H2 Exterior angle at vertex 1 (small polygon):" & Angle_H2'Image);

      New_Line;
   end Test_Exterior_Angles;

   -- ===================================================================
   -- TEST 3: Classical 4-Point Subdivision (E2)
   -- ===================================================================

   procedure Test_Classical_Subdivision is
      type Square_Polygon is new Polygon(E2, 4);
      type Refined_Polygon is new Polygon(E2, 8);

      Input : Square_Polygon := (
        (G => E2, P2 => (0.0, 0.0)),
        (G => E2, P2 => (1.0, 0.0)),
        (G => E2, P2 => (1.0, 1.0)),
        (G => E2, P2 => (0.0, 1.0))
      );
      Output : Refined_Polygon;
   begin
      Put_Line("=== TEST 3: Classical 4-Point Subdivision (E2) ===");
      Put_Line("Input polygon (square):");
      for I in 1 .. 4 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Input(I).P2.X'Image & "," & Input(I).P2.Y'Image & ")");
      end loop;

      -- Apply classical subdivision
      Classical_4Point_Subdivide(Input, Output);

      Put_Line("Output polygon (after subdivision):");
      for I in 1 .. 8 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Output(I).P2.X'Image & "," & Output(I).P2.Y'Image & ")");
      end loop;
      New_Line;
   end Test_Classical_Subdivision;

   -- ===================================================================
   -- TEST 4: Angle-Based Insertion (E2)
   -- ===================================================================

   procedure Test_Angle_Based_Insertion is
      type Triangle_Polygon is new Polygon(E2, 3);
      type Refined_Polygon is new Polygon(E2, 6);

      Input : Triangle_Polygon := (
        (G => E2, P2 => (0.0, 0.0)),
        (G => E2, P2 => (1.0, 0.0)),
        (G => E2, P2 => (0.5, 1.0))
      );
      Output : Refined_Polygon;
      Alpha : Real := 0.1;  -- Test angle
      Q : Universal_Point(E2);
   begin
      Put_Line("=== TEST 4: Angle-Based Insertion (E2) ===");
      Put_Line("Input polygon (triangle):");
      for I in 1 .. 3 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Input(I).P2.X'Image & "," & Input(I).P2.Y'Image & ")");
      end loop;

      -- Insert a vertex at edge 1 with angle Alpha
      Angle_Based_Insert(Input, 1, Alpha, Q);
      Put_Line("Inserted vertex at edge 1 with alpha =" & Alpha'Image & ": (" &
               Q.P2.X'Image & "," & Q.P2.Y'Image & ")");

      -- Apply full subdivision with angle-based insertion
      Classical_4Point_Subdivide(Input, Output, Mu => 0.0);

      Put_Line("Output polygon (after subdivision):");
      for I in 1 .. 6 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Output(I).P2.X'Image & "," & Output(I).P2.Y'Image & ")");
      end loop;
      New_Line;
   end Test_Angle_Based_Insertion;

   -- ===================================================================
   -- TEST 5: Lookup Table Usage
   -- ===================================================================

   procedure Test_Lookup_Table is
      Alpha_Values : Alpha_Array;
      Features : Feature_Vector := (0.1, 0.2, -0.1, 0.05, 1.0, 0.9, 0.0);  -- Example features
      Alpha : Real;
   begin
      Put_Line("=== TEST 5: Lookup Table Usage ===");

      -- Load lookup table
      begin
         Alpha_Values := Load_Alpha_Table;
         Put_Line("Loaded lookup table with" & Alpha_Values'Length'Image & " entries.");

         -- Use a simple hash to select an index (for demo purposes)
         -- In a real implementation, you'd use a proper feature-to-index mapping
         declare
            Index : Integer := 1 + (abs Integer(Features(1) * 1000)) mod Alpha_Values'Length;
         begin
            Alpha := Alpha_Values(Index);
            Put_Line("Selected alpha from lookup table:" & Alpha'Image);
         end;

      exception
         when others =>
            Put_Line("Warning: Could not load lookup table. Using heuristic instead.");
            Alpha := Predict_Alpha(Features);  -- Fallback to heuristic
      end;

      -- Clamp alpha to G1-safe range
      Alpha := Clamp_Alpha(Alpha);
      Put_Line("Clamped alpha:" & Alpha'Image);
      New_Line;
   end Test_Lookup_Table;

   -- ===================================================================
   -- TEST 6: G1 Safety Bound (Theorem 1)
   -- ===================================================================

   procedure Test_G1_Safety is
      Test_Alphas : array (1 .. 5) of Real := (
        -1.0, -0.785, 0.0, 0.785, 1.0  -- Includes values outside safe range
      );
      Clamped_Alpha : Real;
   begin
      Put_Line("=== TEST 6: G1 Safety Bound (Theorem 1) ===");
      for I in Test_Alphas'Range loop
         Clamped_Alpha := Clamp_Alpha(Test_Alphas(I));
         Put_Line("Alpha =" & Test_Alphas(I)'Image & " -> Clamped to" &
                  Clamped_Alpha'Image);
      end loop;
      New_Line;
   end Test_G1_Safety;

   -- ===================================================================
   -- TEST 7: Hyperbolic Geometry (H2)
   -- ===================================================================

   procedure Test_Hyperbolic_Geometry is
      type H2_Polygon is new Polygon(H2, 4);
      Input : H2_Polygon := (
        (G => H2, P2 => (0.1, 0.1)),
        (G => H2, P2 => (0.2, 0.1)),
        (G => H2, P2 => (0.2, 0.2)),
        (G => H2, P2 => (0.1, 0.2))
      );
      Output : Polygon(H2, 8);
      Dist : Real;
   begin
      Put_Line("=== TEST 7: Hyperbolic Geometry (H2) ===");

      -- Test hyperbolic distance
      Dist := C_Hyperbolic_Distance(Input(1).P2, Input(2).P2);
      Put_Line("H2 Distance between vertices 1 and 2:" & Dist'Image);

      -- Test Möbius addition
      declare
         Z, W, Result : Point_2D;
      begin
         Z := Input(1).P2;
         W := Input(2).P2;
         Result := C_Mobius_Add(Z, W);
         Put_Line("Möbius addition of (" & Z.X'Image & "," & Z.Y'Image & ") and (" &
                  W.X'Image & "," & W.Y'Image & ") = (" &
                  Result.X'Image & "," & Result.Y'Image & ")");
      end;

      -- Test classical subdivision in H2
      Put_Line("Input polygon (H2):");
      for I in 1 .. 4 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Input(I).P2.X'Image & "," & Input(I).P2.Y'Image & ")");
      end loop;

      Classical_4Point_Subdivide(Input, Output);
      Put_Line("Output polygon (after subdivision in H2):");
      for I in 1 .. 8 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Output(I).P2.X'Image & "," & Output(I).P2.Y'Image & ")");
      end loop;
      New_Line;
   end Test_Hyperbolic_Geometry;

   -- ===================================================================
   -- TEST 8: Spherical Geometry (S2)
   -- ===================================================================

   procedure Test_Spherical_Geometry is
      type S2_Polygon is new Polygon(S2, 4);
      Input : S2_Polygon := (
        (G => S2, P3 => (1.0, 0.0, 0.0)),
        (G => S2, P3 => (0.0, 1.0, 0.0)),
        (G => S2, P3 => (0.0, 0.0, 1.0)),
        (G => S2, P3 => (0.0, 0.0, -1.0))
      );
      Output : Polygon(S2, 8);
      Dist : Real;
      TX, TY, TZ : Real;
   begin
      Put_Line("=== TEST 8: Spherical Geometry (S2) ===");

      -- Test spherical distance
      Dist := C_Spherical_Distance(
         Input(1).P3.X, Input(1).P3.Y, Input(1).P3.Z,
         Input(2).P3.X, Input(2).P3.Y, Input(2).P3.Z
      );
      Put_Line("S2 Distance between vertices 1 and 2:" & Dist'Image);

      -- Test spherical tangent
      C_Spherical_Tangent(
         Input(1).P3.X, Input(1).P3.Y, Input(1).P3.Z,
         Input(2).P3.X, Input(2).P3.Y, Input(2).P3.Z,
         TX, TY, TZ
      );
      Put_Line("Spherical tangent from vertex 1 to 2: (" &
               TX'Image & "," & TY'Image & "," & TZ'Image & ")");

      -- Test classical subdivision in S2
      Put_Line("Input polygon (S2):");
      for I in 1 .. 4 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Input(I).P3.X'Image & "," & Input(I).P3.Y'Image & "," &
                  Input(I).P3.Z'Image & ")");
      end loop;

      Classical_4Point_Subdivide(Input, Output);
      Put_Line("Output polygon (after subdivision in S2):");
      for I in 1 .. 8 loop
         Put_Line("  Vertex" & I'Image & ": (" &
                  Output(I).P3.X'Image & "," & Output(I).P3.Y'Image & "," &
                  Output(I).P3.Z'Image & ")");
      end loop;
      New_Line;
   end Test_Spherical_Geometry;

begin
   -- Run all tests
   Put_Line("=== Ada SPARK Curve Subdivision Tests ===");
   New_Line;

   Test_Geometry_Operations;
   Test_Exterior_Angles;
   Test_Classical_Subdivision;
   Test_Angle_Based_Insertion;
   Test_Lookup_Table;
   Test_G1_Safety;
   Test_Hyperbolic_Geometry;
   Test_Spherical_Geometry;

   Put_Line("=== All tests completed successfully! ===");
exception
   when others =>
      Put_Line("Error: An exception occurred during testing.");
      raise;
end Main;
