--  alpha_table.ads
--  ===============
--  Specification for loading and using pre-trained insertion angles
--  from the neural network as a lookup table in Ada SPARK.
--
--  This package provides:
--  1. Type definitions for the lookup table
--  2. Function to load the binary table from file
--  3. Fallback heuristic if table is not available
--
--  Usage:
--     Alpha_Values : Alpha_Array := Load_Alpha_Table;
--     Alpha : Real := Alpha_Values(Index);

with Geometry_Subdivision; use Geometry_Subdivision;

package Alpha_Table with
   SPARK_Mode => On
is
   -- ===================================================================
   -- TYPE DEFINITIONS
   -- ===================================================================

   -- Size of the lookup table (must match N_EDGES in extract_weights.py)
   Table_Size : constant := 10_000;

   -- Array type for the lookup table
   type Alpha_Array is array (1 .. Table_Size) of Real;

   -- Feature vector type (matches Python extract_weights.py)
   -- [delta_{j-1}/pi, delta_j/pi, delta_{j+1}/pi, delta_{j+2}/pi,
   --  e_j/ebar, e_{j+1}/ebar, kappa]
   type Feature_Vector is array (1 .. 7) of Real;

   -- ===================================================================
   -- LOOKUP TABLE FUNCTIONS
   -- ===================================================================

   -- Load the lookup table from binary file
   function Load_Alpha_Table return Alpha_Array
     with Post => Load_Alpha_Table'Result'Length = Table_Size;

   -- Get alpha value from lookup table using a simple hash of features
   -- Note: This is a placeholder. In a real implementation, you would:
   -- 1. Quantize the feature space
   -- 2. Use a proper hash function
   -- 3. Handle collisions
   function Get_Alpha_From_Table (
      Features : Feature_Vector;
      Table    : Alpha_Array) return Real
     with Pre => Table'Length = Table_Size;

   -- ===================================================================
   -- FALLBACK HEURISTIC (if lookup table is not available)
   -- ===================================================================

   -- Linear adaptive heuristic (Section 7.1 from the paper)
   -- Alpha = Mu* + alpha*(local_curvature - mean_curvature)
   function Predict_Alpha (Features : Feature_Vector) return Real
     with Post => Predict_Alpha(Features) > -Pi/4 + 0.02 and
                  Predict_Alpha(Features) < Pi/4 - 0.02;

   -- ===================================================================
   -- UTILITY FUNCTIONS
   -- ===================================================================

   -- Simple hash function for feature vectors (for demo purposes)
   function Feature_Hash (Features : Feature_Vector) return Integer
     with Post => Feature_Hash(Features) >= 0 and
                  Feature_Hash(Features) < Table_Size;

private
   -- Private type for the hash state (if needed for more complex hashing)
   type Hash_State is range 0 .. Table_Size - 1;

end Alpha_Table;
