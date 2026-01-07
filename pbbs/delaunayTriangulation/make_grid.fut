import "util"
import "delaunay_triangulation_VD"

-- > :img GPU_STEPS ($loaddata "test_data200.txt")

-- > :img GPU_STEPS_AND_C1 ($loaddata "test_data200.txt")

def GPU_STEPS [n]
    (points : [n][2]f32)  =
    let grid_size = 1024 

    let (grid, unused, scaled_points) = movePointsToGrid points grid_size
    -- let t =trace scaled_points
    -- let (t4, t10) = trace (grid, unused_p_flag)

    let grid' = voronoiDiagram  grid  scaled_points
    let grid''  = removeIslands grid' scaled_points
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    -- G5 and G6
    let triangles = createTriangles voronoi_diagram voronoi_vertices

    in triangleGrid grid voronoi_diagram triangles scaled_points

def GPU_STEPS_AND_C1 [n]
    (points : [n][2]f32)  =
    let grid_size = 1024

    let (grid, unused, scaled_points) = movePointsToGrid points grid_size

    let grid' =  voronoiDiagram grid scaled_points
    let grid' =  removeIslands grid' scaled_points
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid'[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    let triangles = createTriangles voronoi_diagram voronoi_vertices

    let triangles = (fixConvexHull voronoi_diagram scaled_points) ++ triangles

    in triangleGrid grid voronoi_diagram triangles scaled_points
