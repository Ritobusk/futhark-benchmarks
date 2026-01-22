import "util"
import "delaunay_triangulation_VD"

-- > :img GPU_STEPS ($loaddata "test_data10.txt")

-- > :img GPU_STEPS_AND_C1 ($loaddata "test_data10.txt")

def GPU_STEPS [n]
    (points : [n][2]f64)  =
    let grid_size = 512 

    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    -- let t =trace scaled_points
    -- let (t4, t10) = trace (grid, unused_p_flag)

    let grid' = voronoiDiagram2  grid  scaled_points_grid
    let grid''  = removeIslands grid' scaled_points_grid
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    -- G5 and G6
    let triangles = createTriangles voronoi_diagram voronoi_vertices

    -- in voronoiAndPoints grid voronoi_diagram scaled_points_grid
    in triangleGrid grid voronoi_diagram triangles scaled_points_grid

def GPU_STEPS_AND_C1 [n]
    (points : [n][2]f64)  =
    let grid_size = 512

    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size

    let grid' =  voronoiDiagram grid scaled_points_grid
    let grid' =  removeIslands grid' scaled_points_grid
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid'[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    let triangles = createTriangles voronoi_diagram voronoi_vertices

    let triangles = (fixConvexHull voronoi_diagram scaled_points_grid) ++ triangles

    -- in voronoiAndPoints grid voronoi_diagram scaled_points_grid
    in triangleGrid grid voronoi_diagram triangles scaled_points_grid
