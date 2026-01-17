import "delaunay_triangulation_VD"

-- ==
-- entry: g1 g2 g3 g4 g5
-- compiled random input {       [1000][2]f64 } 
-- compiled random input {      [10000][2]f64 } 
-- compiled random input {    [1000000][2]f64 } 
-- compiled random input {   [10000000][2]f64 } 
entry g1 [n]
    (points : [n][2]f64)  =
    let grid_size = 4096
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    in (grid, used, unused, scaled_points, scaled_points_grid)

entry g2 [n]
    (points : [n][2]f64)  =
    let grid_size = 4096
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    let grid' = voronoiDiagram  grid  scaled_points_grid
    in grid'

entry g3 [n]
    (points : [n][2]f64)  =
    let grid_size = 4096
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    let grid' = voronoiDiagram  grid  scaled_points_grid
    let grid''  = removeIslands grid' scaled_points_grid
    in grid''

entry g4 [n]
    (points : [n][2]f64)  =
    let grid_size = 4096
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    let grid' = voronoiDiagram  grid  scaled_points_grid
    let grid''  = removeIslands grid' scaled_points_grid
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 
    in voronoi_vertices 

entry g5 [n]
    (points : [n][2]f64)  =
    let grid_size = 4096
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    let grid' = voronoiDiagram  grid  scaled_points_grid
    let grid''  = removeIslands grid' scaled_points_grid
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 
    let triangles = createTriangles voronoi_diagram voronoi_vertices
    in triangles 

-- def main [n]
--     (points : [n][2]f64)  =
--     let grid_size = 4096  --/ 32
--
--     let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
--
--     let grid' = voronoiDiagram  grid  scaled_points_grid
--     let grid''  = removeIslands grid' scaled_points_grid
--     let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
--     let voronoi_vertices = locateVoronoiVertices voronoi_diagram 
--
--     -- G5 and G6
--     let triangles = createTriangles voronoi_diagram voronoi_vertices
--
--     let triangles = (fixConvexHull voronoi_diagram scaled_points_grid) ++ triangles
--
--     -- in length triangles 
--     --in map (\i -> [triangles[i].0, triangles[i].1,triangles[i].2]) <| indices triangles
--     in  (shiftSites triangles used scaled_points scaled_points_grid)
--     -- in reduce (i64.max) 0 shps
