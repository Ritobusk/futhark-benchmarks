import "delaunay_triangulation_VD"


-- ==
-- entry: main main2 test
-- compiled random input {       [1000][2]f32 } 
-- compiled random input {    [1000000][2]f32 } 
-- compiled random input {   [10000000][2]f32 } 
def main [n]
    (points : [n][2]f32)  =
    let grid_size = 1024
    let (grid, unused, scaled_points) = movePointsToGrid1 points grid_size
    let (grid', unused', scaled_points') = movePointsToGrid2 points grid_size
    in grid

entry main2 [n]
    (points : [n][2]f32)  =
    let grid_size = 1024

    let (grid', unused', scaled_points') = movePointsToGrid2 points grid_size
    in grid'

entry test [n]
    (points : [n][2]f32)  =
    let grid_size = 1024
    let (grid, unused, scaled_points) = movePointsToGrid1 points grid_size
    let (grid', unused', scaled_points') = movePointsToGrid2 points grid_size
    let grid = flatten grid
    let grid' = flatten grid'
    let g_flag = map2 (==) grid grid'
        |> and

    let u_flag = 
        if (length unused) == (length unused') then
            map2 (==) (sized (length unused) unused) (sized (length unused) unused')
            |> and
        else false
    let p_flag = map2 (==) scaled_points scaled_points'
        |> and
    in g_flag && u_flag && p_flag
