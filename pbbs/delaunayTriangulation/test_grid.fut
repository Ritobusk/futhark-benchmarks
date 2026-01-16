import "lib/github.com/diku-dk/sorts/radix_sort"
import "delaunay_triangulation_VD"
import "util"

def movePointsToGrid1 [n] (points : [n][2]f64) (grid_size : i64) : ([grid_size][grid_size]i32, []i64, [n][2]i64) =

    -- G1: Moves points to a grid by translating the points such that 0,0 is the center of mass and scaling the points
    --     If a 2 points fit in the same place in the grid, only one of them is inserted.

    let pointsT = transpose points  
    let mins = map (reduce_comm f64.min f64.highest) pointsT
    let maxs = map (reduce_comm f64.max f64.lowest ) pointsT 
    let largest_min = f64.minimum mins -- Pushes the points towards the first quodrant

    let x_scale = ((f64.i64 grid_size - 1) / (f64.max (f64.abs mins[0]) (f64.abs maxs[0])))
    let y_scale = ((f64.i64 grid_size - 1) / (f64.max (f64.abs mins[1]) (f64.abs maxs[1])))
    let scale   = f64.min x_scale y_scale 

    let scaled_points = map (\cc -> map (\c -> i64.f64 <| f64.round <| (c - largest_min) * scale ) cc) <| transpose pointsT 

    let scaled_points_flat_idx =  map (\cc -> 
        let x = (i64.f64 <| f64.round <| (cc[0] - largest_min) * scale )   
        let y = grid_size * (i64.f64 <| f64.round <| (cc[1] - largest_min) * scale )
        in x + y
    ) <| transpose pointsT

    -- let (grid, unused) = populateGrid (grid_size * grid_size) scaled_points_flat_idx (indices points)
    -- let grid = unflatten grid
    let sorted_ids = radix_sort_int_by_key (\k -> k.0) (i32.i64 <| (log2Int (grid_size**2 ) +1)) i64.get_bit (zip scaled_points_flat_idx (iota (n)))

    let (used, unused) = pack_points_i64 sorted_ids
    let grid      = replicate grid_size (replicate grid_size (-1i32))
    let grid =
        scatter (flatten grid) (map (.0) used) (map (\x -> i32.i64 x.1) used)
        |> unflatten

    in (grid, (map (\x -> x.1) unused), scaled_points)

-- ==
-- entry: main main2 test
-- compiled random input {       [1000][2]f64 } 
-- compiled random input {    [1000000][2]f64 } 
-- compiled random input {   [10000000][2]f64 } 
def main [n]
    (points : [n][2]f64)  =
    let grid_size = 1024
    let (grid, unused, scaled_points) = movePointsToGrid1 points grid_size
    in grid

entry main2 [n]
    (points : [n][2]f64)  =
    let grid_size = 1024

    let (grid', used, unused', scaled_points', spg)  = movePointsToGrid points grid_size
    in grid'

entry test [n]
    (points : [n][2]f64)  =
    let grid_size = 1024
    let (grid, unused, scaled_points) = movePointsToGrid1 points grid_size
    let (grid', used, unused', sp, scaled_points') = movePointsToGrid points grid_size
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
