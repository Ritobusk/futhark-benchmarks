import "delaunay_triangulation_VD"
import "util"

def voronoiDiagram1 [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64)  : ([grid_size][grid_size](f64, i32)) =
    -- G2
    -- I use a tuple grid to represent for each pixel both the distance to the closest encountered site, 
    --   this site's index and the coordinates to the original site which is now referenced
    --   by the specific pixel: (dist, site_index)
    let stencil_1 = stencilK 1
    let plus_1 = tabulate_2d grid_size grid_size 
        (\r c ->
            if grid[r][c] != -1 then (0f64, grid[r][c])
            else
                loop (d, ind) =  (f64.highest, -1) for (i,j) in stencil_1 do
                    let (r', c') = (r + i, c + j)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (d, ind)
                    else 
                        let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 c',f64.i64 r')
                        let ind' = grid[r'][c']
                        in if ind' != -1 && d' < d then
                            (d', ind')
                        else
                            (d, ind)
        )

    in plus_1

def voronoiDiagram2 [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64)  : ([grid_size][grid_size](f64, i32)) =
    let stencil_1 = stencilK 1
    let plus_1 = tabulate_2d grid_size grid_size 
        (\r c ->
            if grid[r][c] != -1 then (0f64, grid[r][c])
            else
                map (\(si, sj) ->
                    let (r', c') = (r + si, c + sj)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (f64.highest, -1)
                    else
                        let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 c',f64.i64 r')
                        let ind' = grid[r'][c']
                        in if ind' != -1 then
                            (d', ind')
                        else
                            (f64.highest, -1)
                    ) stencil_1
                |> reduce (\acc x -> if x.0 < acc.0 then x else acc ) (f64.highest, -1)
        )

    -- in plus_1   
    in loop g = plus_1 for iteration < (log2Int grid_size) + 1 do
        let stencil_it = trace <|stencilK (grid_size / (2**(iteration+1)))
        in tabulate_2d grid_size grid_size 
            (\r c ->
                map (\(si, sj) ->
                    let (r', c') = (r + si, c + sj)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (f64.highest, -1)
                    else
                        let ind' = grid[r'][c']
                        in if ind' != -1 then
                            let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 points[g[r'][c'].1][0], f64.i64 points[g[r'][c'].1][1])
                            in (d', ind')
                        else
                            (f64.highest, -1)
                    ) stencil_it
                |> reduce (\acc x -> if x.0 < acc.0 then x else acc ) (g[r][c].0, g[r][c].1)
            )

-- ==
-- entry: main test
-- compiled random input {       [100][2]f64 } 
def main [n]
    (points : [n][2]f64)  =
    let grid_size = 16
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size

    let grid2 = voronoiDiagram2  grid  scaled_points_grid
    let grid1 = voronoiDiagram  grid  scaled_points_grid
    let grid2' = map (\r -> map (.1) r) grid2 |> flatten 
    let grid1' = map (\r -> map (.1) r) grid1 |> flatten 
    -- let grid1 = map (.1) grid1 |> flatten 
    let g_flag = map2 (==) grid2' grid1'
        

    in (grid2', grid1')
    -- in (g_flag, grid2)
    -- let grid_size = 1024  --/ 32
    --
    -- let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size
    --
    -- let grid' = voronoiDiagram2  grid  scaled_points_grid
    -- in  grid'

entry test [n]
    (points : [n][2]f64)  =
    let grid_size = 16
    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size

    let grid2 = voronoiDiagram2  grid  scaled_points_grid
    let grid1 = voronoiDiagram  grid  scaled_points_grid
    let grid2 = map (\r -> map (.1) r) grid2 |> flatten 
    let grid1 = map (\r -> map (.1) r) grid1 |> flatten 
    -- let grid1 = map (.1) grid1 |> flatten 
    let g_flag = map2 (==) grid2 grid1
        |> and

    in g_flag 
