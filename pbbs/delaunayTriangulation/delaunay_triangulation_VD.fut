--  futhark run delaunay_triangulation_VD.fut < test_data.txt
--  futhark dataset -g 5i32 -g [10][2]f32 > test_data.txt


--Shifts points' center towards (0,0)
def shiftPoints [n][d]
         (points : [n][d]f32) : [n][d]f32  =
     let mass_acc       = replicate (d) 0.0f32
     let mass_sum       = reduce (\acc point -> map2 (+) acc point ) mass_acc points
     let mass           = map (\i -> i / f32.i64 n) mass_sum
     let shifted_Points = map (\point -> map2 (\pv mv -> pv - mv) point mass) points
     in shifted_Points

def dist (p : (f32, f32)) (q: (f32, f32)) =
    f32.sqrt ((p.0 - q.0)**2 + (p.1 - q.1)**2)

def stencilK (k : i64) =
    [(k,0), (k,k),(k,-k),(0,k),(0,-k),(-k,0),(-k,k),(-k,-k)]

def log2Int (n : i64) : i64 =
  let (_, res) =
    loop (n, r) = (n, 0)
      while n > 1 do
        (n >> 1, r+1)
  in res 

def main [n]
    (grid_size : i64)
    (points : [n][2]f32)  =
    let grid_size = trace <| 2 ** (log2Int (grid_size ) + 1) -- To power of 2
    let grid      = replicate grid_size (replicate grid_size (-1i32))


    --G1
    let shiftedP = shiftPoints points

    let pointsT = transpose shiftedP
    let mins = map (reduce_comm f32.min f32.highest) pointsT
    let maxs = map (reduce_comm f32.max f32.lowest ) pointsT 

    let half_grid_size =  ((f32.i64 (grid_size -1) / 2f32) )
    let x_scale = (half_grid_size / (f32.max (f32.abs mins[0]) (f32.abs maxs[0])))
    let y_scale = (half_grid_size / (f32.max (f32.abs mins[1]) (f32.abs maxs[1])))
    let scale   = f32.min x_scale y_scale 
    let scaled_points = map (\cc -> map (\c -> c * scale) cc) pointsT |> transpose

    -- Need to think about this some more.
    let (grid, unused_p_flag) =
        loop (grid' : *[grid_size][grid_size]i32, u_p_f : *[n]i32) = (grid, replicate (n) 0i32) for i < n do 
            let (x, y) = ( f32.round (scaled_points[i][0] + half_grid_size), f32.round (scaled_points[i][1] + half_grid_size))
            let (t5, t6) = trace (x,y)
            let (x, y) = trace (i64.f32 <| x, i64.f32 <| y )
            in if grid'[y][x] == (-1) then
                (grid' with [y,x] = (i32.i64 i), u_p_f)
               else
                (grid', u_p_f with [i] = i32.i64 1)

    let (t4, t10) = trace (grid, unused_p_flag)

    -- G2
    -- I use a tuple grid to represent for each pixel both the distance to the closest encountered site, 
    --   this site's index and the coordinates to the original site which is now referenced
    --   by the specific pixel: (dist, site_index, org_coords)
    -- !! The site index is redundant since you can just read it from the plus_1 grid with the 
    --      coordinates that refrences it.

    -- Beregn et ny grid med tuples, hvor der er distance og idx
    -- I sekvensielt loop
    --   Beregn for hver pixel parallelt
    --   Beregn baggud, da man kender 'reglen' for, hvilken information der skal propageres.

    let stencil_1 = trace <| stencilK 1
    let plus_1 = 
        map (\r -> 
            map (
                \c -> 
                    if grid[r][c] != -1 then
                        (0f32, grid[r][c], (r,c))
                    else
                        loop (d, ind, (o_r, o_c)) =  (f32.highest, -1, (-1,-1)) for (i,j) in stencil_1 do
                            let (r', c') = (r + i, c + j)
                            in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                                (d, ind, (o_r, o_c))
                            else 
                                let d' = dist (f32.i64 c,f32.i64 r) (f32.i64 c',f32.i64 r')
                                let ind' = grid[r'][c']
                                in if ind' != -1 && d' < d then
                                    (d', ind', (r',c'))
                                else
                                    (d, ind, (o_r, o_c))
            ) (iota grid_size)
        ) <| iota grid_size

    let test_grid = map (\i -> map (\j -> plus_1[i][j].1) <| iota grid_size) <| iota grid_size 
    let t21 = trace plus_1
    let t22 = trace <| log2Int grid_size

    let grid = loop g = plus_1 for iteration < (log2Int grid_size) do
        let stencil_it = trace <| stencilK (grid_size / (2**(iteration+1)))
        let g' = map (\r -> 
            map (\c -> 
                loop (d, ind, (o_r, o_c)) =  (g[r][c].0, g[r][c].1, g[r][c].2) for (i,j) in stencil_it do
                    let (r', c') = (r + i, c + j)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (d, ind, (o_r, o_c))
                    else 
                        let ind' = g[r'][c'].1
                        in if ind' != -1  then
                            let d' = dist (f32.i64 c,f32.i64 r) (f32.i64 g[r'][c'].2.1, f32.i64 g[r'][c'].2.0)
                            in if d' < d then
                                (d', ind', g[r'][c'].2)
                            else
                                (d, ind, (o_r, o_c))      
                        else
                            (d, ind, (o_r, o_c))
            ) (iota grid_size)
        ) (iota grid_size)
        let test_grid = map (\i -> map (\j -> g'[i][j].1) <| iota grid_size) <| iota grid_size 
        let t24 = trace g'
        let t23 = trace test_grid
        in g'

    in unused_p_flag
-- In:
-- trace: ([[-1, -1, -1, -1, -1, -1, -1, -1],
--          [ 1, -1, -1, -1, -1, -1, -1, -1],
--          [-1, -1, -1,  0, -1,  8, -1, -1],
--          [-1, -1, -1,  4, -1,  2,  7, -1],
--          [-1, -1,  5, -1, -1,  9, -1, -1],
--          [-1, -1, -1, -1, -1, -1, -1, -1],
--          [-1,  3, -1, -1, -1, -1, -1, -1],
--          [-1, -1, -1, -1, -1, -1,  6, -1]],
--Result:
-- trace: [[1, 1, 0, 0, 8, 8, 8, 8],
--         [1, 1, 0, 0, 8, 8, 8, 7],
--         [1, 1, 0, 0, 8, 8, 7, 7],
--         [1, 5, 5, 4, 2, 2, 7, 7],
--         [5, 5, 5, 5, 9, 9, 9, 7],
--         [3, 3, 5, 5, 9, 9, 9, 6],
--         [3, 3, 3, 3, 6, 6, 6, 6],
--         [3, 3, 3, 3, 6, 6, 6, 6]]
