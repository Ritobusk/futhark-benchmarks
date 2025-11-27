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

type QorA = #left     | #right   | #up       | #down      |
           #leftUp   | #rightUp | #leftDown | #rightDown

def findQuodrantOrAxis (p : (i64, i64)) (q : (i64, i64)) (grid_size : i64) : QorA =
    -- The input coordinates of a point will be (y,x) since its the row/col from the grid
    let (dir_x, dir_y) = ((f64.i64 q.1) - (f64.i64 p.1) ,(f64.i64 (grid_size - q.0)) - (f64.i64 (grid_size -p.0)))
    in if dir_x == 0 || dir_y == 0 then
        if      dir_x == 0 && dir_y > 0 then #up
        else if dir_x == 0 && dir_y < 0 then #down
        else if dir_x < 0 && dir_y == 0 then #left
        else #right
    else 
        if      dir_x > 0 && dir_y > 0 then #rightUp
        else if dir_x < 0 && dir_y > 0 then #leftUp
        else if dir_x < 0 && dir_y < 0 then #leftDown
        else #rightDown



def movePointsToGrid [n] (points : [n][2]f32) (grid_size : i64) : ([grid_size][grid_size]i32, [n]i32) =
    -- G1: Moves points to a grid by translating the points such that 0,0 is the center of mass and scaling the points
    --     If a 2 points fit in the same place in the grid, only one of them is inserted.
    let grid      = replicate grid_size (replicate grid_size (-1i32))
    let shiftedP = shiftPoints points

    let pointsT = transpose shiftedP
    let mins = map (reduce_comm f32.min f32.highest) pointsT
    let maxs = map (reduce_comm f32.max f32.lowest ) pointsT 

    let half_grid_size =  ((f32.i64 (grid_size -1) / 2f32) )
    let x_scale = (half_grid_size / (f32.max (f32.abs mins[0]) (f32.abs maxs[0])))
    let y_scale = (half_grid_size / (f32.max (f32.abs mins[1]) (f32.abs maxs[1])))
    let scale   = f32.min x_scale y_scale 
    let scaled_points = map (\cc -> map (\c -> c * scale) cc) pointsT |> transpose

    -- Need to think about this some more. Currently it is sequential with updates to an array.
    in loop (grid' : *[grid_size][grid_size]i32, u_p_f : *[n]i32) = (grid, replicate (n) 0i32) for i < n do 
        let (x, y) = ( f32.round (scaled_points[i][0] + half_grid_size), f32.round (scaled_points[i][1] + half_grid_size))
        let (x, y) = (i64.f32 <| x, i64.f32 <| y )
        in if grid'[y][x] == (-1) then
            (grid' with [y,x] = (i32.i64 i), u_p_f)
           else
            (grid', u_p_f with [i] = i32.i64 1)


def voronoiDiagram [grid_size] (grid : [grid_size][grid_size]i32) : ([grid_size][grid_size](f32, i32, (i64, i64))) =
    -- G2
    -- I use a tuple grid to represent for each pixel both the distance to the closest encountered site, 
    --   this site's index and the coordinates to the original site which is now referenced
    --   by the specific pixel: (dist, site_index, org_coords)
    -- !! The site index is redundant since you can just read it from the plus_1 grid with the 
    --      coordinates that refrences it.
    let stencil_1 = stencilK 1
    let plus_1 = tabulate_2d grid_size grid_size 
        (\r c ->
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
        )

    let test_grid = map (\i -> map (\j -> plus_1[i][j].1) <| iota grid_size) <| iota grid_size 
    let t21 = trace test_grid

    in loop g = plus_1 for iteration < (log2Int grid_size) do
        let stencil_it = stencilK (grid_size / (2**(iteration+1)))
        in tabulate_2d grid_size grid_size 
            (\r c ->
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
            )

def removeIslands [grid_size] (grid : [grid_size][grid_size](f32, i32, (i64, i64))) : [grid_size][grid_size](f32, i32, (i64, i64))=
    --SPØRGSMÅL: Når man berenger det tæteste point til en island pixel, skal man så finde den mindste distance blandt all points
    --           eller kun blandt dens naboer? Og skal man gå ud fra de 'transformerede' grid koordinater eller input koordinater 
    let (g'', _) = 
        loop (g, cond) = (grid, true) while cond do 
            let g'' = tabulate_2d grid_size grid_size 
                (\r c ->
                    let rule = trace <| findQuodrantOrAxis (r,c) g[r][c].2 grid_size
                    let t31  = trace <| (r,c)
                    let t31  = trace <| g[r][c].2
                    let tmp = 
                        if      rule == #up    && g[r-1][c].1 == g[r][c].1 then g[r][c]
                        else if rule == #down  && g[r+1][c].1 == g[r][c].1 then g[r][c]
                        else if rule == #right && g[r][c+1].1 == g[r][c].1 then g[r][c]
                        else if rule == #left  && g[r][c-1].1 == g[r][c].1 then g[r][c]
                        else (-1, -1, (-1,-1)) --Find new site to associate pixel with
                    in tmp


                )
            in (g, false)
    in g''


def main [n]
    (points : [n][2]f32)  =
    let grid_size = trace <| 2 ** (log2Int (i64.f64 <| 3 * (f64.sqrt <| f64.i64 (n) )) + 1) -- To power of 2

    let (grid, unused_p_flag) = movePointsToGrid points grid_size
    let (t4, t10) = trace (grid, unused_p_flag)

    let grid =  voronoiDiagram grid
    let test_grid = map (\i -> map (\j -> grid[i][j].1) <| iota grid_size) <| iota grid_size 
    let t24 = trace test_grid




    in removeIslands grid
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
