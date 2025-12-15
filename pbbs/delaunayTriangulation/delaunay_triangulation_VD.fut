--  futhark run delaunay_triangulation_VD.fut < test_data.txt
--  futhark dataset -g 5i32 -g [10][2]f32 > test_data.txt

import "util"
import "lib/github.com/diku-dk/sorts/radix_sort"

def movePointsToGrid [n] (points : [n][2]f32) (grid_size : i64) : ([grid_size][grid_size]i32, *[]i32) =
    -- G1: Moves points to a grid by translating the points such that 0,0 is the center of mass and scaling the points
    --     If a 2 points fit in the same place in the grid, only one of them is inserted.
    let grid      = replicate grid_size (replicate grid_size (-1i32))
    let shiftedP = shiftPoints points

    let pointsT = transpose shiftedP
    let mins = map (reduce_comm f32.min f32.highest) pointsT
    let maxs = map (reduce_comm f32.max f32.lowest ) pointsT 

    let half_grid_size =  (f32.i64 (grid_size -1) / 2f32) 
    let x_scale = (half_grid_size / (f32.max (f32.abs mins[0]) (f32.abs maxs[0])))
    let y_scale = (half_grid_size / (f32.max (f32.abs mins[1]) (f32.abs maxs[1])))
    let scale   = f32.min x_scale y_scale 

    let scaled_points = map (\cc -> map (\c -> i64.f32 <| f32.round <| c * scale + half_grid_size) cc) pointsT |> transpose
    let scaled_points_flat_idx = map (\cc -> 
        let x = (i64.f32 <| f32.round <| cc[0] * scale + half_grid_size)   
        let y = grid_size * (i64.f32 <| f32.round <| cc[1] * scale + half_grid_size)
        in x + y
    ) <| transpose pointsT
    let scaled_points_flat_idx = scaled_points_flat_idx ++ [99]
    -- Jeg forstår ikke radix sort...
    let sorted_ids = radix_sort_int_by_key (\k -> k.0) (i32.i64 <| log2Int (grid_size**2)) i64.get_bit (zip scaled_points_flat_idx (iota (n+1)))
    let t11 = trace <| scaled_points_flat_idx
    let t14 = trace <| sorted_ids
    let (used, unused) = pack_points sorted_ids

    let t15 = trace <| used
    let t16 = trace <| unused

    let grid =
        scatter (flatten grid) (map (.0) used) (map (\x -> i32.i64 x.1) used)
        |> unflatten
    in (grid, (map (\x -> i32.i64 x.1) unused))


def voronoiDiagram [grid_size] (grid : [grid_size][grid_size]i32) : ([grid_size][grid_size](f32, i32, (i64, i64))) =
    -- G2
    -- I use a tuple grid to represent for each pixel both the distance to the closest encountered site, 
    --   this site's index and the coordinates to the original site which is now referenced
    --   by the specific pixel: (dist, site_index, org_coords)
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
    let stencil_1 = stencilK 1
    let (g'', _) = 
        loop (g, cond) = (grid, true) while cond do 
            let (g'', island_flag_array) =  unzip <| flatten <| tabulate_2d grid_size grid_size 
                (\r c ->
                    let rule = findQuodrantOrAxis (r,c) g[r][c].2 grid_size
                    let tmp = 
                        if      rule == #center then (g[r][c], false)
                        else if rule == #up    && g[r-1][c].1 == g[r][c].1 then (g[r][c], false)
                        else if rule == #down  && g[r+1][c].1 == g[r][c].1 then (g[r][c], false)
                        else if rule == #right && g[r][c+1].1 == g[r][c].1 then (g[r][c], false)
                        else if rule == #left  && g[r][c-1].1 == g[r][c].1 then (g[r][c], false)
                        else if rule == #leftUp    && (checkQuadrant g[r][c].1 g[r-1][c-1].1 g[r-1][c].1 g[r][c-1].1) then (g[r][c], false)
                        else if rule == #leftDown  && (checkQuadrant g[r][c].1 g[r+1][c-1].1 g[r+1][c].1 g[r][c-1].1) then (g[r][c], false)
                        else if rule == #rightUp && (checkQuadrant g[r][c].1 g[r-1][c+1].1 g[r-1][c].1 g[r][c+1].1) then (g[r][c], false)
                        else if rule == #rightDown  && (checkQuadrant g[r][c].1 g[r+1][c+1].1 g[r+1][c].1 g[r][c+1].1) then (g[r][c], false)

                        --Find new site to associate pixel with
                        else 
                            -- Find the closest point among the neighbours
                            let t33 = trace (rule, r,c, g[r][c])
                            in loop ((d, ind, (o_r, o_c)), island_flag) =  ((f32.highest, -1, (-1,-1)), true) for (i,j) in stencil_1 do
                                let (r', c') = (r + i, c + j)
                                in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                                    ((d, ind, (o_r, o_c)), island_flag)
                                else 
                                    let ind' = g[r'][c'].1
                                    in if ind' != -1  then
                                        let d' = dist (f32.i64 c,f32.i64 r) (f32.i64 g[r'][c'].2.1, f32.i64 g[r'][c'].2.0)
                                        in if d' < d then
                                            ((d', ind', g[r'][c'].2), island_flag)
                                        else
                                            ((d, ind, (o_r, o_c)), island_flag)
                                    else
                                        ((d, ind, (o_r, o_c)), island_flag)
                    in tmp
                )
            -- Check if any islands where found. If so we loop again!
            let cond' = trace <| reduce (\acc f -> f) false (island_flag_array)
            in (unflatten g'', cond')
    in g''

def locateVoronoiVertices [grid_size] (grid : [grid_size][grid_size]i32) : [grid_size][grid_size]i32 = --: [m](i64, i64) =
    -- tabulate_2d (grid_size - 2) (grid_size - 2) 
    tabulate_2d (grid_size ) (grid_size) 
        ( \r c ->
            if r == 0 || c == 0 || r == (grid_size -1) || c == (grid_size - 1) then 0
            else
                let is_corner = classifyVertexI32 grid[r][c+1] grid[r][c] grid[r+1][c] grid[r+1][c+1]
                in is_corner
        )

def locateVoronoiVerticesAndCreateTriangulation [grid_size] (grid : [grid_size][grid_size]i32) : [grid_size][grid_size]((i32, i32, i32),(i32, i32, i32)) = --: [m](i64, i64) =
    -- tabulate_2d (grid_size - 2) (grid_size - 2) 
    let trs =tabulate_2d (grid_size ) (grid_size) 
        ( \r c ->
            if r == 0 || c == 0 || r == (grid_size -1) || c == (grid_size - 1) then ((-1,-1,-1),(-1,-1,-1))
            else
                
                 classifyVertexAndTriangulation grid[r][c+1] grid[r][c] grid[r+1][c] grid[r+1][c+1]
        )

    -- let trs = filter  (\t -> t.0.0 > -1)  <| trs
    in trs
-- G5:
-- Another sgmscan on an empty grid to read the last entry of each segment

--G5 - G6
-- kan man ikke bare gøre alt dette i G4?
--  Når intet findes returnerer man (0, (-1,-1,-1))
--  Når en 3'er findes returnerer man (1, (c1,c2,c3))
--  Når en 4'er findes returnerer man (2, (c1,c2,c3,c4))
-- Og til sidst en filter der fjerner, dem der ikke er trekanter. 

-- 1. prøv scan så man får et array, hvor alle voronoiVertices er indexeret.
-- 2. Derefter alloker et array til trekanterne
-- 3. lav map og scatter til at tilføje trekanterne til 2.


-- > :img main ($loaddata "test_data.txt")
-- gridToGray (tabulate_2d grid_size grid_size (\i j -> i32.bool voronoi_vertices[i][j])) (1)
-- let test_grid' = colours (tabulate_2d grid_size grid_size (\i j -> grid'[i][j].1)) --(i32.i64 <| n-1)

def main [n]
    (points : [n][2]f32)  =
    -- grid is: total_grid_size <= 18n, i.e. O(n)
    -- let grid_size = trace <| 2 ** (log2Int (i64.f64 <| 3 * (f64.sqrt <| f64.i64 (n) )) + 1) -- To power of 2
    let grid_size = 16

    let (grid, unused_p_flag) = movePointsToGrid points grid_size
    let (t4, t10) = trace (grid, unused_p_flag)

    let grid' =  voronoiDiagram grid
    let test_grid = trace <| map (\i -> map (\j -> grid'[i][j].1) <| iota grid_size) <| iota grid_size 
    let test_grid' = colours (tabulate_2d grid_size grid_size (\i j -> grid'[i][j].1)) --(i32.i64 <| n-1)
    let voronoi_vertices = trace <| locateVoronoiVertices test_grid
    let fvv = flatten voronoi_vertices
    let fvv_ids = trace <| scan (+) 0i32 fvv
    let num_ts = trace <| last fvv_ids
    let num_ts' = i64.i32 <| num_ts *2
    let triangles = replicate (num_ts') [-1, -1, -1] 
    let vv_idxs = map (\x -> if x.1 > 0 then (x.0, x.2) else (-1, 0) ) <| zip3 fvv_ids fvv (indices fvv)
    let vv_idxs' = filter (\x -> if x.0 < 0 then false else true) vv_idxs
    let vv_idxs' = trace vv_idxs'
    let a =  filter (\x -> x.0 >= 0) <| map (\i -> 
        let j = vv_idxs'[i/2].1
        let c = j % grid_size
        let r = j / grid_size 
        let t = classifyVertexAndTriangulation test_grid[r][c+1] test_grid[r][c] test_grid[r+1][c] test_grid[r+1][c+1]
        in if i%2==0 then t.1 else t.0

        ) (iota num_ts')

    let t6 = trace <| a

    let grid'' = removeIslands grid'

    in gridToGray (voronoi_vertices) (1)


-- Comments
-- grid_size burde måske ikke afhænge af 'n', da det kan gøre noget ved den asymptotiske køretid.
-- Brug https://futhark-lang.org/examples/removing-duplicates.html   til at få G1 til at være parallel
-- Brug https://futhark-lang.org/examples/literate-basics.html       til at visualiserer gridet
-- Overvej om griddet har brug for de 3 tupler, som bliver lavet i G2 eller om man kan nøjes med kun idx
--    fremfor både idx og org_coords. Man kan nemlig bruge idx til at læse fra et n-langt array med org_coords
--    Vent til at jeg er sidst i processen til at se om det giver en lille speedup eller ej.


