--  futhark run delaunay_triangulation_VD.fut < test_data.txt
--  futhark dataset -g 5i32 -g [10][2]f64 > test_data.txt

import "util"
import "lib/github.com/diku-dk/sorts/radix_sort"

def movePointsToGrid [n] (points : [n][2]f64) (grid_size : i64) : ([grid_size][grid_size]i32, []i64, []i64, [n][2]f64, [n][2]i64) =

    -- G1: Moves points to a grid by translating the points such that 0,0 is the center of mass and scaling the points
    --     If a 2 points fit in the same place in the grid, only one of them is inserted.

    let pointsT = transpose points  
    let mins = map (reduce_comm f64.min f64.highest) pointsT
    let maxs = map (reduce_comm f64.max f64.lowest ) pointsT 
    let smallest_min = f64.minimum mins -- Pushes the points towards the first quodrant

    let x_scale = ((f64.i64 grid_size - 1) / (f64.max (f64.abs mins[0]) (f64.abs maxs[0])))
    let y_scale = ((f64.i64 grid_size - 1) / (f64.max (f64.abs mins[1]) (f64.abs maxs[1])))
    let scale   = f64.min x_scale y_scale 

    let (scaled_points, scaled_points_grid, scaled_points_flat_idx) = transpose pointsT
        |> map (\cc -> 
            let x =(cc[0] - smallest_min) * scale
            let y =(cc[1] - smallest_min) * scale
            let gp = [i64.f64 <| f64.round x, i64.f64 <| f64.round y] 
            in ([x, y], gp, gp[0] + gp[1] * grid_size) 
        ) |> (unzip3)     
        
    let (grid, used, unused) = populateGrid (grid_size * grid_size) scaled_points_flat_idx (indices points)
    let grid = unflatten grid
    in (grid, used, unused, scaled_points, scaled_points_grid)

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
                        let ind' = grid[r'][c']
                        in if ind' != -1 then
                            let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 c',f64.i64 r')
                            in (d', ind')
                        else
                            (f64.highest, -1)
                    ) stencil_1
                |> reduce (\acc x -> if x.0 < acc.0 then x else acc ) (f64.highest, -1)
        )

    -- in plus_1   
    in loop g = plus_1 for iteration < (log2Int grid_size) do
        let stencil_it = stencilK (grid_size / (2**(iteration+1)))
        in tabulate_2d grid_size grid_size 
            (\r c ->
                let (d, ind) = (g[r][c].0, g[r][c].1)
                in map (\(si, sj) ->
                    let (r', c') = (r + si, c + sj)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (d, ind)
                    else
                        let ind' = g[r'][c'].1
                        in if ind' != -1 then
                            let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 points[g[r'][c'].1][0], f64.i64 points[g[r'][c'].1][1])
                            in if d' < d then (d', ind') else (d, ind)
                        else
                            (d, ind)
                    ) stencil_it
                |> reduce (\acc x -> if x.0 < acc.0 then x else acc ) (g[r][c].0, g[r][c].1)
            )

def voronoiDiagram [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64)  : ([grid_size][grid_size](f64, i32)) =
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

    in loop g = plus_1 for iteration < (log2Int grid_size) do
        let stencil_it = stencilK (grid_size / (2**(iteration+1)))
        in tabulate_2d grid_size grid_size 
            (\r c ->
                loop (d, ind) =  (g[r][c].0, g[r][c].1) for (i,j) in stencil_it do
                    let (r', c') = (r + i, c + j)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (d, ind)
                    else 
                        let ind' = g[r'][c'].1
                        in if ind' != -1  then
                            let d' = dist (f64.i64 c,f64.i64 r) (f64.i64 points[g[r'][c'].1][0], f64.i64 points[g[r'][c'].1][1])
                            in if d' < d then
                                (d', ind')
                            else
                                (d, ind)      
                        else
                            (d, ind)
            )

def removeIslands [grid_size] [n] (grid : [grid_size][grid_size](f64, i32)) (points : [n][2]i64) : [grid_size][grid_size](f64, i32) =
    let stencil_1 = stencilK 1
    let (g'', _) = 
        loop (g, cond) = (grid, true) while cond do 
            let (g'', island_flag_array) =  unzip <| flatten <| tabulate_2d grid_size grid_size 
                (\r c ->
                    let (site_r, site_c) = (points[g[r][c].1][1], points[g[r][c].1][0])
                    let rule = findQuodrantOrAxis (r,c) (site_r, site_c) grid_size
                    in 
                        if      rule == #center then (g[r][c], false)
                        else if rule == #up    && g[r-1][c].1 == g[r][c].1                                            then (g[r][c], false)
                        else if rule == #down  && g[r+1][c].1 == g[r][c].1                                            then (g[r][c], false)
                        else if rule == #right && g[r][c+1].1 == g[r][c].1                                            then (g[r][c], false)
                        else if rule == #left  && g[r][c-1].1 == g[r][c].1                                            then (g[r][c], false)
                        else if rule == #leftUp    && (checkQuadrant g[r][c].1 g[r-1][c-1].1 g[r-1][c].1 g[r][c-1].1) then (g[r][c], false)
                        else if rule == #leftDown  && (checkQuadrant g[r][c].1 g[r+1][c-1].1 g[r+1][c].1 g[r][c-1].1) then (g[r][c], false)
                        else if rule == #rightUp   && (checkQuadrant g[r][c].1 g[r-1][c+1].1 g[r-1][c].1 g[r][c+1].1) then (g[r][c], false)
                        else if rule == #rightDown && (checkQuadrant g[r][c].1 g[r+1][c+1].1 g[r+1][c].1 g[r][c+1].1) then (g[r][c], false)

                        --Find new site to associate pixel with
                        else 
                            -- Find the closest point among the neighbours
                            -- let t33 = trace (rule, r,c, g[r][c])
                            loop ((d, ind), island_flag) =  ((f64.highest, -1), true) for (i,j) in stencil_1 do
                                let (r', c') = (r + i, c + j)
                                in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                                    ((d, ind), island_flag)
                                else 
                                    let ind' = g[r'][c'].1
                                    let (site_r, site_c) = (points[ind'][1], points[ind'][0])
                                    in if ind' != -1  then
                                        let d' = dist (f64.i64 r, f64.i64 c) (f64.i64 site_r, f64.i64 site_c)
                                        in if d' < d then
                                            ((d', ind'), island_flag)
                                        else
                                            ((d, ind), island_flag)
                                    else
                                        ((d, ind), island_flag)
                )
            -- Check if any islands where found. If so we loop again!
            -- let cond' = trace <| reduce (\acc f -> f && acc) false (island_flag_array)
            let cond' = reduce (\acc f -> f && acc) false (island_flag_array)
            in (unflatten g'', cond')
    in g''

def locateVoronoiVertices [m] (grid : [m][m]i32) : [m-2][m-2]i32 = 
    tabulate_2d (m -2) (m-2) 
        ( \r c ->
            let (r, c) = (r+1, c+1)
            let is_corner = classifyVertexI32 grid[r][c+1] grid[r][c] grid[r+1][c] grid[r+1][c+1]
            in is_corner
        )

def locateVoronoiVerticesAndCreateTriangulation [grid_size] (grid : [grid_size][grid_size]i32) : [grid_size][grid_size]((i32, i32, i32),(i32, i32, i32)) = --: [m](i64, i64) =
    -- tabulate_2d (grid_size - 2) (grid_size - 2) 
    let trs =tabulate_2d (grid_size) (grid_size) 
        ( \r c ->
            if r == 0 || c == 0 || r == (grid_size -1) || c == (grid_size - 1) then ((-1,-1,-1),(-1,-1,-1))
            else
                
                 classifyVertexAndTriangulation grid[r][c+1] grid[r][c] grid[r+1][c] grid[r+1][c+1]
        )

    -- let trs = filter  (\t -> t.0.0 > -1)  <| trs
    in trs

def createTriangles [m] (voronoi_diagram : [m][m]i32) (voronoi_vertices : [m-2][m-2]i32) : [](i32, i32, i32) =
    -- G5 and G6
    let fvv = flatten voronoi_vertices
    let fvv_ids = scan (+) 0i32 fvv
    let vv_idxs = map (\x -> if x.1 > 0 then (x.0, x.2) else (-1, 0) ) <| zip3 fvv_ids fvv (indices fvv)
    let vv_idxs' = filter (\x -> if x.0 < 0 then false else true) vv_idxs
    let num_ts   =  i64.i32 <| 2*(last fvv_ids)
    in  filter (\x -> x.0 >= 0) 
        <| map (\i -> 
            let j = vv_idxs'[i/2].1
            -- Since I only calculate the voronoi vertices on a (grid_size -2) (grid_size -2) grid I need to adjust the indices a bit
            let column_offset = j / (m - 2)
            let j = j + m + column_offset*2 + 1

            let c = j % m
            let r = j / m 
            let t = classifyVertexAndTriangulation voronoi_diagram[r][c+1] voronoi_diagram[r][c] voronoi_diagram[r+1][c] voronoi_diagram[r+1][c+1]
            in if i%2==0 then t.1 else t.0
        ) (iota num_ts)



def fixConvexHull [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64) : [](i32,i32,i32) =
    let edge = map (\i -> 
        if i < grid_size then grid[0][i]
        else if i < 2 * grid_size then grid[(i % grid_size)][grid_size -1]
        else if i < 3 * grid_size then grid[grid_size -1][grid_size - 1 - (i % grid_size)]
        else grid[grid_size - 1 - (i% grid_size)][0]

        ) (iota (4 * grid_size)) 
    let edge = pack_points_i32 edge 
    let triangles = replicate (length edge) (-1i32, -1i32, -1i32)

    let (_, triangles, len_triangles) = 
        loop (s, triangles, t_i) = ((0, 1, 2), triangles, 0i64) for i < ((length edge) - 2) do 
                if isClockwise points[edge[s.0]] points[edge[s.1]] points[edge[s.2]] grid_size then
                    ((s.1, s.2, s.2 + 1), triangles, t_i)
                else
                    ((s.0, s.2, s.2+1), triangles with [t_i] = (edge[s.0], edge[s.2], edge[s.1]), t_i + 1)
                
    in take len_triangles triangles
                


def shiftSites [t] [p] [n] (triangles : [t](i32, i32, i32)) (used_points : [p]i64) (points : [n][2]f64) (grid_points : [n][2]i64) =
    let test = map (\i -> [i.0, i.1, i.2]) triangles
    -- Seems like the maximum number of triangles in a fan is around 15 for the 1m dataset
    -- This means the loop that should be made might only need 15 iterations.
    let triangle_fans = map (\ts -> [(ts.0, ts), (ts.1, ts), (ts.2, ts)]) triangles
        |> flatten
        |> radix_sort_int_by_key (\(i, _) -> i) ((i32.i64 <| log2Int t) + 2) i32.get_bit

    let flag_arr = map2 (
        \(f, _) i ->  
            if i  > 0  then (f != triangle_fans[i-1].0) |> i32.bool 
            else 0i32 
        ) triangle_fans (indices triangle_fans)
    let f_idx = scan (+) 0 flag_arr

    -- Might be faster with a scatter instead of filter
    --  Then use f_idx with a -1
    let shp = sgmscan (+) (0) flag_arr (replicate (t * 3) 1) 
        |> map2 (\i sh -> 
            if i == (t * 3) - 1 then sh
            else
                if flag_arr[i+1] > 0 then sh
                else 0
                ) (indices flag_arr) 
        |> filter (>0)

    let is_shifted = (replicate (n) false ) --with [0] = true

    -- Compute valid sites
    -- Compute rules for valid sites
    -- Update fan ??
    -- call next iteration
    let (fan, shp, is, r') =
        loop (tfans, shp, is_shifted', r) = (triangle_fans, shp, is_shifted, [])
        while !(reduce (&&) true is_shifted') do

            let sc_shp = scan (+) 0 shp
            let sc_shp_ex = (rotate (-1) sc_shp) with [0] = 0 
            let sites = map (\i -> tfans[i].0) sc_shp_ex

            -- Do heuristic before I compute the rules with hist 
            -- Seems like there is some bug where the shifted '-1' are read from...
            let flat_triangles = map (\tf -> 
                    if is_shifted'[tf.0] then [-1, -1, -1] -- indices that are ignored
                    else map (i64.i32) [tf.1.0,tf.1.1,tf.1.2]
                ) tfans |> flatten

            let f_site_tri  = map (\tf -> replicate 3 tf.0) tfans 
                            |> flatten
            let H = hist i32.min (i32.i64 n) (n) flat_triangles f_site_tri
            let H = trace H

            let valid_sites = map2 (\i site ->
                if is_shifted'[i] then (false, site, i)
                else
                    let s = sc_shp_ex[i] * 3 -- *3 because indx for 
                    let e = sc_shp[i] * 3    --  flat_triangles
                    let indices = flat_triangles[s:e]
                    let test    = map (\ii -> H[ii] == site) indices
                            |> reduce (&&) true
                    in (test, site, i)
                ) (indices sc_shp_ex) sites 
                |> filter (.0)
                |> map (\vs -> (vs.1, vs.2))

            -- Find the rule for each non shifted site
            -- These map to the 3 scenarios from the paper:
            -- 1: in the original fan
            -- 2: In another fan
            -- 3: Outside mesh
            let rules = map (
                \(site, i) -> 
                    let s = sc_shp_ex[i]
                    let e = sc_shp[i]
                    let fan = map (.1) tfans[s:e]
                    let p   = points[site]
                    -- let is_in_fan = map (\tr -> 
                    --         let t1 = if is_shifted'[tr.0] then points[tr.0] else map (f64.i64) grid_points[tr.0]
                    --         let t2 = if is_shifted'[tr.1] then points[tr.1] else map (f64.i64) grid_points[tr.1]
                    --         let t3 = if is_shifted'[tr.2] then points[tr.2] else map (f64.i64) grid_points[tr.2]
                    --         in pointInTriangle t1 t2 t3 p
                    --     ) fan 
                    --     |> reduce (||) false
                    let in_fan_min_dist = map (\tr -> 
                            let t1 = if is_shifted'[tr.0] then points[tr.0] else map (f64.i64) grid_points[tr.0]
                            let t2 = if is_shifted'[tr.1] then points[tr.1] else map (f64.i64) grid_points[tr.1]
                            let t3 = if is_shifted'[tr.2] then points[tr.2] else map (f64.i64) grid_points[tr.2]
                            in pointInTriangleAndDist t1 t2 t3 p [tr.0,tr.1,tr.2]
                        ) fan 
                    let is_in_fan = map (.0) in_fan_min_dist |> reduce (||) false
                    in if is_in_fan then i64.bool is_in_fan
                    else 
                        let closest_point = map (.1) in_fan_min_dist 
                            |> reduce (\acc x -> if x.1 < acc.1 then x else acc) (-1i32, f64.highest)
                            |> (.0)
                        let j = binary_search closest_point sites
                        let s = sc_shp_ex[j]
                        let e = sc_shp[j]
                        let fan = map (.1) tfans[s:e]
                        let p   = points[closest_point]
                        let is_in_fan = map (\tr -> 
                                let t1 = if is_shifted'[tr.0] then points[tr.0] else map (f64.i64) grid_points[tr.0]
                                let t2 = if is_shifted'[tr.1] then points[tr.1] else map (f64.i64) grid_points[tr.1]
                                let t3 = if is_shifted'[tr.2] then points[tr.2] else map (f64.i64) grid_points[tr.2]
                                in pointInTriangle t1 t2 t3 p
                            ) fan 
                            |> reduce (||) false
                        in if is_in_fan then 2i64 else 3i64
                        
                    -- Check if inside triangle fan
                    -- If yes then rule 1
                    -- else check :
                    -- in i64.bool is_in_fan
                ) valid_sites

            let update_vals  = replicate (length rules) true
            let id_to_update = sized (length rules) (map (.1) valid_sites)
            let is_shifted''   = scatter is_shifted' (id_to_update) update_vals
            let rules = trace rules
            in (tfans, shp, is_shifted'', rules)


    in (shp, is, r')


-- ==
-- compiled random input {       [1000][2]f64 } 
-- compiled random input {    [1000000][2]f64 } 
-- compiled random input {   [10000000][2]f64 } 
def main [n]
    (points : [n][2]f64)  =
    let grid_size = 1024  --/ 32

    let (grid, used, unused, scaled_points, scaled_points_grid) = movePointsToGrid points grid_size

    let grid' = voronoiDiagram2  grid  scaled_points_grid
    let grid''  = removeIslands grid' scaled_points_grid
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    -- G5 and G6
    let triangles = createTriangles voronoi_diagram voronoi_vertices

    let triangles = (fixConvexHull voronoi_diagram scaled_points_grid) ++ triangles

    -- in length triangles 
    --in map (\i -> [triangles[i].0, triangles[i].1,triangles[i].2]) <| indices triangles
    in  (shiftSites triangles used scaled_points scaled_points_grid)
    -- in reduce (i64.max) 0 shps
    -- in triangleGrid grid voronoi_diagram triangles scaled_points_grid

-- > :img main ($loaddata "test_data10.txt")

-- Comments/ToDo
-- Tjek G3
-- Gør C1 hurtig


-- Mål:
-- Lav noget parallelt til C2, som muligvis kun håndtere 1 af deres 9 tilfælde.


-- Cs: Jeg burde bruge lave et 'greedy' prøv at løs så mange ting som muligt. Se om der er konflikter og prøv igen approach approach

-- C1: Man kan map |> removeDuplicates, for at få boundry (Skal være i original order.) 
--     Muligvis ikke muligt at tjekke parallelt. I figur 5 se mørkegrøn, grå, lysegrøn.

-- C2: Det ligner jeg for hvert site også burde holde styr på dens edges. Ellers skal jeg søge efter alle trekanter, der indeholder et site.
--     Muligvis kan dette ikke lade sig gøre, da triangulationen kan ændre sig. Det vil betyde at man skal recalculate, nogle af sitesnes 
--      edges. Dette betyder, at hvis man vil have sitesne og deres 'fan' i et flat array, skal man for hver ændring recalculate size arrayet.

-- For hver trekant lav (x, (x,y,z)), (y, (x,y,z)), (z, (x,y,z))
-- derefter sorter efter første coordinat.
-- Så kan man få en flad repræsentation. 

-- C3: I G1 burde jeg returnerer et par af de unused, så man kan referere det ubrugte site til det brugte. 
