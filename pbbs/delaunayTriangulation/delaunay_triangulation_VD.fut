--  futhark run delaunay_triangulation_VD.fut < test_data.txt
--  futhark dataset -g 5i32 -g [10][2]f32 > test_data.txt

import "util"
import "lib/github.com/diku-dk/sorts/radix_sort"



def movePointsToGrid [n] (points : [n][2]f32) (grid_size : i64) : ([grid_size][grid_size]i32, []i64, [n][2]i64) =

    -- G1: Moves points to a grid by translating the points such that 0,0 is the center of mass and scaling the points
    --     If a 2 points fit in the same place in the grid, only one of them is inserted.

    let pointsT = transpose points  
    let mins = map (reduce_comm f32.min f32.highest) pointsT
    let maxs = map (reduce_comm f32.max f32.lowest ) pointsT 
    let largest_min = f32.minimum mins -- Pushes the points towards the first quodrant

    let x_scale = ((f32.i64 grid_size - 1) / (f32.max (f32.abs mins[0]) (f32.abs maxs[0])))
    let y_scale = ((f32.i64 grid_size - 1) / (f32.max (f32.abs mins[1]) (f32.abs maxs[1])))
    let scale   = f32.min x_scale y_scale 

    let scaled_points = map (\cc -> map (\c -> i64.f32 <| f32.round <| (c - largest_min) * scale ) cc) <| transpose pointsT 

    let scaled_points_flat_idx =  map (\cc -> 
        let x = (i64.f32 <| f32.round <| (cc[0] - largest_min) * scale )   
        let y = grid_size * (i64.f32 <| f32.round <| (cc[1] - largest_min) * scale )
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
    -- in (grid, unused, scaled_points)


def voronoiDiagram [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64)  : ([grid_size][grid_size](f32, i32)) =
    -- G2
    -- I use a tuple grid to represent for each pixel both the distance to the closest encountered site, 
    --   this site's index and the coordinates to the original site which is now referenced
    --   by the specific pixel: (dist, site_index)
    let stencil_1 = stencilK 1
    let plus_1 = tabulate_2d grid_size grid_size 
        (\r c ->
            if grid[r][c] != -1 then
                (0f32, grid[r][c])
            else
                loop (d, ind) =  (f32.highest, -1) for (i,j) in stencil_1 do
                    let (r', c') = (r + i, c + j)
                    in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                        (d, ind)
                    else 
                        let d' = dist (f32.i64 c,f32.i64 r) (f32.i64 c',f32.i64 r')
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
                            let d' = dist (f32.i64 c,f32.i64 r) (f32.i64 points[g[r'][c'].1][0], f32.i64 points[g[r'][c'].1][1])
                            in if d' < d then
                                (d', ind')
                            else
                                (d, ind)      
                        else
                            (d, ind)
            )

def removeIslands [grid_size] [n] (grid : [grid_size][grid_size](f32, i32)) (points : [n][2]i64) : [grid_size][grid_size](f32, i32) =
    let stencil_1 = stencilK 1
    let (g'', _) = 
        loop (g, cond) = (grid, true) while cond do 
            let (g'', island_flag_array) =  unzip <| flatten <| tabulate_2d grid_size grid_size 
                (\r c ->
                    let (site_r, site_c) = (points[g[r][c].1][1], points[g[r][c].1][0])
                    let rule = findQuodrantOrAxis (r,c) (site_r, site_c) grid_size
                    let tmp = 
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
                            let t33 = trace (rule, r,c, g[r][c])
                            in loop ((d, ind), island_flag) =  ((f32.highest, -1), true) for (i,j) in stencil_1 do
                                let (r', c') = (r + i, c + j)
                                in if r' < 0 || r' >= grid_size || c' < 0 || c' >= grid_size then
                                    ((d, ind), island_flag)
                                else 
                                    let ind' = g[r'][c'].1
                                    let (site_r, site_c) = (points[ind'][1], points[ind'][0])
                                    in if ind' != -1  then
                                        let d' = dist (f32.i64 r, f32.i64 c) (f32.i64 site_r, f32.i64 site_c)
                                        in if d' < d then
                                            ((d', ind'), island_flag)
                                        else
                                            ((d, ind), island_flag)
                                    else
                                        ((d, ind), island_flag)
                    in tmp
                )
            -- Check if any islands where found. If so we loop again!
            let cond' = trace <| reduce (\acc f -> f || acc) false (island_flag_array)
            in (unflatten g'', cond')
    in g''

def locateVoronoiVertices [m] (grid : [m][m]i32) : [m-2][m-2]i32 = 
    -- tabulate_2d (grid_size ) (grid_size) 
    --     ( \r c ->
    --         if r == 0 || c == 0 || r == (grid_size -1) || c == (grid_size - 1) then 0
    --         else
    --             let is_corner = classifyVertexI32 grid[r][c+1] grid[r][c] grid[r+1][c] grid[r+1][c+1]
    --             in is_corner
    --     )
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



def fixConvexHull [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64) = -- [](i32,i32,i32)
    let edge = map (\i -> 
        if i < grid_size then grid[0][i]
        else if i < 2 * grid_size then grid[(i % grid_size)][grid_size -1]
        else if i < 3 * grid_size then grid[grid_size -1][grid_size - 1 - (i % grid_size)]
        else grid[grid_size - 1 - (i% grid_size)][0]

        ) (iota (4 * grid_size)) 
    -- let tc2 = trace grid
    let edge = pack_points_i32 edge -- It might be slower to remove duplicates

    let (_, triangles) = 
        loop (stack, triangles) = ([edge[0], edge[1]], [])
            for x in edge[2:] do 
                let stack = stack ++ [x]
                in if isClockwise points[stack[0]] points[stack[1]] points[stack[2]] grid_size then
                    (stack[1:], triangles)
                else
                    ([stack[0], stack[2]], triangles ++ [(stack[0], stack[2], stack[1])])
                
    in triangles
                


-- > :img main ($loaddata "test_data200.txt")

-- > :img main2 ($loaddata "test_data200.txt")

-- ==
-- compiled random input {       [1000][2]f32 } 
-- compiled random input {    [1000000][2]f32 } 
-- compiled random input {   [10000000][2]f32 } 
def main [n]
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

    --let triangles = (fixConvexHull voronoi_diagram scaled_points) ++ triangles

    in length triangles 
    --in map (\i -> [triangles[i].0, triangles[i].1,triangles[i].2]) <| indices triangles
    -- in length triangles
    -- in triangleGrid grid voronoi_diagram triangles scaled_points

def main2 [n]
    (points : [n][2]f32)  =
    let grid_size = 1024

    let (grid, unused, scaled_points) = movePointsToGrid points grid_size

    let grid' =  voronoiDiagram grid scaled_points
    let grid' =  removeIslands grid' scaled_points
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid'[i][j].1) 
    let voronoi_vertices = locateVoronoiVertices voronoi_diagram 

    -- G5 and G6
    let triangles = createTriangles voronoi_diagram voronoi_vertices

    let triangles = (fixConvexHull voronoi_diagram scaled_points) ++ triangles

    in triangleGrid grid voronoi_diagram triangles scaled_points

-- Comments/ToDo

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
